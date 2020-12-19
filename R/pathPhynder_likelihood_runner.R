suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(phytools)))
suppressWarnings(suppressPackageStartupMessages(library(phangorn)))
suppressWarnings(suppressPackageStartupMessages(library(pracma)))

tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)
packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')

cat('\n\n',"pathPhynder_likelihood_runner.R", '\n\n\n')

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=5) {
    stop("  Arguments needed.\n
        \tusage
        \tRscript pathPhynder_likelihood_runner.R <input_phylogeny.nwk> <input.vcf> status.txt <results_folder> <prior> <out_prefix>
        ", call.=FALSE)
}



source(paste0(packpwd,"/pathPhynder_likelihood_functions.R"))

tree_file=args[1]
vcf_name<-args[2]
results_folder=args[3]
prior=args[4]
out_prefix=args[5]
calls_file<-paste0("intree_folder/",out_prefix, '.status.txt')


for (testfile in c(tree_file, vcf_name, calls_file)){
    if (!file.exists(testfile)) {
        stop(paste(testfile, "- does this file exist?"))
    }
}

dir.create(results_folder, showWarnings = FALSE)

#read vcf
vcf<-fread(vcf_name, skip = "CHROM", stringsAsFactors=F, na.strings = '.')
vcf$REF[vcf$REF==TRUE]<-"T"
vcf$ALT[vcf$ALT==TRUE]<-"T"
vcf<-data.frame(vcf)

#exclude multiallelic
old_dim_vcf<-dim(vcf)[1]
vcf<-vcf[grep(',',vcf$REF, invert=T),]
vcf<-vcf[grep(',',vcf$ALT, invert=T),]
new_dim_vcf<-dim(vcf)[1]
cat(paste0('excluded ',old_dim_vcf-new_dim_vcf, ' multiallelic sites from analysis.'), '\n')

old_POS<-vcf$POS
vcf<-vcf[nchar(vcf$REF)==1,]
vcf<-vcf[nchar(vcf$ALT)==1,]
new_POS<-vcf$POS


excluded_nonSingleNuc<-old_POS[!old_POS %in% new_POS]
cat(paste0('excluded ',length(old_POS)-length(new_POS), ' sites with base length greater than 1 from analysis.'))
cat('\n') 

cat(paste0("Prior = ",prior), '\n')


#get previously called query samples
getQueries<-function(sample_name){
	# sample_name<-gsub('.*/','',sample_path)
	query<-read.table(paste0('intree_folder/',sample_name,'.status.txt'))
	colnames(query)[1]<-'POS'
	query_name<-paste0('QUERY_', sample_name)
	colnames(query)[6]<-query_name
	query<-query[c('POS', query_name)]
	query<-query[query[[query_name]]!="-9",]
	return(query)
}

finalGenoMatrix<-vcf
query<-NULL
query<-getQueries(out_prefix)
# length(vcf$POS)
#this leaves missing genos in query as NA
finalGenoMatrix<-merge(finalGenoMatrix, query, all=T)
cat(paste(100-round(dim(query)[1]*100/length(vcf$POS),2), '%', 'missingness.'), '\n')

#remove unnecessary info
finalGenoMatrix<-finalGenoMatrix[which(!colnames(finalGenoMatrix) %in% c("X.CHROM","ID","QUAL","FILTER","INFO","FORMAT"))]
finalGenoMatrix[is.na(finalGenoMatrix)]<-'.'

# setwd("Documents/academic/phd/year1/pathphynder/") 
library(phytools); library(phangorn); library(pracma)
tree <- ladderize(read.tree(file=tree_file))
# tree <- (read.tree(file=tree_file))





# data <- read.table("finalGenoMatrix.txt",header=TRUE,stringsAsFactors=FALSE,colClasses=c("character"))
tree$tip.label = make.names(tree$tip.label)

query.name<-colnames(finalGenoMatrix)[grep('QUERY',colnames(finalGenoMatrix))]

data<-finalGenoMatrix
colnames(data)<-make.names(colnames(data))


tip.nucleotides = data[!(colnames(data) %in% c('POS','REF','ALT'))]
prior.type = prior 

# general set up assuming the format above
tree$edge.length = abs(tree$edge.length)


if (length(tree$edge.length[tree$edge.length==0])){
  tree$edge.length[tree$edge.length==0]<-min(tree$edge.length[tree$edge.length>0])
  #if branch length 0, it will crash!
  print('IMPORTANT - setting branches with length 0 to minimum length.')
}


rows.with.missing.data = which(apply(tip.nucleotides, 1, function(r) any(r == ".")))
rows.no.snps = which(apply(tip.nucleotides, 1, function(r) all(r == "0") | all(r=="1")))
rows.to.go = union(rows.with.missing.data,rows.no.snps)
if(length(rows.to.go)>0){tip.nucleotides = tip.nucleotides[-rows.to.go,]}
if(length(rows.to.go)>0){
  refs <- data$REF[-rows.to.go]
} else{
  refs <- data$REF
}
if(length(rows.to.go)>0){
  alts = data$ALT[-rows.to.go]
}else{
  alts = data$ALT
}

temp<-NULL

temp = cbind(tip.nucleotides,refs,alts)
num.samples = ncol(tip.nucleotides); num.sites = nrow(tip.nucleotides); num.edge = nrow(tree$edge)
# set the evolutionary rate model for pml (see description of pml function)
model = "GTR" # make this the default, and make it a user option i think

# likelihood calculation on each edge
fake.fasta = apply(temp,1, function(x) ifelse( x[1:num.samples]==0,x[num.samples+1],x[num.samples+2]) )
phydat = phyDat(fake.fasta)
edge.loglik = function(edge.num){
  tree.with.query = tree.constructor(tree,edge.num,new.edge.length=1e-3*tree$edge.length[edge.num],pos=0.5*tree$edge.length[edge.num])
  return(pml(tree.with.query,phydat,model=model)$logLik)
}

cat(paste0("Sample name: ",query.name),'\n')

cat(paste0(dim(temp)[1], ' SNPs wih no missingness in the data.'), '\n')

begin_time<-Sys.time()
logliks = rep(0,num.edge)
for(i in 1:num.edge){
  logliks[i] = edge.loglik(i)
  visited_edges_message<-paste0('Edge ',i,'/', num.edge,'.',' loglik = ', logliks[i], '.')
  cat("\r",visited_edges_message, ' time',  begin_time-Sys.time())
}
cat("\n")

# calculation of posterior probabilities
prior.dist = calculate.prior(tree,prior.type)
named.logliks = logliks; names(named.logliks) = 1:num.edge
temppost = exp(rev(sort(named.logliks)) - max(named.logliks))
posterior  = temppost*prior.dist[as.numeric(names(temppost))] / sum(temppost * prior.dist[as.numeric(names(temppost))] )
# ^ this is them in order of descending prob. v this is them in order of edge number
posterior.normal.order = unname(posterior[as.character(1:num.edge)])

#choose the edge with highest posterior probability
final.edge.num = as.numeric(names(which(posterior==max(posterior))))


if (length(final.edge.num)>1) {
  result<-data.frame(final_edge=final.edge.num,position=rep(NA, length(final.edge.num)),branch.length=rep(NA, length(final.edge.num)),pprob=max(posterior))
    cat(paste0("We've got ",length(final.edge.num)," multiple equally likely placements (posterior probability=", round(max(posterior),4),"), not estimating position along these edges."), '\n')
    write.table(result,file=paste0(results_folder,'/',out_prefix,".results.",prior,".txt"), row.names=F, quote=F)
    write.table(posterior.normal.order, file=paste0(results_folder,'/',out_prefix,".posteriors.",prior,".txt"), quote=F)
} else {
  branch.length = tree$edge.length[final.edge.num]
  ymin = 1e-3 * branch.length # minimum place on edge
  ymax = (1-1e-3) * branch.length # maximum place on edge
  new.edge.length = 1e-3 * min(tree$edge.length) # length of query branch, effectively 0
  fels = function(y) -pml(tree.constructor( tree,final.edge.num,new.edge.length,pos=y),phydat,model=model)$logLik
  robot = optim(par=branch.length/2,fn=fels,method="L-BFGS-B",
                lower = ymin,upper=ymax,control = list(parscale=ymin))
  new.tree = tree.constructor(tree,final.edge.num,new.edge.length,pos=robot$par)
  newqueryedge = which(new.tree$edge[,2] == which(new.tree$tip.label == query.name))

  # if one result
  result<-data.frame(final_edge=final.edge.num, position=robot$par, branch.length=branch.length, pprob=posterior[1])

  write.table(result,file=paste0(results_folder,'/',out_prefix,".results.",prior,".txt"), row.names=F, quote=F)
  write.table(posterior.normal.order, file=paste0(results_folder,'/',out_prefix,".posteriors.",prior,".txt"), quote=F)
  write.tree(new.tree,file=paste0(results_folder,'/',out_prefix,".new_tree.",prior,".txt"))

  plot_final_tree(new.tree, out_prefix, results_folder, prior)

  cat(paste0("Best edge: ",final.edge.num,". Posterior probability=", round(max(posterior),4)), '\n')
}




#plotting
plot_likes(tree, posterior.normal.order, out_prefix,results_folder, prior)

