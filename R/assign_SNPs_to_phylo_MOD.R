# usage:
# Rscript assign_SNPs_to_phylo.R <input_phylogeny.nwk> <input.vcf> <out prefix>

require(phytools, quietly = TRUE)
#require(phangorn, quietly = TRUE)
suppressPackageStartupMessages(require(svMisc))


cat('\n\n',"Assigning SNPs to branches", '\n\n\n')



args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("\tArguments needed.\n
\tusage:
\tRscript assign_SNPs_to_phylo.R <input_phylogeny.nwk> <input.vcf> <out prefix>", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("assign_SNPs_to_phylo.R", args[1], args[2],args[3]), '\n\n')
}

#read tree
tree<-read.tree(file=args[1])
tree<-ladderize(tree)

#fix names prevent downstream problems withs sample ID
tree$tip.label<-make.names(tree$tip.label)


#read vcf
get_vcf<-function(vcf_name){
	all_content = readLines(vcf_name)
	skip = all_content[-c(grep("CHROM",all_content))]
	vcf <- read.table(textConnection(skip), stringsAsFactors=F)
	header<-unlist(strsplit(all_content[grep("CHROM", all_content)[length(grep("CHROM", all_content))]], '\t'))
	colnames(vcf)<-make.names(header)
	all_content<-NULL
	return(vcf)
}

vcf<-get_vcf(args[2])


#check if any samples are missing and if so, exclude from vcf
miss<- colnames(vcf)[10:length(vcf)][!(colnames(vcf)[10:length(vcf)] %in% tree$tip.label)]

if(length(miss)>0){
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-10),'\n')
    cat(paste0("    Excluded samples in VCF not present in the tree: ", miss),'\n')
    vcf<-vcf[!(colnames(vcf) %in% miss)]
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-10),'\n')
}

#remove if unnecessary
#check if any samples are missing and if so, exclude from vcf
# vcf<-vcf[colnames(genos) %in% tree$tip.label]


if (length(tree$tip.label[!(tree$tip.label %in% colnames(vcf))]) >0){
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-10),'\n')
    cat(paste0("    Excluded ", length(tree$tip.label[!(tree$tip.label %in% colnames(vcf))]), " individuals from the tree not in VCF. (See tree_data/", args[3],".indsremoved.txt)"), '\n')
    write.table(file=paste0("tree_data/", args[3],".indsremoved.txt"),data.frame(tree$tip.label[!(tree$tip.label %in% colnames(vcf[10:dim(vcf)[2]]))]), quote = F, row.names = F, col.names = F, sep='\t')
    cat(paste0("    WARNING: Wrote new tree to tree_data/tree.", args[3],".indsremoved.txt"), '\n')
    cat('\n')
    tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(vcf))])
    write.tree(file=paste0("tree_data/tree.", args[3],".indsremoved.txt"),tree)

}


# > b<-drop.tip(b, b$tip.label[!(b$tip.label %in% colnames(vcf[10:dim(vcf)[2]]))])


#Remove missing data
#it would be good to fix this to impute missingness in a phylogenetically aware way
#I have fixed it, missing data is now being accounted for and used
samples<-colnames(vcf[10:dim(vcf)[2]])

vcf[samples]<-sapply(vcf[samples], as.numeric)


complete_vcf<-vcf
complete_vcf[samples][complete_vcf[samples]=='N']<-NA
complete_vcf<-complete_vcf[complete.cases(complete_vcf[samples]),]

vcf_with_missing<-vcf[!vcf$POS %in% complete_vcf$POS,]



cat(paste0("    Number of SNPs with any missing data: ", dim(vcf)[1]),'\n')
cat(paste0("    Number of SNPs with no missing data: ", dim(complete_vcf)[1]),'\n')




edges<-data.frame(tree$edge)
colnames(edges)<-c('pos1','pos2')
edges$edge<-rownames(edges)
snps<-read.table("snps.txt")


cat("starting non missing genos\n\n")

REFpos<-list()
ALTpos<-list()
total <- length(edges$edge)
# pb <- txtProgressBar(min = 0, max = total, style = 3)

for (edge in edges$edge){
    relevant_node<-edges$pos2[edges$edge==edge]
    desc<-tree$tip.label[getDescendants(tree, relevant_node)][!is.na(tree$tip.label[getDescendants(tree,relevant_node)])]
    nondesc<- samples[!samples %in% desc]
    REFs<-which(rowSums(complete_vcf[desc])==0 & rowSums(complete_vcf[nondesc])==length(nondesc))
    ALTs<-which(rowSums(complete_vcf[desc])==length(desc) & rowSums(complete_vcf[nondesc])==0)
    REFpos[[edge]]<-unique(complete_vcf$POS[REFs])
    ALTpos[[edge]]<-unique(complete_vcf$POS[ALTs])
    allele_count<-(paste0(edge,'/', length(edges$edge),' nodes;    found ' ,"REFs=", length(REFs),' / ', "ALTs=", length(ALTs)))
    cat("\r",allele_count)
# 	snp_info<-sort(as.character(unique(snps$V2[snps$V3 %in% c( unlist(REFpos[[edge]]),unlist(ALTpos[[edge]]) ) ])))
# 	if (length(snp_info)>0){
# 		print(snp_info)
# 	}
# }
}
cat('\n\n\n')



# cat("starting non missing genos\n\n")

# #iterates through tree from the root to tips
# #gets descendants of a given node
# #tests this set of descendants all have 1s and all remaining have 0s (der)
# #tests this set of descendants all have 0s and all remaining have 1s (der)
# #records the positions at each branch which are composed by the ALT (derpos) and REF (ancpos) allele



for (edge in edges$edge){
    relevant_node<-edges$pos2[edges$edge==edge]
    desc<-tree$tip.label[getDescendants(tree, relevant_node)][!is.na(tree$tip.label[getDescendants(tree,relevant_node)])]
    nondesc<- samples[!samples %in% desc]
    vcf_with_missing$na_count_samples <- apply(vcf_with_missing[samples], 1, function(x) sum(is.na(x)))
    vcf_with_missing$na_count_desc <- apply(vcf_with_missing[desc], 1, function(x) sum(is.na(x)))
    vcf_with_missing$na_count_nondesc <- apply(vcf_with_missing[nondesc], 1, function(x) sum(is.na(x)))
    REFs<-which(rowSums(vcf_with_missing[desc], na.rm=T)==0 & rowSums(vcf_with_missing[nondesc], na.rm=T)==length(nondesc)+vcf_with_missing$na_count_nondesc)
    ALTs<-which(rowSums(vcf_with_missing[desc], na.rm=T)==length(desc)+vcf_with_missing$na_count_desc & rowSums(vcf_with_missing[nondesc], na.rm=T)==0)
    REFpos[[edge]]<-c(unlist(REFpos[[edge]]),unique(vcf_with_missing$POS[REFs]))
    ALTpos[[edge]]<-c(unlist(ALTpos[[edge]]),unique(vcf_with_missing$POS[ALTs]))
	allele_count<-(paste0(edge,'/', length(edges$edge),' nodes;    found ' ,"REFs=", length(REFs),' / ', "ALTs=", length(ALTs)))
    cat("\r",allele_count)
	# snp_info<-sort(as.character(unique(snps$V2[snps$V3 %in% c( unlist(REFpos[[edge]]),unlist(ALTpos[[edge]]) ) ])))
	# if (length(snp_info)>0){
	# 	print(snp_info)
	# }
}





dir.create('tree_data', showWarnings = FALSE)



saveRDS(ALTpos, file=paste0('tree_data/',args[3],".derpos.RData"))
saveRDS(REFpos, file=paste0('tree_data/',args[3],".ancpos.RData"))

write.table(unlist(ALTpos), file=paste0('tree_data/',args[3],".derpos.txt"), quote = F, row.names = F, col.names = F, sep='\t')
write.table(unlist(REFpos), file=paste0('tree_data/',args[3],".ancpos.txt"), quote = F, row.names = F, col.names = F, sep='\t')


#make positions
pos_to_call<-as.data.frame(unique(sort(c(unlist(ALTpos), unlist(REFpos)))))
colnames(pos_to_call)<-"pos_to_call"


pos_to_call$chr <- 'Y'
pos_to_call$pos0<-pos_to_call$pos_to_call-1
pos_to_call$REF<-vcf$REF[match(pos_to_call$pos_to_call, vcf$POS)]
pos_to_call$ALT<-vcf$ALT[match(pos_to_call$pos_to_call, vcf$POS)]
pos_to_call$mut<-paste0(pos_to_call$REF,'->', pos_to_call$ALT)
pos_to_call$marker<-NA
pos_to_call$hg<-NA

pos_to_call<-pos_to_call[c('chr', "marker","hg", 'pos_to_call', 'mut', 'REF', 'ALT')]
# chrY    A2607   A00     7029908 C->T    C       T


bedfile_data<-data.frame(pos_to_call$chr, pos_to_call$pos_to_call-1,pos_to_call$pos_to_call)

bedfile_data_w_chr<-bedfile_data
colnames(bedfile_data_w_chr)[1]<-"chr"
bedfile_data_w_chr$chr<-paste0("chr",bedfile_data_w_chr$chr)


#write these in pos to call format and in bed format
write.table(file=paste0('tree_data/',args[3],'.sites.txt'),pos_to_call, quote = F, row.names = F, col.names = F, sep='\t')
write.table(file=paste0('tree_data/',args[3],'.sites.bed'),bedfile_data, quote = F, row.names = F, col.names = F, sep='\t')
write.table(file=paste0('tree_data/',args[3],'.siteschr.bed'),bedfile_data_w_chr, quote = F, row.names = F, col.names = F, sep='\t')

cat(paste0("\t",dim(unique(pos_to_call))[1]," informative positions for variant calling (written to tree_data/", args[3],".sites.bed)"),'\n\n')
cat(paste0("\t",dim(unique(pos_to_call))[1]," informative positions for filtering step (written to tree_data/", args[3],".sites.txt)"),'\n\n')

not_added<-unique(vcf$POS[!vcf$POS %in% pos_to_call$pos_to_call])

write.table(file=paste0('tree_data/',args[3],'.not_added.txt'),not_added, quote = F, row.names = F, col.names = F, sep='\t')
cat(paste0("\t",length(not_added)," positions were not added (written to tree_data/", args[3],".not_added.txt)"),'\n\n')

#Rscript assign_SNPs_to_phylo.R /Users/rm890/Integrating_Y/karmin_data/RAxML_bestTree.karmin.aln.mGTR.50boot_try2 /Users/rm890/Botai_project/Karmin/try_karmin/karmin_w_changes.vcf 50boot_try2




# lens<-NULL
# for (i in a$Edge){
#     tmp<-a[a$Edge==i,]
#     lens[i]<-(length(unique(unlist(strsplit(tmp$positions[1], '\\;')))))
# }
# write.table(a,"numbers.txt")

#read numbers of snps
nums<-read.table("numbers.txt", h=T, stringsAsFactors=F)



make_edge_df<-function(der, anc){
	all_list<-list()
	for (i in 1:length(der)){
		all_list[[i]]<-c(der[[i]], anc[[i]])
	}
    position_counts<-sapply(all_list, length)

	edge_df<-data.frame(tree$edge)
	colnames(edge_df)<-c("Node1","Node2")
	edge_df$Edge<-rownames(edge_df)


	tmp_desc<-NULL
	tmp_pos<-NULL
	known_hg<-NULL
	known_markers<-NULL
	for (i in 1:length(all_list)){
		tmp_pos[i]<-(paste(unique(all_list[[i]]), collapse=";"))
		known_hg[i]<-paste(sort(unique(snps[match(all_list[[i]], snps$V3)[!is.na(match(all_list[[i]], snps$V3))],]$V2)), collapse=';')
		known_markers[i]<-paste(sort(unique(snps[match(all_list[[i]], snps$V3)[!is.na(match(all_list[[i]], snps$V3))],]$V1)), collapse=',')
		tmp_desc[i]<-paste(unique(tree$tip.label[getDescendants(tree,edge_df[edge_df$Edge==i,]$Node2)][!is.na(tree$tip.label[getDescendants(tree,edge_df[edge_df$Edge==i,]$Node2)])]), collapse=';')
        snp_count[i]<-position_counts[i]
	}


	edge_df$positions<-tmp_pos
	edge_df$hg<-known_hg
	edge_df$markers<-known_markers
	edge_df$descendants<-tmp_desc

        edge_df<-edge_df[c('Edge','Node1','Node2','positions','hg','markers','descendants')]
	return(edge_df)

}

edge_df<-make_edge_df(ALTpos, REFpos)

write.table(edge_df,file=paste0("tree_data/",args[3],".edge_df.txt"), quote=F, row.names=F, sep='\t')

cat(paste0("\t","list of markers assigned to each node writtent to written to tree_data/", args[3],".edge_df.txt)"),'\n\n')
