# pathPhynder
# Author: Rui Martiniano
# Contact: rm890 [at] cam.ac.uk
# usage: Rscript assign_SNPs_to_phylo.R <input_phylogeny.nwk> <input.vcf> <out prefix>

require(data.table, quietly = TRUE)
require(phytools, quietly = TRUE)

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


#####################################
#           FUNCTIONS
#####################################


# #' Read in VCF file
# #' 
# #' @param vcf_name
# #' @return a vcf table
# #' @examples
# #' get_vcf('example.vcf')
# get_vcf <- function(vcf_name){
#     all_content = readLines(vcf_name)
#     skip = all_content[-c(grep("CHROM",all_content))]
#     vcf <- read.table(textConnection(skip), stringsAsFactors=F)
#     header <- unlist(strsplit(all_content[grep("CHROM", all_content)[length(grep("CHROM", all_content))]], '\t'))
#     colnames(vcf) <- make.names(header)
#     #if T alleles read as TRUE, convert to character T.
#     vcf$REF[vcf$REF==TRUE]<-"T"
#     vcf$ALT[vcf$ALT==TRUE]<-"T"
#     all_content<-NULL
#     return(vcf)
# }







#' Makes a dataframe with SNPs assigned to branches of the tree
#' 
#' @param der list with ALT alleles
#' @param anc list with REF alleles
#' @return Edges dataframe
#' @examples
#' make_edge_df(ALTlist, REFlist)
make_edge_df <- function(der, anc){
    snp_count <- NULL
    tmp_desc <- NULL
    tmp_pos <- NULL
    known_hg <- NULL
    known_markers <- NULL

    all_list <- list()
    for (i in 1:length(der)){
        all_list[[i]] <- c(der[[i]], anc[[i]])
    }

    position_counts <- sapply(all_list, length)

    edge_df <- data.frame(tree$edge)
    colnames(edge_df) <- c("Node1","Node2")
    edge_df$Edge <- rownames(edge_df)

    for (i in 1:length(all_list)){
        tmp_pos[i] <- (paste(unique(all_list[[i]]), collapse=";"))
        known_hg[i] <- paste(sort(unique(snps[match(all_list[[i]], snps$V3)[!is.na(match(all_list[[i]], snps$V3))],]$V2)), collapse=';')
        known_markers[i] <- paste(sort(unique(snps[match(all_list[[i]], snps$V3)[!is.na(match(all_list[[i]], snps$V3))],]$V1)), collapse=',')
        tmp_desc[i] <- paste(unique(tree$tip.label[getDescendants(tree,edge_df[edge_df$Edge==i,]$Node2)][!is.na(tree$tip.label[getDescendants(tree,edge_df[edge_df$Edge==i,]$Node2)])]), collapse=';')
        snp_count[i] <- position_counts[i]
    }

    edge_df$positions <- tmp_pos
    edge_df$hg <- known_hg
    edge_df$markers <- known_markers
    edge_df$descendants <- tmp_desc
    edge_df$snp_count <- snp_count

    edge_df <- edge_df[c('Edge','Node1','Node2','positions','hg','markers','descendants','snp_count')]
    return(edge_df)

}


#' Makes a long format dataframe with information about each SNP assigned
#' 
#' @param der list with ALT alleles
#' @param anc list with REF alleles
#' @return Long SNP table
#' @examples
#' makeLongSNPtable(der, anc)
makeLongSNPtable <- function(der, anc){
    snp_tab<-data.frame(matrix(ncol=5, nrow=0))
    colnames(snp_tab) <- c('Edge','position','marker','hg','status')

    for (i in 1:length(der)){
        if (length(unique(der[[i]]))>0){
            for (position in unique(der[[i]])){
                Edge <- as.character(i)
                status <- '+'
                snp_tab <-  rbind(snp_tab,data.frame(Edge,position, status))
            }
        }
        if (length(unique(anc[[i]]))>0){           
            for (position in unique(anc[[i]])){
                Edge <- as.character(i)
                status <- '-'
                snp_tab <-  rbind(snp_tab,data.frame(Edge,position, status))
            }
        }
    }

    snp_tab$hg <- as.character(snps$V2[match(snp_tab$position, snps$V3)])
    snp_tab$marker <- as.character(snps$V1[match(snp_tab$position, snps$V3)])

    snp_tab <- (unique(snp_tab))
    snp_tab$REF <- vcf$REF[match(snp_tab$position, vcf$POS)]
    snp_tab$ALT <- vcf$ALT[match(snp_tab$position, vcf$POS)]
    snp_tab$mutation_in_vcf <- paste0(snp_tab$REF,'->', snp_tab$ALT)
    snp_tab$chr <- contig

    snp_tab <- snp_tab[c('chr', 'marker','hg', 'position', 'mutation_in_vcf', 'REF', 'ALT', 'Edge','status')]
    

    write.table(file=paste0('tree_data/',args[3],'.sites.txt'),snp_tab, quote = F, row.names = F, col.names = F, sep='\t')
    
    cat(paste0("\t",length(unique(snp_tab$position))," informative positions (auxiliary information written to tree_data/", args[3],".sites.txt)"),'\n\n')

    return(snp_tab)
}


#' Makes bed files for SNP calling
#' 
#' @return Long SNP table
#' @examples
#' writeBed()
writeBed <- function(){
    if (dim(LongSNPtable)[1]==0){
        stop('\n\n', '\tNo positions in the VCF were assigned. Confirm that your vcf and tree match the requirements: - vcf needs to be haploid and biallelic','\n\n')
    } else {
        bedfile_data<-data.frame(LongSNPtable$chr, LongSNPtable$position-1,LongSNPtable$position)
        colnames(bedfile_data)<-c("chr","pos1", "pos2")

        bedfile_data<-unique(bedfile_data[order(bedfile_data$pos1),])

        bedfile_data_w_chr<-bedfile_data
        bedfile_data_w_chr$chr<-paste0("chr",bedfile_data_w_chr$chr)
    
        write.table(file=paste0('tree_data/',args[3],'.sites.bed'),bedfile_data, quote = F, row.names = F, col.names = F, sep='\t')
        write.table(file=paste0('tree_data/',args[3],'.siteschr.bed'),bedfile_data_w_chr, quote = F, row.names = F, col.names = F, sep='\t')

        cat(paste0("\t",dim(unique(bedfile_data))[1]," informative positions for variant calling (written to tree_data/", args[3],".sites.bed)"),'\n\n')

    }
}


#' Makes a report of SNPs which were not assigned to branches
#' 
#' @return Long SNP table
#' @examples
#' makeReport()
makeReport<-function(){
    # number of snps not added and why
    # Reasons:
    # 1) observed alleles are incompatible with tree topology
    # 2) missing data
    # 3) SNPs are monomorphic

    snps_not_added<-unique(vcf$POS[!vcf$POS %in% LongSNPtable$position])
    snps_not_added<-unique(c(snps_not_added,excluded_multiallelic))
    if (length(snps_not_added)>0){

        monomorphic_snps_ALT<-vcf$POS[which(rowSums(vcf[10:dim(vcf)[2]])==length(colnames(vcf)[10:dim(vcf)[2]]))]
        monomorphic_snps_REF<-vcf$POS[which(rowSums(vcf[10:dim(vcf)[2]])==0)]
        monomorphic_snps <- c(monomorphic_snps_ALT, monomorphic_snps_REF)


        report_df<-data.frame(position=snps_not_added)
        report_df$reason<-'incompatible with tree'

        if (length(vcf_with_missing$POS)>0){
            report_df$reason[report_df$position %in% vcf_with_missing$POS]<-'missing data'
        }

        if (length(monomorphic_snps)>0){
            report_df$reason[report_df$position %in% monomorphic_snps]<-'monomorphic snp'
        }

        write.table(file=paste0('tree_data/',args[3],'.report_not_added_SNPs.txt'),report_df, quote = F, row.names = F, col.names = F, sep='\t')
        cat(paste0("\t",dim(report_df)[1]," positions were not added (report written to tree_data/", args[3],".report_not_added_SNPs.txt)"),'\n\n')
    }
}




#####################################
#              MAIN
#####################################


dir.create('tree_data', showWarnings = FALSE)

#read tree
tree<-read.tree(file=args[1])
tree<-ladderize(tree)

#fix names prevent downstream problems withs sample ID
tree$tip.label<-make.names(tree$tip.label)

vcf_name<-args[2]

#Read in VCF file
vcf<-fread(vcf_name, skip = "CHROM", stringsAsFactors=F, na.strings = '.')
vcf<-data.frame(vcf)

#get contig
contig<-unique(vcf$X.CHROM)

#if more than one contig, error
if (length(unique(vcf$X.CHROM))>1){
    stop('Please provide a vcf with a single chromosome/contig of interest.')
}




#annoying thing where R reads Ts as bool
vcf$REF[vcf$REF=="TRUE"]<-'T'
vcf$ALT[vcf$ALT=="TRUE"]<-'T'

#exclude multiallelic

old_POS<-vcf$POS
old_dim_vcf<-dim(vcf)[1]
vcf<-vcf[grep(',',vcf$REF, invert=T),]
vcf<-vcf[grep(',',vcf$ALT, invert=T),]
new_POS<-vcf$POS

excluded_multiallelic<-old_POS[!old_POS %in% new_POS]
cat(paste0('excluded ',length(old_POS)-length(new_POS), ' multiallelic sites from analysis.'))
cat('\n') 

old_POS<-vcf$POS
vcf<-vcf[nchar(vcf$REF)==1,]
vcf<-vcf[nchar(vcf$ALT)==1,]
new_POS<-vcf$POS


excluded_nonSingleNuc<-old_POS[!old_POS %in% new_POS]
cat(paste0('excluded ',length(old_POS)-length(new_POS), ' sites with base length greater than 1 from analysis.'))
cat('\n') 



samples<-colnames(vcf[10:dim(vcf)[2]])

#check if any samples are missing and if so, exclude from vcf
miss<- colnames(vcf)[10:length(vcf)][!(colnames(vcf)[10:length(vcf)] %in% tree$tip.label)]

if(length(miss)>0){
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-9),'\n')
    cat(paste0("    Excluded samples in VCF not present in the tree: ",paste(miss, collapse=';')),'\n')
    if (length(miss)>=length(samples)-2){
        stop("There are no individuals left for analysis. Check the sample IDs in your VCF file and tree.\n")
    } else {
        vcf<-vcf[!(colnames(vcf) %in% miss)]
        cat(paste0("    Number of individuals: ", dim(vcf)[2]-9),'\n')
        samples<-colnames(vcf[10:dim(vcf)[2]])

}}

if (length(tree$tip.label[!(tree$tip.label %in% colnames(vcf))]) >0){
    cat(paste0("    Number of individuals: ", dim(vcf)[2]-9),'\n')
    cat(paste0("    Excluded ", length(tree$tip.label[!(tree$tip.label %in% colnames(vcf))]), " individuals from the tree not in VCF. (See tree_data/", args[3],".indsremoved.txt)"), '\n')
    write.table(file=paste0("tree_data/", args[3],".indsremoved.txt"),data.frame(tree$tip.label[!(tree$tip.label %in% colnames(vcf[10:dim(vcf)[2]]))]), quote = F, row.names = F, col.names = F, sep='\t')
    cat(paste0("    WARNING: Wrote new tree to tree_data/tree.", args[3],".indsremoved.nwk"), '\n')
    cat(paste0("    Use this new tree for further analyses or fix the input tree."), '\n')
    cat('\n')
    tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(vcf))])
    write.tree(file=paste0("tree_data/tree.", args[3],".indsremoved.nwk"),tree)

}


#get rows with complete data only
vcf[samples]<-sapply(vcf[samples], as.numeric)
complete_vcf<-vcf
complete_vcf[samples][complete_vcf[samples]=='N']<-NA
complete_vcf<-complete_vcf[complete.cases(complete_vcf[samples]),]


#get vcf positions with missing data
vcf_with_missing<-vcf[!vcf$POS %in% complete_vcf$POS,]


cat(paste0("    Number of SNPs with missing data: ", dim(vcf_with_missing)[1]),'\n')
cat(paste0("    Number of SNPs with no missing data: ", dim(complete_vcf)[1]),'\n')

tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)

packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')


edges<-data.frame(tree$edge)
colnames(edges)<-c('pos1','pos2')
edges$edge<-rownames(edges)
snps<-read.table(paste0(packpwd,"/../data/snps_isogg2019.txt"))



cat(paste0("\t","Assigning SNPs to branches."),'\n\n')


ref_alleleCountTracker<-0
alt_alleleCountTracker<-0

REFpos<-list()
ALTpos<-list()

for (edge in edges$edge){
    relevant_node<-edges$pos2[edges$edge==edge]
    desc<-tree$tip.label[getDescendants(tree, relevant_node)][!is.na(tree$tip.label[getDescendants(tree,relevant_node)])]
    nondesc<- samples[!samples %in% desc]
    REFs<-which(rowSums(complete_vcf[desc])==0 & rowSums(complete_vcf[nondesc])==length(nondesc))
    ALTs<-which(rowSums(complete_vcf[desc])==length(desc) & rowSums(complete_vcf[nondesc])==0)
    REFpos[[edge]]<-unique(complete_vcf$POS[REFs])
    ALTpos[[edge]]<-unique(complete_vcf$POS[ALTs])
    
    ref_alleleCountTracker<-length(REFs)+ref_alleleCountTracker
    alt_alleleCountTracker<-length(ALTs)+alt_alleleCountTracker
    
    allele_count_update_message<-(paste0(edge,'/', length(edges$edge),' nodes;    found ' , ref_alleleCountTracker+alt_alleleCountTracker, ' branch defining alleles ' ,'(REFs=', ref_alleleCountTracker,' / ', 'ALTs=', alt_alleleCountTracker,')'))
    cat("\r",allele_count_update_message)
}

cat('\n\n\n')


if ((ref_alleleCountTracker+alt_alleleCountTracker)==0){
    stop("There are no informative SNPs in your data. Check your VCF file.\n")
}

# saveRDS(ALTpos, file=paste0('tree_data/',args[3],".derpos.RData"))
# saveRDS(REFpos, file=paste0('tree_data/',args[3],".ancpos.RData"))

write.table(unlist(ALTpos), file=paste0('tree_data/',args[3],".derpos.txt"), quote = F, row.names = F, col.names = F, sep='\t')
write.table(unlist(REFpos), file=paste0('tree_data/',args[3],".ancpos.txt"), quote = F, row.names = F, col.names = F, sep='\t')


edge_df<-make_edge_df(ALTpos, REFpos)

write.table(edge_df,file=paste0("tree_data/",args[3],".edge_df.txt"), quote=F, row.names=F, sep='\t')

cat(paste0("\t","table with markers assigned to each branch written to tree_data/", args[3],".edge_df.txt)"),'\n\n')


LongSNPtable <- makeLongSNPtable(ALTpos, REFpos)

writeBed()

makeReport()





