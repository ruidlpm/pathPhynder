#!/usr/bin/env Rscript


require(data.table, quietly = TRUE)
require(phytools, quietly = TRUE)

cat('\n\n',"Preparing sites for likelihood estimation", '\n\n\n')


args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=2) {
    stop("\tArguments needed.\n
    \tusage:
    \tprepare_sites_likelihood.R <input.vcf> <out prefix>", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("prepare_sites_likelihood.R", args[1], args[2]), '\n\n')
}


dir.create('tree_data', showWarnings = FALSE)


# for likes
makeVCFsites<-function(vcf){
    site_info<-data.frame(chr=gsub('chr','',vcf$X.CHROM),NA,NA, vcf$POS, paste0(vcf$REF,'->', vcf$ALT), vcf$REF, vcf$ALT)
    bedfile_data<-data.frame(chr=gsub('chr','',vcf$X.CHROM), vcf$POS-1, vcf$POS)
    bedfile_data_w_chr<-data.frame(chr=vcf$X.CHROM, vcf$POS-1, vcf$POS)
    write.table(file=paste0('tree_data/',args[2],'.sites.bed'),bedfile_data, quote = F, row.names = F, col.names = F, sep='\t')
    write.table(file=paste0('tree_data/',args[2],'.siteschr.bed'),bedfile_data_w_chr, quote = F, row.names = F, col.names = F, sep='\t')
    write.table(file=paste0('tree_data/',args[2],'.sites.txt'),site_info, quote = F, row.names = F, col.names = F, sep='\t')
    cat(paste0("written ",dim(unique(bedfile_data))[1]," positions for variant calling."))
    cat('\n')
}


if(args[2]=='.likes'){
    stop('Please provide an output prefix -p')
}

vcf_name<-args[1]

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

nucleotides=c('A','T','C','G')
old_POS<-vcf$POS
vcf<-vcf[vcf$REF %in% nucleotides,]
vcf<-vcf[vcf$ALT %in% nucleotides,]
new_POS<-vcf$POS


excluded_nonSingleNuc<-old_POS[!old_POS %in% new_POS]
cat(paste0('excluded ',length(old_POS)-length(new_POS), ' non ATCG sites from analysis.'))
cat('\n') 


makeVCFsites(vcf)


