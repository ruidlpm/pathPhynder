
get_vcf <- function(vcf_name){
    all_content = readLines(vcf_name)
    skip = all_content[-c(grep("CHROM",all_content))]
    vcf <- read.table(textConnection(skip), stringsAsFactors=F)
    header <- unlist(strsplit(all_content[grep("CHROM", all_content)[length(grep("CHROM", all_content))]], '\t'))
    colnames(vcf) <- make.names(header)
    #if T alleles read as TRUE, convert to character T.
    vcf$REF[vcf$REF==TRUE]<-"T"
    vcf$ALT[vcf$ALT==TRUE]<-"T"
    all_content<-NULL
    return(vcf)
}




createNAs <- function (x, pctNA = 0.1) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  NAloc[sample.int(n * p, floor(n * p * pctNA))] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}





vcf<-get_vcf("simulated_data.vcf")

set.seed(2)

vcf[10: dim(vcf)[2]] <- createNAs(vcf[10: dim(vcf)[2]], 0.1)

vcf[is.na(vcf)] <-'.'

write.table(vcf, "simulated_data.miss0.1.vcf", quote=F, row.names=F, sep='\t')


