get_vcf<-function(vcf_name){
	all_content = readLines(vcf_name)
	skip = all_content[-c(grep("CHROM",all_content))]
	vcf <- read.table(textConnection(skip), stringsAsFactors=F)
	header<-unlist(strsplit(all_content[grep("CHROM", all_content)[length(grep("CHROM", all_content))]], '\t'))
    colnames(vcf)<-make.names(header)
    #if T alleles read as TRUE, convert to character T.
    vcf$REF[vcf$REF==TRUE]<-"T"
    vcf$ALT[vcf$ALT==TRUE]<-"T"
	all_content<-NULL
	return(vcf)
}











choose_best_node<-function(x){
  return(best_node)
}

estimate_probability<-function(table_counts){
  return(table_w_prop)
}

plot_tree_ancder_count<-function(x){
  plot()
}

make_output<-function(table_w_prop){
  return(out)
}
