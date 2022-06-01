
cat('\n\n',"Preparing files for sample placement", '\n\n\n')


args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=3) {
    stop("\tArguments needed.\n
    \tusage:
    \tRscript make_vcf.R <folder with intree files> <chromosome name> <output name>", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("make_vcf", args[1], args[2],args[3]), '\n\n')
}


files<-list.files(pattern='intree\\.txt', path=args[1])

contig<-args[2] 
outfile<-args[3] 

counter=0
to_excl<-NA
for (i in files){
	sample_name<-gsub('.intree.txt','',i)
	cat('processing',sample_name,'\n') 
	tmp<-try(read.table(paste0(args[1],'/', i)))

	if(class(tmp)!='try-error') {
		if (dim(tmp)[1]>0){
			colnames(tmp)<-c('POS', 'REF', 'ALT','anc_counts', 'der_counts', sample_name)
			tmp[[sample_name]][tmp[[sample_name]]=='-9']<-'.'
			tmp<-tmp[c('POS','REF', 'ALT',sample_name)]
			if (counter==0){
				counter=counter+1
				main_file<-tmp
				tmp<-NULL
			} else if (counter>0) {
				main_file<-merge(main_file, tmp,all=T)
				main_file[is.na(main_file)]<-'.'
			}
		}
	} else {
		to_excl<-c(to_excl,sample_name)
	}
}


sample_names<-colnames(main_file)[4:length(colnames(main_file))]
sample_names<-sample_names[!sample_names %in% to_excl]


main_file$CHROM<-contig
main_file$ID<-'.'
main_file$QUAL<-'.'
main_file$FILTER<-'.'
main_file$INFO<-'.'
main_file$FORMAT<-'GT'

main_file<-main_file[c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',sample_names)]



system(paste('echo \\#\\#fileformat=VCFv4.2 > ',outfile))
system(paste('echo \\#\\#FORMAT=\\<ID=GT,Number=1,Type=String,Description=\\"Genotype\\"\\> >>',outfile))
system(paste("echo", "\\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", paste(sample_names, collapse="\t"), " | tr -s ' ' |tr ' ' '\t' >>", outfile,  sep='\t'))

write.table(file=outfile,main_file, append=T, sep='\t', col.names=F, row.names=F, quote=F)

# system(paste('bgzip -f', outfile))
# system(paste('tabix -p vcf',  paste0(outfile, '.gz')))

cat('VCF file with', length(sample_names), 'samples and', dim(main_file)[1] ,'variants written to',outfile,'\n\n')
