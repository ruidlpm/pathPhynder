# pileup_and_filter.R
# runs pileup and filters bases, outputs intree files

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=9) {
    stop("  Arguments needed.\n
	\tusage
	\tRscript pileup_and_filter.R <bam_list> <data_prefix> <folder_out> <refgen_path> <mode>['conservative'/'relaxed'] <chromosome_name>>['chrY'/'Y'] <pileup_read_mismatch_threshold>[default 0.7] <bam file> or <list of bams> <out_prefix>
    	", call.=FALSE)
}


tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)

packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')


sites_data<-args[2] 
intree_folder<-args[3]
refgen_path<-args[4]
mode<-args[5]
chromosome_name<-args[6]
pileup_read_mismatch_threshold<-args[7]
input_type<-args[8]
out_prefix<-args[9]

dir.create(intree_folder, showWarnings = FALSE)

if (chromosome_name=="chrY"){
	sites_var<-paste0(sites_data,".siteschr.bed")
} else if (chromosome_name=="Y") {
	sites_var<-paste0(sites_data,".sites.bed")
}


for (testfile in c(args[1], sites_var, args[4])){
	if (!file.exists(testfile)) {
		stop(paste(testfile, "- does this file exist?"))
	}
}


if (input_type=="bam_file"){
		file_path<-args[1]

		sample_name<-unlist(strsplit(file_path,'\\/'))[as.numeric(length(unlist(strsplit(file_path,'\\/'))))]
		if (out_prefix!="sample"){
			sample_name<-out_prefix
		}
		if (!file.exists(file_path)) {
			print(paste(file_path, "- does this file exist?"))
		} else {
			pos_var<-paste0(sites_data,".sites.txt")
			
			cmd=paste0("samtools mpileup ",
				paste0(file_path),
				" --ignore-RG --positions ",
				sites_var,
				paste0(" -f ",refgen_path),
				paste0(">", intree_folder,'/',sample_name,".pileup"))
			
			system(cmd, wait=T)
		
			cmd2=paste0("python3 ",gsub("R$","",packpwd), "inst/python/call_bases_chrY_v2.1.py ",
				paste0(" -i ", intree_folder,'/',sample_name,".pileup"),
				paste0(" -m ", mode),
				paste0(" -o ", intree_folder,'/',sample_name,".intree.txt"),
				paste0(" -t ", pos_var),
				paste0(" -c ", pileup_read_mismatch_threshold)
				)
			system(cmd2, wait=T)
		}	
}	
