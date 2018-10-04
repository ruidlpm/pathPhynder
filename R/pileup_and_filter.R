

# pileup_and_filter.R
# runs pileup and filters bases, outputs intree files

cat('\n\n',"pileup_and_filter", '\n\n\n')



args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)!=6) {
    stop("  Arguments needed.\n
	\tusage
	\tRscript pileup_and_filter.R <bam_list> <data_prefix> <folder_out> <refgen_path> <mode>['conservative'/'relaxed'] <chromosome_name>>['chrY'/'Y']
    	", call.=FALSE)
} else {
    cat("   Command used:",'\n\n')
    cat(paste("pileup_and_filter.R", args[1], args[2],args[3]),args[4],args[5],args[6], '\n\n')
}


sites_data<-args[2] 
intree_folder<-args[3]
mode<-args[5]
chromosome_name<-args[6]

dir.create(intree_folder, showWarnings = FALSE)



if (chromosome_name=="chrY"){
	sites_var<-paste0("tree_data/",sites_data,".siteschr.bed")
} else if (chromosome_name=="Y") {
	sites_var<-paste0("tree_data/",sites_data,".sites.bed")
}


for (testfile in c(args[1], sites_var, args[4])){
	if (!file.exists(testfile)) {
		stop(paste(testfile, "- does this file exist?"))
	}
}

bam_list<-read.table(args[1])
refgen_path<-args[4]






for(samp in bam_list$V1){
	file_path<-samp
	sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]
	
	if (!file.exists(file_path)) {
		print(paste(file_path, "- does this file exist?"))
	} else {

		
		pos_var<-paste0("tree_data/",sites_data,".sites.txt")
		
		cmd=paste0("samtools mpileup ",
			paste0(file_path),
			" --ignore-RG --positions ",
			sites_var,
			paste0(" -f ",refgen_path),
			paste0(">", intree_folder,'/',sample_name,".pileup"))
		
		print(cmd)
		system(cmd, wait=T)
	
		
		cmd2=paste0("python3 call_bases_chrY_v2.1.py ",
			paste0(" -i ", intree_folder,'/',sample_name,".pileup"),
			paste0(" -m ", mode),
			paste0(" -o ", intree_folder,'/',sample_name,".intree.txt"),
			paste0(" -t ", pos_var))
		
		print(cmd2)
		system(cmd2, wait=T)
	}	
}
