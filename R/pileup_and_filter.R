# pileup_and_filter.R
# runs pileup and filters bases, outputs intree files

args = commandArgs(trailingOnly=TRUE)

# test if no args are given
if (length(args)<10) {
    stop("  Arguments needed.\n
	\tusage
	\tRscript pileup_and_filter.R <bam_list> <data_prefix> <folder_out> <refgen_path> <mode>['default'/'no-filter'/'transversions'] <chromosome_name>>['chrY'/'Y'] <pileup_read_mismatch_threshold>[default 0.7] <bam file> or <list of bams> <out_prefix> haplogroups basequal
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
haplogroups<-args[10]
base_qual<-args[11]



dir.create(intree_folder, showWarnings = FALSE)

if (length(grep('chr',chromosome_name))==1){
	sites_var<-paste0(sites_data,".siteschr.bed")
} else if (length(grep('chr',chromosome_name))==0){
	sites_var<-paste0(sites_data,".sites.bed")
}

# sites_var<-paste0(sites_data,".sites.bed")

for (testfile in c(args[1], sites_var, args[4])){
	if (!file.exists(testfile)) {
		stop(paste(testfile, "- does this file exist?"))
	}
}

print(refgen_path)



if (args[10]!='none'){
	haplogroups<-args[10]
	a<-read.table(haplogroups, stringsAsFactors=F)
	a$V3<-as.numeric(a$V3)
	a<-a[!is.na(a$V3),]
	if (dim(a)[2]!=5){
		stop(paste0('expected 5 columns in haplogroup table, found ', dim(5)[2]))
	}
	hg_call_df<-data.frame(chromosome_name,a$V3-1, a$V3)
	write.table(file='hg_call_df.txt',hg_call_df, quote=F, row.names=F, col.names=F, sep='\t')
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
			
			cat('calling SNPs for tree placement\n')
			cmd=paste0("samtools mpileup -q 25 ", file_path,
				paste0(" --min-BQ ", base_qual), ' ' ,
				" --no-BAQ --ignore-RG --positions ",
				sites_var,
				paste0(" -f ",refgen_path),
				paste0(">", intree_folder,'/',sample_name,".pileup"))
			
			system(cmd, wait=T)
	print(cmd)
		
			if (args[10]!='none'){
				
				cat('calling haplogroup defining SNPs\n')
				cmd=paste0("samtools mpileup ",
					paste0(file_path),
					paste0(" --min-BQ ", base_qual),
					" --no-BAQ --ignore-RG --positions hg_call_df.txt",
					paste0(" -f ",refgen_path),
					paste0(">", intree_folder,'/',sample_name,".calls_hgs.pileup"))

				system(cmd, wait=T)
				# system('rm hg_call_df.txt')

			}
			#for likelihood estimation
			#if 'likes' in names of sites variable, then output with name 'status.txt'
			if (length(grep('likes', sites_var))>0){
				cmd_likes=paste0("python3 ",gsub("R$","",packpwd), "inst/python/call_bases_chrY_v2.1.py ",
					paste0(" -i ", intree_folder,'/',sample_name,".pileup"),
					paste0(" -m ", mode),
					paste0(" -o ", intree_folder,'/',sample_name,".status.txt"),
					paste0(" -t ", pos_var),
					paste0(" -c ", pileup_read_mismatch_threshold)
					)

				system(cmd_likes, wait=T)
			} else {
				#for best path estimation
				cmd_bestpath=paste0("python3 ",gsub("R$","",packpwd), "inst/python/call_bases_chrY_v2.1.py ",
					paste0(" -i ", intree_folder,'/',sample_name,".pileup"),
					paste0(" -m ", mode),
					paste0(" -o ", intree_folder,'/',sample_name,".intree.txt"),
					paste0(" -t ", pos_var),
					paste0(" -c ", pileup_read_mismatch_threshold),
					paste0(" -g ", args[10])
					)
				system(cmd_bestpath, wait=T)
		}	
	}
}
