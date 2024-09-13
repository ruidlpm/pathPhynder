#!/usr/bin/env Rscript
# pathPhynder
# Author: Rui Martiniano
# Contact: ruidlpm [at] gmail.com

suppressWarnings(suppressPackageStartupMessages(library("optparse")))
suppressWarnings(suppressPackageStartupMessages(library("phytools")))

current_version<-'Version: 1.2.2'

tmparg <- commandArgs(trailingOnly = F)  
scriptPath <- normalizePath(dirname(sub("^--file=", "", tmparg[grep("^--file=", tmparg)])))
# tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)
packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',scriptPath))),'/R')

# source(paste0(packpwd,'/pathPhynder_likelihood_functions.R'))

st=format(Sys.time(), "%Y-%m-%d_%H:%M")
logname<-paste("log.",st, ".txt", sep = "")

sink(logname,append = T, type = c("output", "message"),
          split = FALSE)

# # pathPhynder  -i ~/software/pathPhynder/data/BigTree_Y/bigtree_annotated_V1.nwk -p ~/Ycap_project/Analysis/placement_bamq25/tree_data/BigTree_Y_data -b ../BAMs/GL107A.merged.sorted.rmdup.q25.rg.realigned.bam -s all -t 100 -G ~/software/pathPhynder/data/200803.snps_isogg.txt -r ~/software/pathPhynder/data/reference_sequences/hs37d5_Y.fa


option_list <- list(
    
   # \t\t\t- assign - assigns SNPs to branches.

    make_option(c("-s", "--step"), default="all", help="Specifies which step to run. Options:
    \t\t\t- prepare - prepares files for running pathPhynder.
    \t\t\t- all - runs all steps to map ancient SNPs to branches (1,2,3).
    \t\t\t- 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
    \t\t\t- 2 or chooseBestPath - finds the best branch/node of the tree for each sample.
    \t\t\t- 3 or addAncToTree - adds ancients samples to tree.
    \t\t\t[default %default]"),

    make_option(c("-i","--input_tree"), help = "Input tree in Newick format. [required]"),

    # make_option(c("-v","--input_vcf"), 
    #     help = "Input haploid vcf. Needed for SNP to branch assignment and if estimating likelihoods."),
    make_option(c("-f","--branches_file"), 
        help = "branches.snp file - SNP placement file created with phynder."),

    make_option(c("-p","--prefix"),
        help = "Prefix for the data files associated with the tree.
        \tThese were previously generated in the branch assignment step. [required]"),

    make_option(c("-b","--bam_file"), 
        help = "Input bam file. [required]"),

    make_option(c("-l","--list_of_bam_files"), 
        help = "List of paths to bam files. [required]"),

    make_option(c("-r","--reference"), default=paste0(packpwd,"/../data/reference_sequences/hs37d5_Y.fa.gz"),
        help = "Reference genome (fasta format). [default \"%default\"]"),

    make_option(c("-m","--filtering_mode"), default="default", 
        help = "Mode for filtering pileups.  Options: default, no-filter or transversions. [default \"%default\"]"),

    make_option(c("-t", "--maximumTolerance"), type="numeric", default=3, 
        help="Maximum number of ALT alleles tolerated while traversing the tree.
                If exceeded, the algorithm stops and switches to the next path. [default %default]"),
    
    make_option(c("-q", "--baseQuality"), type="numeric", default=20, 
        help="Minimum base quality for samtools mpileup. [default %default]"),

    make_option(c("-c", "--pileup_read_mismatch_threshold"), type="numeric", default=0.7, 
        help = "Mismatch threshold for accepting a variant (for cases where reads for both alleles are present in pileup).
        \tFor a variant to pass filtering, reads containing the most frequent allele have to occur at least  
        \tat x proportion of the total reads. 1 is the most stringent, 0.5 is the most relaxed. [default %default]"),

    make_option(c("-o", "--output_prefix"), type="character", default="bamFileName", 
        help = "Sample name. This only works if a single bam file is used as an input. [default %default]"),

    make_option(c("-G", "--haplogroups"), type="character", default='none', 
        help = "List of known haplogroup-defining SNPs")

)



# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))

base_qual<-opt$baseQuality

#check if bams are ok
checkBamIntegrity<-function(bam_file){
    #test if bam files are bams
    if (file_test("-f", bam_file)==F){
        stop(paste0(bam_file," bam file does not exist."))
    } else {
        con=gzfile(bam_file, 'r')
        magic = readChar(con, 4)
        if (!identical(magic, 'BAM\1')){
            close(con)
            stop(paste0(bam_file," is not a valid bam file."))
        }
    }
    close(con)
}


checkBamListIntegrity<-function(list_bams){
    if (file_test("-f", list_bams)==F){
        stop(paste0(list_bams," bam file list does not exist."))
    } else {
        bam_file_list<-read.table(list_bams)
        for (bam in bam_file_list$V1){
            checkBamIntegrity(bam)
        }
    }
}




# #test if arguments were given
if (is.null(opt$input_tree)){
    stop("Please provide the necessary arguments. Pass the -h parameter for help.\n
        Usage: pathPhynder -s all -i tree.nwk -p tree_data/<prefix_output> -b <sample.bam>\n\n")
} else if (is.null(read.tree(opt$input_tree))){
    stop("Please provide valid tree Newick file (-i).")
} else if (opt$step != "prepare" & is.null(opt$bam_file) & is.null(opt$list_of_bam_files)){
    stop("Please provide either a bam file (-b) or a list of bam files (-l).")
} else if (opt$step != "prepare" & length(opt$bam_file)>0 & length(opt$list_of_bam_files)>0){
    stop("Please provide either a bam file (-b) or a list of bam files (-l), not both.")
} else if (is.null(opt$prefix)){
    stop("Please provide a tree data prefix.")
}

if (opt$filtering_mode!="default" & opt$filtering_mode!="no-filter" & opt$filtering_mode!="transversions"){
    stop("The --filtering_mode parameter needs to be either default, no-filter or transversions.")
}


#decide on the type of input, bam or bam file
if (opt$step != "prepare" & is.null(opt$list_of_bam_files)){
    if (opt$step != "prepare" & file_test("-f", opt$bam_file)==F){
        stop("Please provide an existing bam file")
    } else {
        checkBamIntegrity(opt$bam_file)
    }
    input_type="bam_file"
}

if (opt$step != "prepare" & is.null(opt$bam_file)){
    if (opt$step != "prepare" & file_test("-f", opt$list_of_bam_files)==F){
        stop("Please provide an existing list of bam files")
    } else {

        checkBamListIntegrity(opt$list_of_bam_files)
    }
    input_type="bam_list"
}

if (opt$step != "prepare" & file_test("-f", opt$reference)==F){
    stop("Please provide an existing reference genome fasta file.")
} else if (opt$step != "prepare" & file_test("-f", paste0(opt$prefix,'.sites.bed'))==F){
    stop("Please provide a prefix for existing tree data.")
} 



# #decide whether to use "chrY" or "Y"
# if (opt$step != "assign" & length(grep("hg19", opt$reference))>0){
#     chromosome_name<-"chrY"
# } else if (opt$step != "assign" & length(grep("hs37d5", opt$reference))>0){
#     chromosome_name<-"Y"
# } else {
#     if (opt$step != "assign") {
#         stop("Your reference genome needs to be named hg19 or hs37d5")
#     }
# }


if (opt$step != "prepare" & opt$haplogroups!='none'){
    if (file_test("-f", opt$haplogroups)==F){
        stop("When passing the parameter -G/--haplogroups, please provide an existing known haplogroups file")
    } else {
        haplogroups_file=opt$haplogroups
        call_haps=T
    }
}


#print parameters into terminal
cat("\n", paste0("Program: pathPhynder","\n"))
cat("", paste0("Version: ", current_version,"\n"))
cat("", paste0("Contact: ruidlpm@gmail.com","\n", "\n\tParameters:\n"))
cat("\n\tBest path mode")
cat("\n\t--input_tree ", opt$input_tree)
cat("\n\t--prefix ", opt$prefix)

if (opt$step != "prepare") {
   if (input_type=="bam_file"){
      cat("\n\t--bam_file ", opt$bam_file)
    } else if (input_type=="bam_list"){
      cat("\n\t--list_of_bam_files ", opt$list_of_bam_files)
    } 
    cat("\n\t--reference ", opt$reference)
    cat("\n\t--filtering_mode ", opt$filtering_mode)
    cat("\n\t--maximumTolerance ", opt$maximumTolerance)
    cat("\n\t--pileup_read_mismatch_threshold ", opt$pileup_read_mismatch_threshold)

    if (opt$haplogroups!='none'){
        cat("\n\t--haplogroups", opt$haplogroups, "\n")
    }

    if (input_type=="bam_file"){
        if (opt$output_prefix=="bamFileName"){
            sample_name<-unlist(strsplit(opt$bam_file,'\\/'))[as.numeric(length(unlist(strsplit(opt$bam_file,'\\/'))))]
            cat("\n\t--output_prefix ", sample_name,'\n\n')
        } else {
            cat("\n\t--output_prefix ", opt$output_prefix,'\n\n')
            sample_name<-opt$output_prefix
        }
    }

}



cat("\n")

if (opt$step != "prepare" ){
    incon <- gzcon(file(opt$reference,open="rb"))
    chromosome_name<-gsub('^>','',readLines(incon,1))
    cat(paste0('Contig name is ',chromosome_name, '\n'))
    if (length(chromosome_name) == 0) {
        stop("Please fix your reference genome.")
    }
}

if(opt$step == "prepare") {
    cat("Preparing files for running pathPhynder \n")
    if (is.null(opt$branches_file) ){
        stop('Please provide branches.snp file created with phynder.')
    }
    if (file_test("-f", opt$branches_file)==F){
        stop('Please provide branches.snp file created with phynder.')
    }
    
    system(paste("Rscript", paste0(packpwd,"/prep_files.R"), opt$input_tree, opt$branches_file, opt$prefix, opt$haplogroups))




} else if( opt$step == "all") {
    
    cat("All steps.\n")

    if (input_type=="bam_file"){


        cat(paste0('Sample name:', opt$bam_file, '\n'))

        system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_file,opt$prefix,'intree_folder', opt$reference, opt$filtering_mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file' , sample_name, opt$haplogroups,base_qual))

        system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',opt$bam_file,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name, opt$haplogroups ))

    } else if (input_type=="bam_list"){

        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

        for(samp in bam_list$V1){

            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]
            cat(paste0('Sample name: ', sample_name, '\n'))

            system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),samp,opt$prefix,'intree_folder', opt$reference, opt$filtering_mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file', sample_name ,opt$haplogroups,base_qual))
            
            system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',sample_name,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name , opt$haplogroups ))
      
        }
    }

    system(paste("Rscript", paste0(packpwd,"/addAncToTree.R"),opt$input_tree, 'results_folder',opt$prefix))


# # pathPhynder  -i ~/software/pathPhynder/data/BigTree_Y/bigtree_annotated_V1.nwk -p ~/Ycap_project/Analysis/placement_bamq25/tree_data/BigTree_Y_data -b ../BAMs/GL107A.merged.sorted.rmdup.q25.rg.realigned.bam -s 1 -t 100 -G ~/software/pathPhynder/data/200803.snps_isogg.txt -r ~/software/pathPhynder/data/reference_sequences/hs37d5_Y.fa


} else if(opt$step == "pileup_and_filter" | opt$step == 1) {

    cat("Running pileup_and_filter\n\n")

    if (input_type=="bam_file"){
        cat(paste0('Sample name:', opt$bam_file, '\n'))

        system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_file,opt$prefix,'intree_folder', opt$reference, opt$filtering_mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file', sample_name, opt$haplogroups,base_qual))
    
    } else if (input_type=="bam_list"){

        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

        for(samp in bam_list$V1){

            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]
            cat(paste0('Sample name: ', sample_name, '\n'))

            system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),samp,opt$prefix,'intree_folder', opt$reference, opt$filtering_mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file', sample_name ,opt$haplogroups,base_qual))
        
        }
    }

} else if(opt$step == "chooseBestPath" | opt$step == 2) {

    cat("Running chooseBestPath\n")

    if (input_type=="bam_file"){
        cat(paste0('Sample name:', opt$bam_file, '\n'))

        system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',opt$bam_file,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name , opt$haplogroups))
    
    } else if (input_type=="bam_list"){

        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

            for(samp in bam_list$V1){

            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]
            cat(paste0('Sample name: ', sample_name, '\n'))

            system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',sample_name,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name , opt$haplogroups))
      
        }
    }
    
} else if(opt$step == "addAncToTree" | opt$step == 3) {

    cat("Running addAncToTree\n")
    system(paste("Rscript", paste0(packpwd,"/addAncToTree.R"),opt$input_tree, 'results_folder',opt$prefix))

} else {
    
    stop("Choose the step you would like to run.
        Options:
        \t- prepare - prepares files for running pathPhynder based on the SNP placement from phynder.
        \t- all - runs all steps of the analysis.
        \t- 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
        \t- 2 or chooseBestPath - finds the best branch/node of the tree for each sample.
        \t- 3 or addAncToTree - adds ancients samples to tree.")

}

sink()


cat("\n")


