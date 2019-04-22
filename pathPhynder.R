# pathPhynder
# Author: Rui Martiniano
# Contact: rm890 [at] cam.ac.uk

suppressWarnings(suppressPackageStartupMessages(library("optparse")))
suppressWarnings(suppressPackageStartupMessages(library("phytools")))


option_list <- list(
    make_option(c("-s", "--step"), default="all", help="Specifies which step to run. Options:
    \t\t\t- assign - assigns SNPs to branches.
    \t\t\t- all - runs all steps to map ancient SNPs to branches (1,2,3).
    \t\t\t- 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
    \t\t\t- 2 or chooseBestPath - finds the best branch/node of the tree for each sample.
    \t\t\t- 3 or addToTree - adds ancients samples to tree.
    \t\t\t[default %default]"),

    make_option(c("-i","--input_tree"),
        help = "Input tree in Newick format. [required]"),
    make_option(c("-v","--input_vcf"), 
        help = "Input vcf. Only needed for SNP to branch assignment. Needs to be haploid."),
    make_option(c("-p","--prefix"),
        help = "Prefix for the data files associated with the tree.
        \tThese were previously generated in the branch assignment step. [required]"),

    make_option(c("-b","--bam_file"), 
        help = "Input bam file. [required]"),

    make_option(c("-l","--list_of_bam_files"), 
        help = "List of paths to bam files. [required]"),

    make_option(c("-r","--reference"), default="~/in_development/pathPhynder/data/reference_sequences/hs37d5_Y.fa.gz", 
        help = "Reference genome (fasta format). [default \"%default\"]"),

    make_option(c("-m","--mode"), default="conservative", 
        help = "Mode for filtering pileups.  Options: relaxed or conservative. [default \"%default\"]"),

    make_option(c("-t", "--maximumTolerance"), type="numeric", default=3, 
        help="Maximum number of ALT alleles tolerated while traversing the tree.
                If exceeded, the algorithm stops and switches to the next path. [default %default]"),

    make_option(c("-c", "--pileup_read_mismatch_threshold"), type="numeric", default=0.7, 
        help = "Mismatch threshold for accepting a variant (for cases where reads for both alleles are present in pileup).
        \tFor a variant to pass filtering, reads containing the most frequent allele have to occur at least  
        \tat x proportion of the total reads. 1 is the most stringent, 0.5 is the most relaxed. [default %default]"),

    make_option(c("-o", "--output_prefix"), type="character", default="bamFileName", 
        help = "Sample name. This only works if a single bam file is used as an input. [default %default]")

)


tmpstr<-system('bash -l',input=c("shopt -s expand_aliases","type pathPhynder"), intern=T)

packpwd<-paste0(gsub('pathPhynder.R','',gsub('\'','',gsub('.*.Rscript ','',tmpstr))),'R')


# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))


# #test if arguments were given
if (is.null(opt$input_tree)){
    stop("Please provide the necessary arguments. Pass the -h parameter for help.")
} else if (is.null(read.tree(opt$input_tree))){
    stop("Please provide valid tree Newick file (-i).")
} else if (opt$step != "assign" & is.null(opt$bam_file) & is.null(opt$list_of_bam_files)){
    stop("Please provide either a bam file (-b) or a list of bam files (-l).")
} else if (opt$step != "assign" & length(opt$bam_file)>0 & length(opt$list_of_bam_files)>0){
    stop("Please provide either a bam file (-b) or a list of bam files (-l), not both.")
} else if (is.null(opt$prefix)){
    stop("Please provide a tree data prefix.")
}

if (opt$mode!="conservative" & opt$mode!="relaxed"){
    stop("The --mode parameter needs to be either conservative or relaxed.")
}



#test if bam files are bams
checkBamIntegrity<-function(bam_file){
    #test if bam files are bams
    if (file_test("-f", bam_file)==F){
        stop(paste0(bam_file," bam file does not exist."))
    } else {
        magic = readChar(gzfile(bam_file, 'r'), 4)
        if (!identical(magic, 'BAM\1')){
            stop(paste0(bam_file," is not a valid bam file."))
        }
    }
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


#############
#test if files exist
#############

#decide on the type of input
if (opt$step != "assign" & is.null(opt$list_of_bam_files)){
    if (opt$step != "assign" & file_test("-f", opt$bam_file)==F){
        stop("Please provide an existing bam file")
    } else {
        checkBamIntegrity(opt$bam_file)
    }
    input_type="bam_file"
}

if (opt$step != "assign" & is.null(opt$bam_file)){
    if (opt$step != "assign" & file_test("-f", opt$list_of_bam_files)==F){
        stop("Please provide an existing list of bam files")
    } else {

        checkBamListIntegrity(opt$list_of_bam_files)
    }
    input_type="bam_list"
}

if (opt$step != "assign" & file_test("-f", opt$reference)==F){
    stop("Please provide an existing reference genome fasta file.")
} else if (opt$step != "assign" & file_test("-f", paste0(opt$prefix,'.sites.bed'))==F){
    stop("Please provide a prefix for existing tree data.")
} 




#decide whether to use "chrY" or "Y"
if (opt$step != "assign" & length(grep("hg19", opt$reference))>0){
    chromosome_name<-"chrY"
} else if (opt$step != "assign" & length(grep("hs37d5", opt$reference))>0){
    chromosome_name<-"Y"
} else {
    if (opt$step != "assign") {
        stop("Your reference genome needs to be named hg19 or hs37d5")
    }
}



#print parameters into terminal
cat("\n\tpathPhynder v.0.0 \n", "\n\tParameters:\n")
cat("\n\t--input_tree ", opt$input_tree)
cat("\n\t--prefix ", opt$prefix)

if (opt$step != "assign") {
   if (input_type=="bam_file"){
      cat("\n\t--bam_file ", opt$bam_file)
    } else if (input_type=="bam_list"){
      cat("\n\t--list_of_bam_files ", opt$list_of_bam_files)
    } 
    cat("\n\t--reference ", opt$reference)
    cat("\n\t--mode ", opt$mode)
    cat("\n\t--maximumTolerance ", opt$maximumTolerance)
    cat("\n\t--pileup_read_mismatch_threshold ", opt$pileup_read_mismatch_threshold)


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




if( opt$step == "assign") {
    cat("assign SNPs to branches \n")
    system(paste("Rscript", paste0(packpwd,"/assign_noNA.R"), opt$input_tree, opt$input_vcf, opt$prefix))


} else if( opt$step == "all") {
    
    cat("All steps.\n")


    if (input_type=="bam_file"){

        system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_file,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file' , sample_name))

        system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',opt$bam_file,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))

    } else if (input_type=="bam_list"){

        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

        for(samp in bam_list$V1){
            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]

            system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),samp,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file', sample_name ))
            system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',sample_name,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))

        }
    }


#add step 3 here

} else if(opt$step == "pileup_and_filter" | opt$step == 1) {

    cat("Running pileup_and_filter\n\n")

    if (input_type=="bam_file"){
        system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_file,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file', sample_name ))
    } else if (input_type=="bam_list"){
        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

        for(samp in bam_list$V1){
            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]

            system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),samp,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold, 'bam_file', sample_name ))
        }
    }

} else if(opt$step == "chooseBestPath" | opt$step == 2) {

    cat("Running chooseBestPath\n")

    if (input_type=="bam_file"){
        system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',opt$bam_file,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))
    } else if (input_type=="bam_list"){
        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

            for(samp in bam_list$V1){
            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]

            system(paste("Rscript", paste0(packpwd,"/chooseBestPath.R"),opt$input_tree,opt$prefix,paste0('intree_folder/',sample_name,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))
        }
    }
    
#add step 3 here

# else if(opt$step == "addToTree" | opt$step == 2) {

#     cat("Running addToTree\n")
#     system(paste("Rscript", paste0(packpwd,"/addToTree.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))

# }


} else if(opt$step == "addToTree" | opt$step == 2) {

    cat("Running addToTree\n")

    if (input_type=="bam_file"){
        system(paste("Rscript", paste0(packpwd,"/workInProgres"),opt$input_tree,opt$prefix,paste0('intree_folder/',opt$bam_file,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))
    } else if (input_type=="bam_list"){
        bam_list<-read.table(opt$list_of_bam_files, stringsAsFactors=F)

            for(samp in bam_list$V1){
            sample_name<-unlist(strsplit(samp,'\\/'))[as.numeric(length(unlist(strsplit(samp,'\\/'))))]

            system(paste("Rscript", paste0(packpwd,"/workInProgres"),opt$input_tree,opt$prefix,paste0('intree_folder/',sample_name,'.intree.txt'), 'results_folder',  opt$maximumTolerance, sample_name ))
        }
    }
    
#add step 3 here

# else if(opt$step == "addToTree" | opt$step == 2) {

#     cat("Running addToTree\n")
#     system(paste("Rscript", paste0(packpwd,"/addToTree.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))

# }


} else {
    stop("Choose the step you would like to run.
        Options:
        \t- all - runs all steps of the analysis.
        \t- 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
        \t- 2 or chooseBestPath - finds the best branch/node of the tree for each sample.
        \t- 3 or addToTree - adds ancients samples to tree.")
}



cat("\n")



