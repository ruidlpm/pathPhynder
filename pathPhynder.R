# pathPhynder
# Author: Rui Martiniano
# Contact: rm890 [at] cam.ac.uk

suppressPackageStartupMessages(library("optparse"))


option_list <- list(
    make_option(c("-s", "--step"), default="all", help="Specifies which step to run. Options:
    \t\t\t- all - runs all steps of the analysis.
    \t\t\t- 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
    \t\t\t- 2 or ancient_SNPs_to_branches - adds anc/der alleles to the input tree.
    \t\t\t- 3 or decide - finds the best branch/node of the tree for each sample.
    \t\t\t[default %default]"),

    make_option(c("-i","--input_tree"),
        help = "Input tree in Newick format. [required]"),
    # make_option(c("-v","--input_vcf"), default="data/test_tree.vcf", 
    #     help = "Input vcf. Needs to be haploid. [default \"%default\"]"),
    make_option(c("-p","--prefix"),
        help = "Prefix for the data files associated with the tree.
        \tThese were previously generated in the branch assignment step. [required]"),

    make_option(c("-b","--bam_list"), 
        help = "List of paths to bam files. [required]"),

    make_option(c("-r","--reference"), default="data/reference_sequences/hs37d5_Y.fa.gz", 
        help = "Reference genome (fasta format). [default \"%default\"]"),

    make_option(c("-m","--mode"), default="conservative", 
        help = "Mode for filtering pileups.  Options: relaxed or conservative. [default \"%default\"]"),

    make_option(c("-t", "--maximumTolerance"), type="integer", default=3, 
        help="Maximum number of ALT alleles tolerated while traversing the tree.
                If exceeded, the algorithm stops and switches to the next path. [default %default]"),

    make_option(c("-c", "--pileup_read_mismatch_threshold"), type="numeric", default=0.7, 
        help = "Mismatch threshold for accepting a variant (for cases where reads for both alleles are present in pileup).
        \tFor a variant to pass filtering, reads containing the most frequent allele have to occur at least  
        \tat x proportion of the total reads. 1 is the most stringent, 0.5 is the most relaxed. [default %default]")

)



packpwd<-("~/in_development/pathPhynder/R")

# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))


#decide whether to use "chrY" or "Y"
if (length(grep("hg19", opt$reference))>0){
    chromosome_name<-"chrY"
} else if (length(grep("hs37d5", opt$reference))>0){
    chromosome_name<-"Y"
} else {
    stop("Your reference genome needs to be named hg19 or hs37d5")
}


cat("\n\tpathPhynder v.0.0 \n", "\n\tParameters:\n")
cat("\n\t--input_tree ", opt$input_tree)
cat("\n\t--prefix ", opt$prefix)
cat("\n\t--bam_list ", opt$bam_list)
cat("\n\t--reference ", opt$reference)
cat("\n\t--mode ", opt$mode)
cat("\n\t--maximumTolerance ", opt$maximumTolerance)
cat("\n\t--pileup_read_mismatch_threshold ", opt$pileup_read_mismatch_threshold,'\n\n')




#test if files exist
if (file_test("-f", opt$input_tree)==F){
    stop("Please provide an existing tree Newick file.")
} else if (file_test("-f", opt$bam_list)==F){
    stop("Please provide an existing bam list.")
} else if (file_test("-f", opt$reference)==F){
    stop("Please provide an existing reference genome fasta file.")
} else if (file_test("-f", paste0(opt$prefix,'.sites.bed'))==F){
    stop("Please provide a prefix for existing tree data.")
} 


#test if bam files are bams
bam_list<-read.table(opt$bam_list)
for (bam in bam_list$V1){
    if (file_test("-f", bam)==F){
        stop(paste0(bam," bam file does not exist. Fix the sample list."))
    } else {
        magic = readChar(gzfile(bam, 'r'), 4)
        if (!identical(magic, 'BAM\1')){
            stop(paste0(bam," is not a valid bam file."))
        }
    }
}




# do some operations based on user input
if( opt$step == "all") {
    
    cat("All steps.\n")

    system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold ))

    system(paste("Rscript", paste0(packpwd,"/ancient_SNPs_to_branches.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))

    system(paste("Rscript", paste0(packpwd,"/decide_len.R"),opt$input_tree, 'results_folder', opt$maximumTolerance, opt$prefix))

} else if(opt$step == "pileup_and_filter" | opt$step == 1) {

    cat("Running pileup_and_filter\n\n")
    system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold))

} else if(opt$step == "ancient_SNPs_to_branches" | opt$step == 2) {

    cat("Running ancient_SNPs_to_branches\n")
    system(paste("Rscript", paste0(packpwd,"/ancient_SNPs_to_branches.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))

} else if(opt$step == "decide" | opt$step == 3) {

    cat("Running decide\n")
    system(paste("Rscript", paste0(packpwd,"/decide_len.R"),opt$input_tree, 'results_folder', opt$maximumTolerance, opt$prefix))

} 

cat("\n")



