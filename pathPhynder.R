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
    
    make_option(c("-i","--input_tree"), default="data/test_tree.nwk",
        help = "Input tree in Newick format. [default \"%default\"]"),
    # make_option(c("-v","--input_vcf"), default="data/test_tree.vcf", 
    #     help = "Input vcf. Needs to be haploid. [default \"%default\"]"),
    make_option(c("-p","--prefix"), default="data/test_tree", 
        help = "Prefix for the data files associated with the tree.
        \tThese were previously generated in the branch assignment step. [default \"%default\"]"),

    make_option(c("-b","--bam_list"), default="sample_list.txt", 
        help = "List of paths to bam files. [default \"%default\"]"),

    make_option(c("-r","--reference"), default="data/hg19.fa", 
        help = "Reference genome (fasta format). [default \"%default\"]"),

    make_option(c("-m","--mode"), default="conservative", 
        help = "Mode for filtering pileups.  Options: relaxed or conservative. [default \"%default\"]"),

    make_option(c("-t", "--maximumTolerance"), type="integer", default=3, 
        help="Maximum number of ALT alleles tolerated while traversing the tree.
                If exceeded, the algorithm stops and switches to the next path. [default %default]"),

    make_option(c("-c", "--pileup_read_mismatch_threshold"), type="numeric", default=0.7, 
        help = "Mismatch threshold for accepting a variant (for cases where reads for both alleles are present in pileup).
        \tFor a variant to pass filtering, reads containing the most frequent allele have to occur at least  
        \tat x proportion of the total reads. [default %default]")

)









packpwd<-("~/in_development/pathPhynder/R")

# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))

if (length(grep("hg19", opt$reference))>0){
    chromosome_name<-"chrY"
} else if (length(grep("hs37d5", opt$reference))>0){
    chromosome_name<-"Y"
} else {
    stop("Your reference genome needs to be named hg19 or hs37d5")
}


print(opt)
# do some operations based on user input
if( opt$step == "all") {
    cat("All steps.\n")

    system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold ))

    system(paste("Rscript", paste0(packpwd,"/ancient_SNPs_to_branches.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))

    system(paste("Rscript", paste0(packpwd,"/decide_len.R"),opt$input_tree, 'results_folder', opt$maximumTolerance, opt$prefix))

} else if(opt$step == "pileup_and_filter" | opt$step == 1) {
    cat("Running pileup_and_filter\n")
    system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, chromosome_name,opt$pileup_read_mismatch_threshold))
} else if(opt$step == "ancient_SNPs_to_branches" | opt$step == 2) {
    cat("Running ancient_SNPs_to_branches\n")
    system(paste("Rscript", paste0(packpwd,"/ancient_SNPs_to_branches.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))
} else if(opt$step == "decide" | opt$step == 3) {
    cat("Running decide\n")
    system(paste("Rscript", paste0(packpwd,"/decide_len.R"),opt$input_tree, 'results_folder', opt$maximumTolerance, opt$prefix))
} 
cat("\n")


# bam_list.txt





# input_phylogeny
# prefix
# intree_folder
# results_folder
# input.vcf
# bam_list
# refgen_path
# mode
# chromosome_name
# <Maximum Tolarance (int)>