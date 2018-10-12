 
suppressPackageStartupMessages(library("optparse"))


# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

# debug(getopt)


option_list <- list(
    make_option(c("-i","--input_tree"), default="data/test_tree.nwk",
        help = "Input tree in Newick format. [default \"%default\"]"),
    # make_option(c("-v","--input_vcf"), default="data/test_tree.vcf", 
    #     help = "Input vcf. Needs to be haploid. [default \"%default\"]"),
    make_option(c("-p","--prefix"), default="data/test_tree", 
        help = "Prefix for the data files associated with the tree. [default \"%default\"]"),
    make_option(c("-b","--bam_list"), default="data/bam_list.txt", 
        help = "List of paths to bam files. [default \"%default\"]"),
    make_option(c("-r","--reference"), default="data/hg19.fa", 
        help = "Reference genome. Needs to be in fasta format. [default \"%default\"]"),
    make_option(c("-m","--mode"), default="conservative", 
        help = "Mode for filtering pileups.  Options: relaxed or conservative. [default \"%default\"]"),
    make_option(c("-t", "--maximumTolerance"), type="integer", default=3, 
        help="Maximum number of ALT alleles tolerated while traversing the tree.
                If exceeded, the algorithm stops and switches to the next path. [default %default]"),
    make_option(c("-s", "--step"), default="all", 
        help="Specifies which step to run. [default %default] Options:
                 - all - runs all steps (pileup_and_filter, ancient_SNPs_to_branches and decide).
                 - 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
                 - 2 or ancient_SNPs_to_branches - adds anc/der alleles to the input tree.
                 - 3 or decide - finds the best branch/node of the tree for each sample.")
    )



wd<-getwd()
packpwd<-("~/in_development/pathPhynder/R")

# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))


# do some operations based on user input
if( opt$step == "all") {
    cat("All steps.\n")
    # system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, 'chrY'))

    system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, 'Y'))

    system(paste("Rscript", paste0(packpwd,"/ancient_SNPs_to_branches.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))

    system(paste("Rscript", paste0(packpwd,"/decide_len.R"),opt$input_tree, 'results_folder', opt$maximumTolerance))

} else if(opt$step == "pileup_and_filter" | opt$step == 1) {
    cat("Running pileup_and_filter\n")
    system(paste("Rscript", paste0(packpwd,"/pileup_and_filter.R"),opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, 'Y'))
} else if(opt$step == "ancient_SNPs_to_branches" | opt$step == 2) {
    cat("Running ancient_SNPs_to_branches\n")
    system(paste("Rscript", paste0(packpwd,"/ancient_SNPs_to_branches.R"),opt$input_tree,opt$prefix,'intree_folder', 'results_folder'))
} else if(opt$step == "decide" | opt$step == 3) {
    cat("Running decide\n")
    system(paste("Rscript", paste0(packpwd,"/decide_len.R"),opt$input_tree, 'results_folder', opt$maximumTolerance))
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