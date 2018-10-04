 
suppressPackageStartupMessages(library("optparse"))


# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")



option_list <- list(
    make_option(c("-i","--input_tree"), default="data/test_tree.nwk",
        help = "Input tree in Newick format. [default \"%default\"]"),
    # make_option(c("-v","--input_vcf"), default="data/test_tree.vcf", 
    #     help = "Input vcf. Needs to be haploid. [default \"%default\"]"),
    make_option(c("-p","--prefix"), default="data/test_tree", 
        help = "Prefix for the data files associated with the tree. [default \"%default\"]"),
    make_option(c("-b","--bam_list"), default="data/bam_list.txt", 
        help = "List with paths to bam files. [default \"%default\"]"),
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
                 - pileup_and_filter - runs pileup in ancient bam files and filters bases.
                 - ancient_SNPs_to_branches - adds anc/der alleles to the input tree.
                 - decide - finds the best branch/node of the tree for each sample.")
    )




# get command line options, if help option encountered print help and exit,
opt <- parse_args(OptionParser(option_list=option_list))

# do some operations based on user input
if( opt$step == "all") {
    cat("All steps.\n")
    system(paste("Rscript ~/in_development/pathPhynder/R/pileup_and_filter.R",opt$bam_list,opt$prefix,'intree_folder', opt$reference, opt$mode, 'chrY'))

    # Rscript pileup_and_filter.R sample_list.txt karmin_tree_prob_snps_excl intree_folder /Users/rm890/software/y_leaf/index_hg19/hg19.fa 'conservative' 'chrY'

    # Rscript ancient_SNPs_to_branches.R  karmin_tree_prob_snps_excl.nwk karmin_tree_prob_snps_excl intree_folder results_folder

    # Rscript decide_len.R  karmin_tree_prob_snps_excl.nwk karmin_tree_prob_snps_excl intree_folder results_folder 3

} else {
    cat(paste(do.call(opt$generator), collapse="\n"))
}
cat("\n")


# bam_list.txt

# Rscript ~/in_development/pathPhynder/pathPhynder.R -s all \
#      -p sole_test \
#      -b bam_list.txt \
#      -i tree_sole_only.nwk \
#      -r /Users/rm890/software/y_leaf/index_hg19/hg19.fa \
#      -m "conservative" \
#      -t 3





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