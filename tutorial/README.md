# Tutorial

In this tutorial, we are going to identify informative SNPs in a collection of 2014 Y-chromosome samples and ~120,000 SNPs and then add ancient samples into the tree.

The data required for this dataset is in the 'data/BigTree_Y/' folder

## pathPhynder best path method

```bash
#the data required for this dataset is in the 'data/BigTree_Y/' folder
# change the path to files according your present directory
# 1) Assign SNPs to branches of the tree 
phynder -B -o branches.snp $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk $PATHPHYNDER_DATA/BigTree_Y/BigTree.Y.201219.vcf.gz

# 2) Prepare sites (writes bed files for variant calling and other files for phylogenetic placement).
#The -G parameter is optional and in this case adds ISOGG haplogroup information to each variant.
pathPhynder -s prepare -i $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk -p BigTree_Y_data -f branches.snp -G $PATHPHYNDER_DATA/210513.snps_isogg_curated.txt 

# 3) Run pathPhynder best path, call variants, place samples, plot results (the -G can be used to identify haplogroups and it is optional)
pathPhynder  -i $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk   -p tree_data/BigTree_Y_data -l bam.list -s all -t 100 -G $PATHPHYNDER_DATA/210513.snps_isogg_curated.txt 

NOTE: If your bam files are aligned to GRCH38 or if you are working with non-human genomes, then you must pass the appropriate reference genome using the -r option.

```


## phynder maximum likelihood method

```bash
#the data required for this dataset is in the 'data/BigTree_Y/' folder

# 1) Assign SNPs to branches of the tree, as above. Skip if you have done this before.
phynder -B -o branches.snp $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk $PATHPHYNDER_DATA/BigTree_Y/BigTree.Y.201219.vcf.gz

# 2) Prepare sites (writes bed files for variant calling and other files for phylogenetic placement). Skip if you have done this before.
#The -G parameter is optional and in this case adds ISOGG haplogroup information to each variant.
pathPhynder -s prepare -i $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk -p BigTree_Y_data -f branches.snp -G $PATHPHYNDER_DATA/200803.snps_isogg.txt

#convert calls to vcf
make_vcf.R intree_folder/ chrY ancient_calls.vcf

#place samples with phynder
phynder -q ancient_calls.vcf -p 0.01 -o query.phy $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk $PATHPHYNDER_DATA/BigTree_Y/BigTree.Y.201219.vcf.gz

#plot results
plot_likes.R $PATHPHYNDER_DATA/BigTree_Y/bigtree_annotated_V1.nwk query.phy results_folder

```


##Test data
Adding ancient samples to the tree - example with 27 samples from Haak et al. (2015) You can get them from this link https://drive.google.com/open?id=1utVO5pua6dk8eUtW3CiEup5ouomwN6QT


##Description of the output files.

```
<sample>.pileup.txt - pileup for a given sample at branch defining SNPs
(see this link for description of thsi format http://samtools.sourceforge.net/pileup.shtml#)
Y       2650701 G       2       ..      JJ
Y       2657214 G       1       ,       J
Y       2657247 G       1       ,       J

<sample>.intree.txt - status at each marker after filtering
position,ref base,alt base,ref_count,alt_count,genotype
23973594 G T 0 2 1
23980238 C A 1 0 0
23987274 C T 0 1 -9
(Note: -9 indicates missing data, or that the genotype was filtered out)

<sample>.best_path.pdf - contains tree and best path for a given sample. Red and green circles indicate ancestral and derived alleles, respectivelly.

<sample>.best_path_report.txt - counts of markers supporting and in conflict with query/aDNA sample membership to branches along the best path.

<sample>.all_paths_report.txt - ccounts of markers supporting and in conflict with query/aDNA sample membership to branches along all paths from root to tips.

<sample>.best_node_info.txt - node and position along branch where the ancient sample belongs

<sample>.hg_in_tree_status.txt - counts of markers at haplogroup defining branches (markers in the tree data only)

<sample>.hg_in_tree_status_derived_only.txt - counts of derived markers at haplogroup defining branches (markers in the tree data only)

<sample>.posteriors.pdf - results of sample placement using phynder's maximum likelihood method.

final_tree.nwk - newick tree with ancient samples added into it.

final_tree.pdf - figure of the tree with ancient samples added into it.

```


Additional comments:

This is an useful tool for exploration of Y-chromosome and other haploid data and ancient DNA sample affinity to present-day lineages.
I recommend exploring the different parameters, especially the stringency of the filtering ('default', 'no-filter' or  'transversions'),
the number of ancestral markers tolerated at a branch before switching to a different path, and the pileup mismatch threshold.
These can give slightly different results, and different parameters will definitely work better for some samples than for others,
depending on the nature of the data (shotgun sequencing vs capture, UDG/USER treatment vs. untreated, coverage, etc...) 

