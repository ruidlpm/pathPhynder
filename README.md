# pathPhynder
Title: A workflow for integrating ancient lineages into present-day Y-chromosome phylogenies.

Description: Ancient DNA data is characterized by deamination and low-coverage sequencing, which results in a high fraction of missing data and erroneous calls. These factors affect the estimation of phylogenetic trees with modern and ancient DNA, especially in when dealing with many ancient samples sequenced to lower coverage. Furthermore, most ancient DNA analyses usually rely on known markers on the Y-chromosome, but additional variation will continuously emerge as more data is generated. This workflow offers a solution for integrating ancient and present-day Y-chromososome data, first by identifiying informative Y-chromosome markers in a high coverage dataset, second, by calling and filtering these SNPs in ancient samples and lastly, by traversing the tree and evaluate the number of derived and ancestral markers in the ancients to find the most likely clade where it belongs.


Aims of this method:
  - Provide a tool that allows integrating ancient DNA data from multiple sources (target capture and shotgun sequencing) into present-day phylogenies, which is challenging with current available methods;
  - Use all available Y-chromosomal variants rather than a subset of known SNPs (this provides a greater chance of detecting informative SNPs in low coverage aDNA samples);
  - Provide a visualization tool for ancestral and derived allele state at each branch of the tree, aiding haplogroup determination;
  - Generate a database of ancient DNA Y-chromosome variability which allows observing the main trends of Y-chromosome affinity throughout time and geography.

_________________________________________________

### Installation

1) Download pathPhynder to your computer and install the following packages:
 - phytools (example:conda install -c bioconda r-phytools)
 - scales
 - optparse

2) For running pathPhynder, you can either specify the whole path:

```
Rscript ~/<path_to_pathPhynder_folder>/pathPhynder.R
```

or add the following line to your ~/.bash_profile. Replace <path_to_pathPhynder_folder> with the location of the pathPhynder folder in your system.
```
alias pathPhynder="Rscript <path_to_pathPhynder_folder>"
```
For example, if you have downloaded the folder to your ~/software/ directory, then you would add the following lines to ~/.bash_profile.
```
alias pathPhynder="Rscript ~/software/pathPhynder/pathPhynder.R"
```
and then:
```
source ~/.bash_profile
```

3) Test if it works.

```
pathPhynder -h
```
_________________________________________________

### Workflow

0) Generate an accurate Y-chromosome phylogeny from a vcf file (Use raXML or MEGA, for example, running for several iterations). The quality of the tree has a major impact on SNP assignment to branches and therefore all on downstream analyses. At the moment, pathPhynder does not handle well high numbers of missing genotypes in modern samples, so remove poorly genotyped individuals and SNPs with high missingness across individuals, or try imputing your vcf.

1) Assign informative SNPs to tree branches.

```bash
#will output tables with information about which SNPs map to each branch of the tree and a bed file for snp calling.
pathPhynder -s assign -i <tree>.nwk -v <tree_data>.vcf -p <prefix_output>
```


2) Run pathPhynder to call those SNPs in a given dataset of ancient samples and find the best path and branch where these can be mapped in the tree.

```bash
#To run all the steps at once, including variant calling, choosing the best path and adding ancient samples to the tree
pathPhynder -s all -i <tree>.nwk -p path_to/<prefix_output> -l <sample.list>

#(Optional) Note that you can also run the analysis on a single bam file
pathPhynder -s all -i <tree>.nwk -p path_to/<prefix_output> -b <sample.bam>
```

If you want to run each step individually:
```
pathPhynder -s <1 or pileup_and_filter> -i <tree>.nwk -p path_to/<prefix_output> -l <sample.list>

pathPhynder -s <2 or chooseBestPath> -i <tree>.nwk -p path_to/<prefix_output> -l <sample.list>

pathPhynder -s <3 or addAncToTree> -i <tree>.nwk -p path_to/<prefix_output> -l <sample.list>
```

There are also a few additional parameters that can be adjusted according to the user's needs.
pathPhynder -h.
```
Options:
	-s STEP, --step=STEP
		Specifies which step to run. Options:
    			- assign - assigns SNPs to branches.
    			- all - runs all steps to map ancient SNPs to branches (1,2,3).
    			- 1 or pileup_and_filter - runs pileup in ancient bam files and filters bases.
    			- 2 or chooseBestPath - finds the best branch/node of the tree for each sample.
    			- 3 or addAncToTree - adds ancients samples to tree.
    			[default all]

	-i INPUT_TREE, --input_tree=INPUT_TREE
		Input tree in Newick format. [required]

	-v INPUT_VCF, --input_vcf=INPUT_VCF
		Input vcf. Only needed for SNP to branch assignment. Needs to be haploid.

	-p PREFIX, --prefix=PREFIX
		Prefix for the data files associated with the tree.
        	These were previously generated in the branch assignment step. [required]

	-b BAM_FILE, --bam_file=BAM_FILE
		Input bam file. [required]

	-l LIST_OF_BAM_FILES, --list_of_bam_files=LIST_OF_BAM_FILES
		List of paths to bam files. [required]

	-r REFERENCE, --reference=REFERENCE
		Reference genome (fasta format). [default "~/in_development/pathPhynder/data/reference_sequences/hs37d5_Y.fa.gz"]

	-m MODE, --mode=MODE
		Mode for filtering pileups.  Options: relaxed or conservative. [default "conservative"]

	-t MAXIMUMTOLERANCE, --maximumTolerance=MAXIMUMTOLERANCE
		Maximum number of ALT alleles tolerated while traversing the tree.
                If exceeded, the algorithm stops and switches to the next path. [default 3]

	-c PILEUP_READ_MISMATCH_THRESHOLD, --pileup_read_mismatch_threshold=PILEUP_READ_MISMATCH_THRESHOLD
		Mismatch threshold for accepting a variant (for cases where reads for both alleles are present in pileup).
        	For a variant to pass filtering, reads containing the most frequent allele have to occur at least
        	at x proportion of the total reads. 1 is the most stringent, 0.5 is the most relaxed. [default 0.7]

	-o OUTPUT_PREFIX, --output_prefix=OUTPUT_PREFIX
		Sample name. This only works if a single bam file is used as an input. [default bamFileName]

	-h, --help
		Show this help message and exit
```


Tutorial:

https://github.com/ruidlpm/pathPhynder/tree/master/tutorial



________________________________________________________________________


To do:

 - [ ] fix header in output vcf after imputation
- [ ] SNP table produced by assign_noNA.R should contain Anc/Derived information
- [ ] need to have a path with gr37 hg19 and gr38 fa and respective indexes
- [ ] need to have an independent snps.txt file, which does not rely on poznik. Ideally for both ISOGG2016 and 18
- [ ] plotting requires a list of hgs of ancient samples and modern samples. This should be an optional requirement.
- [ ] The main clades could be added to the plots.
- [ ] need to add warnings all over the code, if the files are not present, if they do not conform with the specs, etc...
- [ ] could output a summary for each sample:
 		Read 13044 positions from file M3397.chrY.realigned.calmd.bam.intree.txt
 		SNP count:  381 derived and 11847 ancestral.
- [ ] Ideally, I shouldn't rely on a previous haplogroup assignment, I should produce one myself... Supported by the graph produced.
- [ ] Instead of neighbour joining, it would be good to have properly calibrated trees. Could use the counts at each branch to estimate coalescence.
- [ ] My ancient Y chr capture samples could be really useful for improving calibration times.


- Jacknife/Bootstrap procedure
- Update tree as each ancient sample gets added. (Might be worth it to start with aDNA samples sorted by number of variants overlappin$
- can we explore at all rare variation on the Y chromosome?
- Haplogroup determination - try Yfull tree?
