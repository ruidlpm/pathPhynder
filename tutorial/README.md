# Tutorial

In this tutorial, we are going to identify informative SNPs in the Karmin et al., 2015 Y-chromosome dataset and then add ancient samples into the tree.

You can find the original paper here
https://genome.cshlp.org/content/25/4/459.full


step 1 - Assign SNPs to branches of the tree

```
#download the Karmin VCF file
wget http://evolbio.ut.ee/chrY/PLINKandVCF/Karmin_et_al_2015_chrY_454samples_CGplusIlluminaA00_filtered1_noAustralians.vcf


#assign SNPs to branches of the tree
time pathPhynder -s assign -i karmin.nwk -v Karmin_et_al_2015_chrY_454samples_CGplusIlluminaA00_filtered1_noAustralians.vcf -p karmin_data



	table with markers assigned to each branch written to tree_data/karmin_data.edge_df.txt)

	31239 informative positions (auxiliary information written to tree_data/karmin_data.sites.txt)

	31239 informative positions for variant calling (written to tree_data/karmin_data.sites.bed)

	6834 positions were not added (report written to tree_data/karmin_data.report_not_added_SNPs.txt)



real	5m46.280s
user	4m54.475s
sys	0m50.063s


Description of the output files of step 1.

```
karmin_data.edges_df.txt - contains information about SNPs found at each branch.


karmin_data.sites.bed - contains the coordinates of branch defining varaints which will be called in the ancient samples
Y	2650700	2650701
Y	2650852	2650853
Y	2651026	2651027
Y	2651324	2651325
...

karmin_data.sites.txt - contains information about each informative variant, including if they define any known haplogroup and whether the derived allele is ALT (+) of REF (-). It also states to which branch they have been assigned to (706, 707, 901, etc).
Y	U1	G2a2b2a1a1a	15058878	G->A	G	A	706	-
Y	Z6390	G2a2b2a1a1a2	2830607	C->A	C	A	707	+
Y	NA	NA	28715314	C->T	C	T	901	+
...

karmin_data.report_not_added_SNPs.txt
2650414	monomorphic snp
2657176	incompatible with tree
2661371	missing data
...
```

step 2 - Adding ancient samples to the tree - example with 27 samples from Haak et al. (2015)

```
pathPhynder -s all -i karmin.nwk -p tree_data/karmin_data -l sample.list


real	2m44.761s
user	2m33.428s
sys	0m6.820s
```

