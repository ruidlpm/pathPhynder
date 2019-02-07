# pathPhynder
Title: A workflow for integrating ancient lineages into present-day Y-chromosome phylogenies.

Description: Ancient DNA data is characterized by deamination and low-coverage sequencing, which results in a high fraction of missing data and erroneous calls. These factors affect the estimation of phylogenetic trees with modern and ancient DNA, especially in when dealing with many ancient samples sequenced to lower coverage. Furthermore, most ancient DNA analyses usually rely on known markers on the Y-chromosome, but additional variation will continuously emerge as more data is generated. This workflow offers a solution for integrating ancient and present-day Y-chromososome data, first by identifiying informative Y-chromosome markers in a high coverage dataset, second, by calling and filtering these SNPs in ancient samples and lastly, by traversing the tree and evaluate the number of derived and ancestral markers in the ancients to find the most likely clade where it belongs.


________________________________________________________________________


To do:

Urgent
 - [ ] currently the python program which filters alleles sets the position as missing in case of mismatches. I'm losing ~7% of variants because of this. It should be straighforward to recover them! 
 
 Basically if p_der<0.5, assign anc, if p_der>0.5, assign der. If p_der==0.5, then leave as missing.
________


 - [ ] fix header in output vcf after imputation
- [ ] SNP table produced by assign_noNA.R should contain Anc/Derived information
- [ ] assign_noNA.R should automatically produce "siteschr" files
- [ ] need to have a path with gr37 hg19 and gr38 fa and respective indexes
- [ ] need to have an independent snps.txt file, which does not rely on poznik. Ideally for both ISOGG2016 and 18
- [ ] the plotting is very fragile. Need a function which can be reused.
- [ ] pathPhynder can deal with 1000 Genomes type data. But the HGDP has 60 Gb which is not possible at the - [ ]ment. Should implement something to split vcf files and run assign to branches
- [ ] plotting requires a list of hgs of ancient samples and modern samples. This should be an optional requirement.
- [ ] The main clades could be added to the plots.
- [ ] Fix also a problem with sample names and anchgs file.
			[1] "TV32032extra.chrY.realigned.calmd.bam___NA___NA"
- [ ] need to add warnings all over the code, if the files are not present, if they do not conform with the specs, etc...
- [ ] formula to estimate the size of the plot automattically
- [ ] sample list only works for first step. Need to fix this to allow processing only one sample or a subset
- [ ] could output a summary for each sample:
 		Read 13044 positions from file M3397.chrY.realigned.calmd.bam.intree.txt
 		SNP count:  381 derived and 11847 ancestral.
- [ ] output a table with best node for sample, this would allow plotting only a chose subset of ancient individuals.
- [ ] Ideally, I shouldn't rely on a previous haplogroup assignment, I should produce one myself... Supported by the graph produced.
- [ ] plotting needs to highlight the chosen path, and not only the branches containing derived alleles.
- [ ] also input could be a single bam or a list of bams
- [ ] Instead of neighbour joining, it would be good to have properly calibrated trees. Could use the counts at each branch to estimate coalescence.
- [ ] My ancient Y chr capture samples could be really useful for improving calibration times.

