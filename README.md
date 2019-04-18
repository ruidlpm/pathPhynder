# pathPhynder
Title: A workflow for integrating ancient lineages into present-day Y-chromosome phylogenies.

Description: Ancient DNA data is characterized by deamination and low-coverage sequencing, which results in a high fraction of missing data and erroneous calls. These factors affect the estimation of phylogenetic trees with modern and ancient DNA, especially in when dealing with many ancient samples sequenced to lower coverage. Furthermore, most ancient DNA analyses usually rely on known markers on the Y-chromosome, but additional variation will continuously emerge as more data is generated. This workflow offers a solution for integrating ancient and present-day Y-chromososome data, first by identifiying informative Y-chromosome markers in a high coverage dataset, second, by calling and filtering these SNPs in ancient samples and lastly, by traversing the tree and evaluate the number of derived and ancestral markers in the ancients to find the most likely clade where it belongs.





### Background:

 - Ancient DNA characteristics such as low coverage, post-mortem deamination and contamination make accurate Y-chromosome phylogeny estimation very challenging. Furthermore, the majority of ancient DNA Y-chromosome analysis restricts itself to a set of known SNPs used for the identification of Y-chromosome lineages, however, many more variants exist in present-day datasets which can be used to inform aDNA paternal affinities. 
 
 Aims of this method:
  - Provide a tool that allows integrating ancient DNA data from multiple sources (target capture and shotgun sequencing) into present-day phylogenies, which is challenging with current available methods;
  - Use all available Y-chromosomal variants rather than a subset of known SNPs (this provides a greater chance of detecting informative SNPs in low coverage aDNA samples);
  - Provide a visualization tool for ancestral and derived allele state at each branch of the tree, aiding haplogroup determination;
  - Generate a database of ancient DNA Y-chromosome variability which allows observing the main trends of Y-chromosome affinity throughout time and geography.
 


Workflow

![alt text](https://github.com/ruidlpm/Integrating_aDNA_Y/blob/master/figures/workflow_poster.png)


0) Generate an accurate Y-chromosome phylogeny from a vcf file (Use raXML or MEGA, for example, running for several iterations). The quality of the tree has a major impact on SNP assignment to branches and therefore all on downstream analyses. At the moment, pathPhynder does not handle well high numbers of missing genotypes in modern samples, so remove poorly genotyped individuals and SNPs with high missingness across individuals, or try imputing your vcf.



1) Assign informative SNPs to tree branches.

```bash
#will output Rdata file with information about which SNPs map to branches and a bed file for snp calling
Rscript assign_SNPs_noNA.R <input_phylogeny.nwk> <input.vcf> <out prefix>
```


2) Run pathPhynder to call those SNPs in a given dataset of ancient samples and find the best path and branch where these can be mapped in the tree.

```bash
#call positions using samtools mpileup
samtools mpileup <bam> \
--ignore-RG \
--positions <bedfile> \
-f <refgen> > <out.pileup>


#convert to intree format, which contains info about anc der markers
python call_bases_chrY.py \
-i <out.pileup>  \
-m conservative  \
-o <ancient.out.intree.txt>
```




3) Add samples to tree and see results with micro


    ![alt text](https://github.com/ruidlpm/Integrating_aDNA_Y/blob/master/figures/micro_poster.png)



Tutorial:

https://github.com/ruidlpm/Integrating_aDNA_Y/tree/master/tutorial/tutorial_data




TODO:

- Add parameters to continue and stop <PRIORITY>
- Jacknife/Bootstrap procedure
- Update tree as each ancient sample gets added. (Might be worth it to start with aDNA samples sorted by number of variants overlapping the tree)
- can we explore at all rare variation on the Y chromosome?
- HGDP data
- Haplogroup determination - try Yfull tree?
  - IMPUTE missing markers in moderns, for example they may have some R1b defining SNPs, but not all. So if derived for some, impute derived for all. Same for ancestral markers.













Requirements:
 - phytools
 - scales





























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

