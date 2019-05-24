shopt -s expand_aliases
source ~/.bash_profile

#assign SNPs to branches
pathPhynder -s assign -i small_example_tree.nwk -v small_example.vcf -p small_example
