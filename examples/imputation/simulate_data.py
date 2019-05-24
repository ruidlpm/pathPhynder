import msprime


tree_sequence = msprime.simulate(
sample_size=10, Ne=1e4, length=1e5, recombination_rate=0,
mutation_rate=2e-8)
#, random_seed=10)


outname=str("simulated_data.vcf")

with open(outname, "w") as vcf_file:
		tree_sequence.write_vcf(vcf_file, ploidy=1, contig_id='Y')
		vcf_file.close()


tree = tree_sequence.first()
print(tree.draw(format="unicode"))


outtreename=str("simulated_data.nwk")


tree = tree_sequence.first()
print(tree.draw(format="unicode"))

with open(outtreename, "w") as tree_file:
	nwk_tree=tree.newick()
	tree_file.write(nwk_tree)


