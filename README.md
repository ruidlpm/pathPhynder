<p align="center">
  <img src="figures/logo.svg" width="180" alt="pathPhynder logo">
</p>

# pathPhynder v1.3.1

[![install with 
bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pathphynder/README.html)

A workflow for ancient DNA placement into reference phylogenies.

Ancient DNA is characterized by post-mortem damage and often low sequencing coverage, resulting in high levels of missing data and increased uncertainty in variant calls. These challenges complicate standard phylogenetic analysis of modern and ancient genomes, particularly when dealing with large numbers of low coverage ancient samples.

For haplogroup determination in haploid systems such as the Y chromosome and mtDNA, most ancient DNA studies currently rely on previously identified phylogenetically informative SNPs. However, as additional high-coverage ancient and modern genomes become available, new variation continues to be discovered, requiring approaches that can incorporate newly identified markers into phylogenetic analyses.

pathPhynder provides a framework for integrating ancient and present-day haploid datasets by: (1) identifying informative markers from high-coverage reference datasets; (2) calling and filtering these variants in ancient samples while accounting for ancient DNA-specific sources of error; and (3) traversing the phylogenetic tree to evaluate the distribution of derived and ancestral states, enabling the assignment of ancient samples to their most likely phylogenetic position.

_________________________________________________

### Documentation

 - [Installation](docs/installation.md) - conda and manual install, plus installing phynder
 - [Workflow overview](docs/workflow.md) - the placement workflow, step by step, and the full command-line reference
 - [Tutorial](tutorial) - quick start with the Y-chromosome reference dataset
 - [Changelog](CHANGELOG.md) - version history

### Citation

Rui Martiniano, Bianca De Sanctis, Pille Hallast, Richard Durbin,
Placing Ancient DNA Sequences into Reference Phylogenies,
Molecular Biology and Evolution, Volume 39, Issue 2, February 2022, msac017
https://doi.org/10.1093/molbev/msac017
