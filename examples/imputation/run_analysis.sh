
python simulate_data.py

Rscript fix_tree_labels.R 

pathPhynder -s assign -i simulated_data.nwk -v simulated_data.vcf -p test

Rscript plot_tree.R

open *pdf

Rscript simulate_missing.R


Rscript ../../R/Ympute.R simulated_data.nwk simulated_data.miss0.1.vcf simulated_data.miss0.1.imputed.vcf

