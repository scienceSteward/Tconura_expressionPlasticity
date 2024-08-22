We used site allele frequency files to calculate Fst and dxy between the CH and CO source populations (CHST, COSK) and to estimate nucleotide diversity (π) and Tajima’s D over 50kb non-overlapping windows throughout the 1.9G *T. conura* genome. The two-dimensional folded site frequency spectrum (2DSFS) was used to calculate Fst using the Bhatia estimator (`Tconura_ANGSD_fstBhatia_ancTconura_oct2023.sh` for all population pairs). 

To calculate dxy, we first re-calculated the 2DSFS for every 50kb window (`TconuraSfs2d_windowed.sh`), then calculated dxy using a modified version of dxy_wsfs.py script from D. Marques (https://github.com/marqueda/PopGenCode/blob/master/dxy_wsfs.py, accessed Nov. 2023) modified to run in R (`dxy_wsfs_windowed.R`) which was run in the command line (`Tconura_dxy_from_2dsfs_windowed.sh`). 

We also used ANGSD to calculate π and Tajima’s D over 50kb windows (`Tconura_ANGSD_theta_ancTconura_oct2023.sh`).
