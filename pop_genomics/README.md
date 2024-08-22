We used site allele frequency files to calculate Fst and dxy between the CH and CO source populations (CHST, COSK) and to estimate nucleotide diversity (π) and Tajima’s D over 50kb non-overlapping windows throughout the 1.9G T. conura genome. The two-dimensional folded site frequency spectrum (2DSFS) was used to calculate Fst using the Bhatia estimator (`Tconura_ANGSD_fstBhatia_ancTconura_oct2023.sh` for all population pairs). 

To calculate dxy, we first re-calculated the 2DSFS for every 50kb window, then calculated dxy using a modified version of dxy_wsfs.py script from D. Marques (https://github.com/marqueda/PopGenCode/blob/master/dxy_wsfs.py, accessed Nov. 2023) modified to run in R (dxy_wsfs_windowed.R). 

Rscript --vanilla /proj/snic2021-6-323/Projects/Tconura/working/Rachel/bin/dxy_wsfs.R --wsfs /proj/snic2021-6-323/Projects/Tconura/working/Rachel/popgen_Tconura/results/05_dxy_oct2023/01_sfs_oct2023/Sfs2d/windowed/CHES.CHFI_norep-sites.saf.windowed.sfs --npop1 12 --npop2 12 --out /proj/snic2021-6-323/Projects/Tconura/working/Rachel/popgen_Tconura/results/05_dxy_oct2023/CHES.CHFI_norep-sites.saf.windowed.dxy.txt



We also used ANGSD to calculate π and Tajima’s D over 50kb windows.
