#!/bin/bash
#SBATCH -A naiss2023-5-41
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=rachel.steward@biol.lu.se
#SBATCH --mail-type=ALL
#SBATCH -e results/07_theta_oct2023/theta_popList_ancTconura_oct2023.err
#SBATCH -J theta_popList_ancTconura_oct2023
#SBATCH -D /proj/snic2021-6-323/Projects/Tconura/working/Rachel/popgen_Tconura

module load bioinfo-tools 

TOOL_PATH=/proj/snic2020-6-222/bin/angsd
IN_PATH=/crex/proj/snic2020-6-222/Projects/Tconura/working/zach/Oct2023redo/results/safs
OUT_PATH=/proj/snic2021-6-323/Projects/Tconura/working/Rachel/popgen_Tconura/results

SFS_OUT=$OUT_PATH/01_sfs_oct2023/Sfs1d
mkdir -p $SFS_OUT

THETA_OUT=$OUT_PATH/07_theta_oct2023
mkdir -p $THETA_OUT

saf_suffix=norep-sites.saf

while read -r first; do

pop1="$first"

#Prepare files for analysis

$TOOL_PATH/misc/realSFS $IN_PATH/${pop1}_${saf_suffix}.idx \
   -P 20 -fold 1 > $SFS_OUT/${pop1}_${saf_suffix}.sfs

# calculate per-site theta
$TOOL_PATH/misc/realSFS saf2theta $IN_PATH/${pop1}_${saf_suffix}.idx \
    -sfs $SFS_OUT/${pop1}_${saf_suffix}.sfs \
    -outname $THETA_OUT/${pop1} -P 20 -fold 1

# Calculate theta with a 50kb window
$TOOL_PATH/misc/thetaStat do_stat $THETA_OUT/${pop1}.thetas.idx \
    -win 50000 -step 10000 -outnames $THETA_OUT/${pop1}_theta.50_10.gz

done < $1

# Tconura_ANGSD_theta_ancTconura_oct2023.sh

