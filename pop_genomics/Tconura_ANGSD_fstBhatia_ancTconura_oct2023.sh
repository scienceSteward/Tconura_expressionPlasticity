#!/bin/bash
#SBATCH -A naiss2023-22-412
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH --mail-user=rachel.steward@biol.lu.se
#SBATCH --mail-type=ALL
#SBATCH -e results/02_Fst_oct2023/fst_popList_ancTconura_oct2023.err
#SBATCH -J fst_popList_ancTconura_oct2023
#SBATCH -D /proj/snic2021-6-323/Projects/Tconura/working/Rachel/popgen_Tconura

module load bioinfo-tools 

TOOL_PATH=/proj/snic2020-6-222/bin/angsd
IN_PATH=/crex/proj/snic2020-6-222/Projects/Tconura/working/zach/Oct2023redo/results/safs
OUT_PATH=/proj/snic2021-6-323/Projects/Tconura/working/Rachel/popgen_Tconura/results

# /proj/snic2022-6-377/Julio/dec2022redo/results/saf
SFS_OUT=$OUT_PATH/01_sfs_oct2023/Sfs2d
mkdir -p $SFS_OUT

FST_OUT=$OUT_PATH/02_Fst_oct2023
mkdir -p $FST_OUT




while read -r first second; do

pop1="$first"
pop2="$second"

#Prepare files for analysis

# $TOOL_PATH/misc/realSFS $IN_PATH/${pop1}_${saf_suffix}.idx \
#     $IN_PATH/${pop2}_${saf_suffix}.idx \
#     -P 20 -fold 1 > $SFS_OUT/${pop1}.${pop2}_${saf_suffix}.ml

$TOOL_PATH/misc/realSFS fst index $IN_PATH/${pop1}_${saf_suffix}.idx \
    $IN_PATH/${pop2}_${saf_suffix}.idx \
    -sfs $SFS_OUT/${pop1}.${pop2}_${saf_suffix}.ml \
    -fstout $FST_OUT/${pop1}.${pop2}_${saf_suffix} -P 20 -whichFst 1 -fold 1
q
# Calculate FST with a 10kb window
$TOOL_PATH/misc/realSFS fst stats2 $FST_OUT/${pop1}.${pop2}_${saf_suffix}.fst.idx  \
    -win 50000 -step 50000 -type 2 > $FST_OUT/${pop1}.${pop2}_${saf_suffix}.windowed

done < $1

# modified 22 Oct 2023
# SAF = norep-sites
# folded sfs -fold 1
# bhatia estimator -whichFst 1

#Tconura_ANGSD_fstBhatia_ancTconura_oct2023.sh 
