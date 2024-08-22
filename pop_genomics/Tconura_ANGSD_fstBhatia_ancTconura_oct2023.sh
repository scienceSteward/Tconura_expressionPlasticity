#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 4-00:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -e results/02_Fst_oct2023/fst_popList_ancTconura_oct2023.err
#SBATCH -J fst_popList_ancTconura_oct2023
#SBATCH -D Rachel/popgen_Tconura

module load bioinfo-tools 
WRK_dir=[redacted]
TOOL_PATH=path/to/angsd
IN_PATH=$WRK_dir/zach/Oct2023redo/results/safs
OUT_PATH=$WRK_dir/Rachel/popgen_Tconura/results

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

# Calculate FST with a 10kb window
$TOOL_PATH/misc/realSFS fst stats2 $FST_OUT/${pop1}.${pop2}_${saf_suffix}.fst.idx  \
    -win 50000 -step 50000 -type 2 > $FST_OUT/${pop1}.${pop2}_${saf_suffix}.windowed

done < $1

# modified 22 Oct 2023
# SAF = norep-sites
# folded sfs -fold 1
# bhatia estimator -whichFst 1

# sbatch Tconura_ANGSD_fstBhatia_ancTconura_oct2023.sh inputFile
# input file is list of population pairs
