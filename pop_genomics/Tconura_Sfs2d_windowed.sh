#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -e results/01_sfs_oct2023/Sfs2d/windowedSFS_popList_ancTconura_oct2023.err
#SBATCH -J windowedSFS_popList_ancTconura_oct2023
#SBATCH -D Rachel/popgen_Tconura


module load bioinfo-tools gnuparallel

# path to angsd
TOOL_PATH=path/to/angsd
export TOOL_PATH

# path to safs
IN_PATH=$WRK_dir/zach/Oct2023redo/results/safs
export IN_PATH

saf_suffix=norep-sites.saf
export saf_suffix

# path to results directory
OUT_PATH=$WRK_dir/Rachel/popgen_Tconura/results
SFS_OUT=$OUT_PATH/01_sfs_oct2023/Sfs2d/windowed
mkdir -p $SFS_OUT

#Prepare files for analysis
# bedtools makewindows -g pt_042_hifiasm20201214.primary.fasta.fai -w 50000 > pt_042_hifiasm20201214.primary.50kb.50kb.bed
# awk ' {printf "%s:%s-%s\n", $1, $2+1, $3 } ' pt_042_hifiasm20201214.primary.50kb.50kb.bed > $SFS_OUT/50kb_regions.txt
saf_suffix="norep-sites.saf" # Assuming this is defined somewhere in your script

# Precompute the output file name to avoid recomputation in the loop

# Define a function to process each line
process_line() {
    echo "$(echo $1 | tr ":" " " | tr "-" " " ) $($TOOL_PATH/misc/realSFS "$IN_PATH/${pop1}_${saf_suffix}.idx" "$IN_PATH/${pop2}_${saf_suffix}.idx" -r "$1" -fold 1)"
}
export -f process_line


while read -r first second; do

pop1="$first"
pop2="$second"

export pop1
export pop2

output_file="${pop1}.${pop2}_${saf_suffix}.windowed.sfs"
export output_file

# Use GNU Parallel to run the process_line function on each line of the input file
# -j 20 specifies the number of cores to use
# --keep-order ensures that the output is in the same order as the input

parallel -j 20 --keep-order process_line ::: $(cat "$SFS_OUT/50kb_regions.txt") >> "$SFS_OUT/$output_file"

done < $1

# where $1 is a list of tab separated comparisons
#scripts/Tconura_Sfs2d_windowed.sh
