#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 16
#SBATCH -t 2-00:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -e results/05_dxy_oct2023/dxy_from_sfs2d_ancTconura_oct2023.err
#SBATCH -J dxy_from_sfs2d_ancTconura_oct2023
#SBATCH -D Rachel/popgen_Tconura

module load bioinfo-tools
module load R_packages
module load gnuparallel

mydir=[redacted]
script_path=$mydir/bin


OUT_PATH=$mydir/popgen_Tconura/results/05_dxy_oct2023
SFS_OUT=$mydir/popgen_Tconura/results/01_sfs_oct2023/Sfs2d/windowed

rm $mydir/popgen_Tconura/scripts/Tconura_dxy_from_2dsfs_jobs.txt

while read line; do
echo "Rscript --vanilla $script_path/dxy_wsfs.R --wsfs $SFS_OUT/$line --npop1 12 --npop2 12 --out $OUT_PATH/${line%.sfs}.dxy.txt" >> $mydir/popgen_Tconura/scripts/Tconura_from_2dsfs_jobs.txt
done < $1

parallel -j 16 < Tconura_from_2dsfs_jobs.txt 
