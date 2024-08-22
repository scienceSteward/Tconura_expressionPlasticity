#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J decorate_braker
#SBATCH -e decorate_braker.err
#SBATCH -o decorate_braker.out

transcripts_merged=transcripts_merged.gff
braker=Tcon_braker3_1/braker.adj.gtf
braker_with_utrs=Tcon_braker3_1/braker.adj.UTRs.gtf

python3 ./intervalTree/stringtie2utr.py -g $braker -s $transcripts_merged -o $braker_with_utrs

# decorate_braker.sh
