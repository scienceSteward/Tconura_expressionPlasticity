#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J braker_busco
#SBATCH -e braker_busco.err


module load bioinfo-tools
module load blast_databases
module load BUSCO/5.3.1

annot=$1

# Use UPPMAX's script to setup the Busco environment.
source $BUSCO_SETUP


database=diptera_odb10
run_BUSCO.py -i $annot \
    -o ${annot%.gtf}.${database}.out \
    -l $BUSCO_LINEAGE_SETS/$database \
    -m prot \
    -c 10 -f

database=insecta_odb10
run_BUSCO.py -i $annot \
    -o ${annot%.gtf}.${database}.out \
    -l $BUSCO_LINEAGE_SETS/$database \
    -m prot \
    -c 10 -f 

date

# make sure to have a current augustus config in your path, before running this script
# module load bioinfo-tools augustus
# source $AUGUSTUS_CONFIG_COPY
