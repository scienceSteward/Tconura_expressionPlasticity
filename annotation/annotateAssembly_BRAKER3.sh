#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J braker3Annot
#SBATCH -e logs/braker3Annot.err
#SBATCH -o logs/braker3Annot.out

# Modified Nov 5, 2023 by R. Steward
# from script by katharina.hoff@uni-greifswald.de written Jan 12th 2023

mydir=[redacted]
cd $mydir/annotation_Tconura/braker
BRAKER_SIF=./braker3.sif

# Check whether braker3.sif is available

if [[ -z "${BRAKER_SIF}" ]]; then
    echo ""
    echo "Variable BRAKER_SIF is undefined."
    echo "First, build the sif-file with \"singularity build braker3.sif docker://teambraker/braker3:latest\""
    echo ""
    echo "After building, export the BRAKER_SIF environment variable on the host as follows:"
    echo ""
    echo "export BRAKER_SIF=\$PWD/braker3.sif"
    echo ""
    echo "You will have to modify the export statement if braker3.sif does not reside in \$PWD."
    echo ""
    exit 1
fi

# Check whether singularity exists

if ! command -v singularity &> /dev/null
then
    echo "Singularity could not be found."
    echo "On some HPC systems you can load it with \"module load singularity\"."
    echo "If that fails, please install singularity."
    echo "Possibly you misunderstood how to run this script. Before running it, please copy it to the directory where you want to execute it by e.g.:"
    echo "singularity exec -B \$PWD:\$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh ."
    echo "Then execute on the host with \"bash test1.sh\"".
    exit 1
fi

# remove output directory if it already exists

wd=Tcon_braker3_2

if [ -d $wd ]; then
    rm -r $wd
fi

genome_path=$mydir/repeats_Tconura/repeatmasker
genome=pt_042_hifiasm20201214.primary.UPPER.fasta.masked

prots_path=$mydir/annotation_Tconura/braker/prots
prots=Arthropoda.fa

bams_path=$mydir/transcriptomics_Tconura/01_RNA/03_mapped/HISAT2

# remove poorly mapping samples
ls -1 $bams_path/*.bam | grep -v "129\|149\|158\|176" | tr '\n' ',' > P18653.sample.bams 

bams=`cat P18653.sample.bams` 

singularity exec -B ${PWD}:${PWD} \
${BRAKER_SIF} braker.pl \
    --genome=$genome_path/$genome \
    --prot_seq=$prots_path/$prots \
    --bam=$bams \
    # --GENEMARK_PATH=${ETP}/gmes \
    --JAVA_PATH=. \
    --softmasking \
    --workingdir=${wd} \
    --threads 20 


#Tcon_braker3.sh
