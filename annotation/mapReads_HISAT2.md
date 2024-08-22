
# Mapping trimmed reads for braker3 using hisat2

## Environment
```
EXPR_home=$WRK_dir/transcriptomics_Tconura

RNA_trimmed=$EXPR_home/01_RNA/02_trimmed/trimgalore
RNA_mapped=$EXPR_home/01_RNA/03_mapped/HISAT2

REF_dir=$WRK_dir/refs_Tconura
REF=pt_042_hifiasm20201214.primary.fasta
ANNOT_dir=$WRK_dir/annotation_Tconura
ANNOT=Tconura.all.maker.features.gffreadConvert.gtf
```


## Hisat2 genome index
Create the genome index in the $EXPR_home folder.

```
cd $EXPR_home
```

### Generate splice site and exon lists, and build index. 
Within an interactive session, set up the environmental variables:
```
REF_dir=$WRK_dir/refs_Tconura
REF=pt_042_hifiasm20201214.primary.fasta
ANNOT_dir=$WRK_dir/annotation_Tconura
ANNOT=Tconura.all.maker.features.gffreadConvert.gtf
module load bioinfo-tools HISAT2/2.2.1
```
#### Extract splice sites:
```
extract_splice_sites.py $ANNOT_dir/$ANNOT > ${ANNOT%.gtf}_splicesites.tsv
```
#### Extract exons: 
```
hisat2_extract_exons.py $ANNOT_dir/$ANNOT > ${ANNOT%.gtf}_exons.tsv
```

### Build Hisat2 index
```
#!/bin/sh
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J hisat_index
#SBATCH -e hisat_index.err
#SBATCH -o hisat_index.out

module load bioinfo-tools HISAT2/2.2.1

#Set environmental variables
WRK_dir=[redacted]
REF_dir=$WRK_dir/refs_Tconura
REF=pt_042_hifiasm20201214.primary.fasta
ANNOT_dir=$WRK_dir/annotation_Tconura
ANNOT=Tconura.all.maker.features.gffreadConvert.gtf

# Set threads 
THREADS=20 

# build index (do not run!)
hisat2-build -p $THREADS --ss ${ANNOT%.gtf}_splicesites.tsv\
  --exon ${ANNOT%.gtf}_exons.tsv $REF_dir/$REF ${REF} 2> \
  ${REF%.fa}_${ANNOT%.gtf}_hisat-build.log

date
```

## Hisat2 mapping: parallel jobs

### First pass 
Create text file of first-pass mapping jobs
```
cd $RNA_trimmed
rm $EXPR_home/Scripts/hisat2_parallelScript.txt

THREADS=20

for reads in $RNA_trimmed/*_R1_001_val_1.fq.gz; 
  do
    file=$(basename $reads R1_001_val_1.fq.gz);
    in1=$RNA_trimmed/${file}R1_001_val_1.fq.gz; 
    in2=$RNA_trimmed/${file}R2_001_val_2.fq.gz;

    # align reads with hisat2
   echo "hisat2 -p $THREADS -x $EXPR_home/${REF} -1 $in1 -2 $in2 -S $RNA_mapped/${file}Aligned_out.sam --summary-file $RNA_mapped/${file}summaryfile.txt --dta --rna-strandness FR" >> $EXPR_home/Scripts/hisat2_parallelScript.txt;
  done

head -1 $EXPR_home/Scripts/hisat2_parallelScript.txt
```
  
## Hisat2 mapping: batch script
```
cd $EXPR_home
nano Scripts/hisat2_parallelScript_run.sh
```

```
#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J STAR_mapping
#SBATCH -e STAR_mapping.err
#SBATCH -o STAR_mapping.out

# Load modules
module load bioinfo-tools HISAT2 gnuparallel samtools

# Set environmental variables
WRK_dir=[redacted]
EXPR_home=$WRK_dir/transcriptomics_Tconura
RNA_mapped=$EXPR_home/01_RNA/03_mapped/HISAT2

# Begin Hisat2 first pass
echo "Hisat mapping in progress" 

# Move into 1st_SJ directory

parallel -j1 < $EXPR_home/Scripts/hisat2_parallelScript.txt

echo "hisat mapping finished"

cd $RNA_mapped
for reads in *.sam; 
  do
    file=$(basename $reads .sam);
    
    samtools view -@ $THREADS -bS ${file}.sam | \
      samtools sort -@ $THREADS - -o ${file}_sorted.bam; 
    
    samtools index -@ $THREADS ${file}_sorted.bam;
    
    echo "${file} bam sorted, indexed"; 
  done

```
