# Map to genome-derived transcriptome

```
WRK_dir=[redacted]
EXPR_home=$WRK_dir/transcriptomics_Tconura

RNA_trimmed=$EXPR_home/01_RNA/02_trimmed/trimgalore
RNA_mapped=$EXPR_home/01_RNA/03_mapped

REF_dir=$WRK_dir/refs_Tconura
REF=pt_042_hifiasm20201214.primary.fasta

ANNOT_dir=$WRK_dir/annotation_Tconura/braker/Tcon_braker3_1
ANNOT=braker.adj.UTR.mod.gtf	
TRANSCRIPTOME=braker.adj.UTR.mod
```

Create the transcriptome from the annotation + genome. Make sure to set the environment within the interactive session. 
```
interactive -A naiss2023-5-41 

module load bioinfo-tools gffread

# gffread $ANNOT_dir/$ANNOT -F -w $ANNOT_dir/${TRANSCRIPTOME}.aa -y $ANNOT_dir/${TRANSCRIPTOME}.codingseq  -g $REF_dir/$REF
```

## Decoys
From https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/. 

Generate list of decoy sequences:
```

mkdir $REF_dir/transcriptome
cd $REF_dir/transcriptome

REF_masked=$WRK_dir/repeats_Tconura/repeatmasker/${REF%.fasta}.UPPER.fasta.hardmasked
cat $REF_masked | grep "^>" | cut -d " " -f 1 > decoys.txt
head decoys.txt
```
Concatenate transcriptome and genome
```
TRANSCRIPTOME=braker.adj.codingseq
cat $ANNOT_dir/$TRANSCRIPTOME $REF_dir/$REF > Tconura.cat.braker.transcriptome.genome.fa
```
# Salmon Indexing

Create and run the salmon indexing script
```
nano salmon_indexing.sh
```

```
#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J salmon_indexing
#SBATCH -e salmon_indexing.err
#SBATCH -o salmon_indexing.out

module load bioinfo-tools Salmon

salmon index -t Tconura.cat.transcriptome.genome.fa -d decoys.txt -p 20 -i salmon_index --gencode 
```

# Salmon alignment + quantification
From https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon
```
SALMON_dir=$EXPR_home/01_RNA/07_Salmon
mkdir -p $SALMON_dir

cd $RNA_trimmed
rm $EXPR_home/Scripts/salmon_map_quant_ISF_parallelScript.txt

for reads in $RNA_trimmed/*_R1_001_val_1.fq.gz; 
  do
    file=$(basename $reads R1_001_val_1.fq.gz);
    
    index=$REF_dir/transcriptome/salmon_index
    in1=$RNA_trimmed/${file}R1_001_val_1.fq.gz; 
    in2=$RNA_trimmed/${file}R2_001_val_2.fq.gz;
    out=$SALMON_dir/${file}salmon_quant_ISF;

    echo "salmon quant -i $index -l ISF -1 $in1 -2 $in2 --validateMappings -o $out" -p 20 --gcBias >> $EXPR_home/Scripts/salmon_map_quant_ISF_parallelScript.txt;
  done
head -1 $EXPR_home/Scripts/salmon_map_quant_ISF_parallelScript.txt
```
Create and run batch script
```
nano $EXPR_home/Scripts/salmon_map_quant_ISF_parallelScript_run.sh
```

```
#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J salmon_quant
#SBATCH -e Logs/salmon_quant_ISF.err
#SBATCH -o Logs/salmon_quant_ISF.out

module load bioinfo-tools Salmon gnuparallel
WRK_dir=[redacted]
EXPR_home=$WRK_dir/transcriptomics_Tconura
cd $EXPR_home
parallel -j1 < `head -1 $EXPR_home/Scripts/salmon_map_quant_ISF_parallelScript.txt`

```
Run batch script
```
cd $EXPR_home
sbatch $EXPR_home/Scripts/salmon_map_quant_ISF_parallelScript_run.sh
```


