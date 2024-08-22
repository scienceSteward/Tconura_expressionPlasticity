# Trimming raw reads 
```
EXPR_home=$WRK_dir/transcriptomics_Tconura
RNA_raw=$EXPR_home/01_RNA/01_raw

RNA_data=[path redacted]
SEQ_group=P18653
SAMPLE_names=exprSampleList_01-2023.txt
THREADS_trim=2

cd $RNA_data/$SEQ_group
ls -1 | grep -v \\. | grep $SEQ_group > $EXPR_home/$SAMPLE_names

cd $EXPR_home
mkdir -p $RNA_raw
```

## Move RNA to 01_RNA/01_raw
```
cat $SAMPLE_names | while read line; 
do cp $RNA_data/$SEQ_group/$line/02-FASTQ/210512_A00187_0484_AH5G5TDSX2/*.fastq.gz ./01_RNA/01_raw/.;
done
```

## FastQC raw reads
```
mkdir -p $RNA_raw/fastQC_raw
rm $EXPR_home/Scripts/fastQC_raw_parallelScript.txt

cd $RNA_raw

module load bioinfo-tools FastQC/0.11.9
ls -1 *.fastq.gz | while read line; 
    do echo "fastqc  $RNA_raw/$line  --outdir $RNA_raw/fastQC_raw" >> $EXPR_home/Scripts/fastQC_raw_parallelScript.txt;
    done
```
Create the sbatch script:
```
cd $EXPR_home 
nano Scripts/fastQC_raw_parallelScript_run.sh
```
Paste: 
```
#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mail-user=[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J fastqc_parallel
#SBATCH -e fastqc_parallel.err
#SBATCH -o fastqc_parallel.out

module load bioinfo-tools FastQC/0.11.9 gnuparallel/20180822 
parallel -j 20 < Scripts/fastQC_raw_parallelScript.txt"
```
Run the script
```
cd $EXPR_home
sbatch Scripts/fastQC_raw_parallelScript_run.sh
```
Run multiqc in an interactive session:
```
interactive -A naiss2023-22-412
module load bioinfo-tools MultiQC/1.10.1
multiqc .
```

## Trimming with TrimGalore 
Prepare the parallel script

```
RNA_trimmed=$EXPR_home/01_RNA/02_trimmed/trimgalore
mkdir -p $RNA_trimmed

rm $EXPR_home/Scripts/trimgalore_parallelScript.txt

module load bioinfo-tools TrimGalore

cd $RNA_raw
ls -1 *R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//g' | while read line; 
do echo "trim_galore -q 10 --phred33 --fastqc --stringency 1 -e 0.1 \
        --length 36 --paired -o $RNA_trimmed --cores 4 \
        $RNA_raw/${line}_R1_001.fastq.gz $RNA_raw/${line}_R2_001.fastq.gz" \
        >> $EXPR_home/Scripts/trimgalore_parallelScript.txt;
done 

cd $EXPR_home
head -1 $EXPR_home/Scripts/trimgalore_parallelScript.txt
```
- `-q/--quality <INT> `: Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimming will be performed first, and adapter trimming is carried in a second round. Other files are quality and adapter trimmed in a single pass. The algorithm is the same as the one used by BWA (Subtract INT from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal). Default Phred score: 20.
- `--stringency <INT>`: Overlap with adapter sequence required to trim a sequence. Defaults to a very stringent setting of 1, i.e. even a single bp of overlapping sequence will be trimmed off from the 3' end of any read.
- `-e <ERROR RATE>`: Maximum allowed error rate (no. of errors divided by the length of the matching region) (default: 0.1)
- `--length <INT>`: Discard reads that became shorter than length INT because of either quality or adapter trimming. A value of '0' effectively disables this behaviour. Default: 20 bp. For paired-end files, both reads of a read-pair need to be longer than <INT> bp to be printed out to validated paired-end files (see option --paired). If only one read became too short there is the possibility of keeping such unpaired single-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.
- `--paired `: This option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files. To pass the validation test, both sequences of a sequence pair are required to have a certain minimum length which is governed by the option `--length` (see above). If only one read passes this length threshold the other read can be rescued (see option --retain_unpaired). Using this option lets you discard too short read pairs without disturbing the sequence-by-sequence order of FastQ files which is required by many aligners. Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g. file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... .

Make the sbatch script:
```
cd $EXPR_home
nano Scripts/trimgalore_parallelScript_run.sh
```
Paste:
```
#!/bin/bash
#SBATCH -A [redacted]
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mail-user=r[redacted]
#SBATCH --mail-type=ALL
#SBATCH -J trimgalore_parallel
#SBATCH -e trimgalore_parallel.err
#SBATCH -o trimgalore_parallel.out

module load bioinfo-tools TrimGalore gnuparallel/20180822 

parallel -j 5 < Scripts/trimgalore_parallelScript.txt"
```
Run the script
```
sbatch Scripts/trimgalore_parallelScript_run.sh
```
Run interactive session:
```
cd $RNA_trimmed
interactive -A snic2022-5-25
module load bioinfo-tools MultiQC/1.10.1
multiqc .
```
