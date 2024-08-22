# process braker
| Annotation        | Details |
| ------------- |:-------------:| 
| braker.gtf      | default braker output | 
| braker.adj.gtf      | braker output with tsebra intron support at 0.25, prefix 'Tcon' |
| braker.adj.UTR.gtf | UTRs added with intervaltree/stringtie2utr.py      |  

## Get statistics
agat_sp_statistics.pl --gff braker.gtf --output braker.statistics.txt
agat_sp_statistics.pl --gff braker.adj.gtf --output braker.adj.statistics.txt
agat_sp_statistics.pl --gff braker.adj.UTR.gtf --output braker.adj.UTR.statistics.txt

## Run busco
sbatch braker3_busco.sh braker.adj.aa

## Modify for downstream analyses 
awk '($3 != "gene" && $3 != "transcript" )'  braker.adj.UTR.gtf > braker.adj.UTR.mod.gtf 
awk '($3 != "gene" && $3 != "transcript" )'  braker.adj.gtf > braker.adj.mod.gtf 

## Keep longest isoform
agat_sp_keep_longest_isoform.pl --gff braker.adj.UTR.mod.gtf --output braker.adj.UTR.mod.longIso.gff

## extract amino acid sequences for longest isoform
genome_path=$mydir/repeats_Tconura/repeatmasker
genome=pt_042_hifiasm20201214.primary.UPPER.fasta.softmasked

gff_file=braker.adj.UTR.mod.longIso.gff
cds_outfile=braker.adj.UTR.mod.longIso.codingseq
prot_outfile=braker.adj.UTR.mod.longIso.aa
gffread "$gff_file" -g "$genome_path/$genome" -J -w "$cds_outfile" -y "$prot_outfile"

cat $prot_outfile | sed 's/.*gene=/>/' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > ${prot_outfile%.aa}_2line.aa


