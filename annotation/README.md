# Annotation workflow

1. Trimmed reads were mapped with HISAT2: `mapReads_HISAT2.md`
2. A new TSEBRA config file (`Tcon_tsebra.cfg`) was made with lower support needed for intron inclusion in order to reduce the number of single exon genes.
3. The *T. conura* genome was annotated with BRAKER3 using the braker3.sif: `annnotateAssembly_BRAKER3.sh`. Make sure to recover the transcripts_merged.gtf before it is deleted. 
5. The *T. conura* merged annotation was modified with a lower intron inclusion level (`processBraker.sh`) and untranslated regions were added with intervalTree (`decorate_braker.sh`).
5. Annotation completeness was assessed with BUSCO against diptera and insecta orthodatabases (v. 10). (`checkAnnotation_BUSCO.sh`)
