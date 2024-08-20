library(tximport)
library(tximportData)
library(GenomicFeatures)
library(readr)
library(magrittr)
library(DESeq2)


mydir <- [redacted]
annot_dir <- paste0(mydir,"annotation_Tconura/braker/Tcon_braker3_1/")
annot <- "braker.adj.UTR.mod.gtf"

salmon_dir <- paste0(mydir, "transcriptomics_Tconura/01_RNA/07_Salmon/")
salmon_suffix <- "_salmon_quant_ISF"
samples_file <- paste0(mydir, "/transcriptomics_Tconura/00_metaData/Tconura_RNAsampleMetadata.txt")

#### Load TxDb from annotation and subset tx2gene ####
Tcon_txdb <- GenomicFeatures::makeTxDbFromGFF(
  file = paste0(annot_dir, annot),
  format=c("gtf"),
  dataSource= "Tcon_braker")

k <- keys(Tcon_txdb, keytype = "TXNAME")
tx2gene <- select(Tcon_txdb, k, "GENEID", "TXNAME")

#### Load samples data frame ####
samples_df <- read_tsv(samples_file, col_names = T) %>%
  dplyr::filter(Stage == "Larva") %>%
  dplyr::select(Run, Sample, Name = `Submitted Name`, Host_form, CF, Treatment, Family, Notes)

#### Import transcript counts ####
files <- paste0(salmon_dir, samples_df$Run, salmon_suffix, "/quant.sf")
names(files) <- samples_df$Sample

txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene)

#### Construct DESeqDataSet ####
ddsTxi <- DESeqDataSetFromTximport(
  txi,
  colData = samples_df,
  design = ~ Treatment)

save(ddsTxi,file= paste0(salmon_dir, "counts/salmon_braker_txi_dds_cfLarvae.Rdata"))

# mapping rates
cat fileList.txt | while read line;
  do name=$(echo $line | cut -d "_" -f 1,2); 
    tot_reads=$(grep "total reads" ${line}/logs/salmon_quant.log |cut -d " " -f 6  | tr -d ",");
    map_rate=$(grep "Mapping rate" ${line}/logs/salmon_quant.log |cut -d " " -f 8  | tr -d "%");
  echo $name $tot_reads $map_rate >> counts/mapping_rate.txt
  done
