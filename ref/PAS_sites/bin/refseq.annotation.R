
packages <- c("rtracklayer", "tidyverse")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == "rtracklayer") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("rtracklayer")
    } else {
      install.packages(pkg)
    }
  }
}


library(rtracklayer)
library(tidyverse)

argvs <- commandArgs(trailingOnly=T)

gtf.file1 <- argvs[1]
gtf.gr1 = rtracklayer::import.gff(con = gtf.file1, format = 'gtf')

## Autosomal genes
sub_gtf <- subset(gtf.gr1, type=="gene")
sub_gtf <- subset(sub_gtf, gene_biotype%in%c("protein_coding", "lncRNA")) # , "pseudogene", "transcribed_pseudogene"
table(sub_gtf$gene_biotype)
table(GenomicRanges::seqnames(sub_gtf))
sub_gtf <- as.data.frame(sub_gtf)
sub_gtf <- subset(sub_gtf, seqnames%in%paste0("chr", c(1:22, "X", "Y")))
table(sub_gtf$seqnames)

output_info <- sub_gtf %>% select(seqnames, start, end, strand, gene_id, gene_biotype)

write.table(output_info, argvs[2], 
            quote = F, sep = "\t", row.names = F, col.names = F)

