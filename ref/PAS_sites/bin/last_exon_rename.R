
library(tidyverse)
argvs <- commandArgs(trailingOnly=T)

merge_bed <- read.table(argvs[1])
gene_transcript <- read.table(argvs[2])

## multiple genes
merge_bed_rename <- apply(merge_bed, 1, function(x){
    if(grepl(",", x[4])){
        trans_list <- strsplit(x[7], split = ",")[[1]]
        sub_gene_trans <- subset(gene_transcript, V7%in%trans_list)
        tmp_list <- strsplit(x[4], split = ",")[[1]]
        if(all(grepl("-", tmp_list))){
            x[4] <- tmp_list[1]
        }else if(any(grepl("-", tmp_list))){
            x[4] <- tmp_list[!grepl("-", tmp_list)][1]
        }else{
            if(all(grepl("LOC", tmp_list))){
                x[4] <- tmp_list[1]
            }else if(any(grepl("LOC", tmp_list))){
                x[4] <- tmp_list[!grepl("LOC", tmp_list)][1]
            }else{
                x[4] <- tmp_list[1]
            }
        }
        x[7] <- paste0(sub_gene_trans$V7[sub_gene_trans$V4==x[4]], collapse = ",")
    }
    return(x)
}) %>% t() %>% as.data.frame()
merge_bed_rename$V2 <- as.numeric(merge_bed_rename$V2)
merge_bed_rename$V3 <- as.numeric(merge_bed_rename$V3)

## order
merge_bed_rename <- merge_bed_rename %>% arrange(V1, V2)
merge_bed_rename <- merge_bed_rename %>% group_by(V4) %>% 
            mutate(
                order=case_when(
                    "+"%in%V6 ~ order(V2, decreasing = T), 
                    "-"%in%V6 ~ order(V2))
                )
merge_bed_rename$UTRname <- paste0(merge_bed_rename$V4, ":I", merge_bed_rename$order)
write.table(merge_bed_rename[, c("V1", "V2", "V3", "UTRname", "V5", "V6", "V7")], argvs[3], quote = F, sep = "\t", row.names = F, col.names = F)

