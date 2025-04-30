
library(tidyverse)
library(bedtoolsr)
library(pbmcapply)

argvs <- commandArgs(trailingOnly=T)

merge_bed <- read.table(argvs[1])
colnames(merge_bed) <- c("chr", "start", "end", "Name", "score", "strand", "transcript")

merge_bed_list <- split(merge_bed, merge_bed$Name)
merge_bed_res <- pbmclapply(merge_bed_list, function(x){
    bt.merge(x, c="4,5,6,7", o="distinct,distinct,distinct,distinct")
}, mc.cores = 10, mc.preschedule = F) %>% bind_rows()

merge_bed_intersect <- bt.intersect(a=merge_bed ,b=merge_bed, wo = T, s = T)
merge_bed_intersect <- subset(merge_bed_intersect, V4!=V11)

merge_bed1 <- subset(merge_bed, Name%in%merge_bed_intersect$V4)
merge_bed1 <- split(merge_bed1, merge_bed1$Name)
merge_bed1 <- lapply(merge_bed1, function(x){
    bt.merge(x, c="4,5,6,7", o="distinct,distinct,distinct,distinct")
}) %>% bind_rows()

merge_bed2 <- subset(merge_bed, !Name%in%merge_bed_intersect$V4)
merge_bed2 <- bt.merge(merge_bed2, s = T, c="4,5,6,7", o="distinct,distinct,distinct,distinct")


## order
merge_bed_rename <- rbind(merge_bed1, merge_bed2) %>% arrange(V1, V2)
merge_bed_rename <- merge_bed_rename %>% group_by(V4) %>% 
            mutate(
                order=case_when(
                    "+"%in%V6 ~ order(V2, decreasing = T), 
                    "-"%in%V6 ~ order(V2))
                )
merge_bed_rename$UTRname <- paste0(merge_bed_rename$V4, ":I", merge_bed_rename$order)
write.table(merge_bed_rename[, c("V1", "V2", "V3", "UTRname", "V5", "V6", "V7")], argvs[2], quote = F, sep = "\t", row.names = F, col.names = F)
