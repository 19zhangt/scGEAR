
library(tidyverse)
library(pbmcapply)
library(vroom)

argvs <- commandArgs(trailingOnly=T)
tmp_files <- list.files(sprintf("%s/site_report/", argvs[1]))

read_tail_res <- function(x){
    read_support <- 5
    tmp_df <- vroom::vroom(sprintf("%s/site_report/%s", argvs[1], x)) %>% as.data.frame()
    if(!any(grepl("chr", tmp_df$chr))){
        tmp_df$chr <- paste0("chr", tmp_df$chr)
    }
    tmp_df <- subset(tmp_df, chr%in%paste0("chr", c(1:22, "X", "Y")))
    tmp_df$batch <- gsub("_polyA_positions.csv", "", x)
    tmp_clean_df <- tmp_df %>% mutate(total=sum(count)) %>% 
        mutate(CPM=count/total*10^6) %>% 
        ungroup() %>% as.data.frame()
    tmp_clean_df <- subset(tmp_clean_df, count >= read_support)
    tmp_clean_df
}

tmp_res_list <- pbmclapply(tmp_files, read_tail_res, 
                         mc.cores = min(30, length(tmp_files)), 
                         mc.preschedule = TRUE) %>% bind_rows()

tmp_aggr_res <- tmp_res_list %>% group_by(chr, strand, coord) %>% 
    summarise(count=sum(count), support=length(unique(batch)), 
              type=paste0(sort(unique(type)), collapse = ":"), CPM=sum(CPM))

## read supported at least 5 samples
total_source_support <- min(max(tmp_aggr_res$support), 5)
tmp_aggr_clean <- subset(tmp_aggr_res, support >= total_source_support)

## Output result
tmp_bed <- tmp_aggr_clean %>% 
    mutate(start=coord-1, end=coord, 
           name=paste0(chr, ":", coord, ":", strand)) %>% 
    select(chr, start, end, name, CPM, strand)
tmp_bed <- tmp_bed %>% arrange(chr, as.numeric(start))
tmp_bed$start <- gsub(" ", "", format(tmp_bed$start, scientific = F))
tmp_bed$end <- gsub(" ", "", format(tmp_bed$end, scientific = F))

colnames(tmp_bed) <- c("chr", "start", "end", "name", "CPM", "strand")
out_bed <- tmp_bed %>% arrange(chr, as.numeric(start))
out_bed$start <- gsub(" ", "", format(out_bed$start, scientific = F))
out_bed$end <- gsub(" ", "", format(out_bed$end, scientific = F))

dim(out_bed)
write.table(out_bed, sprintf("%s/findtail_sites.bed", argvs[1]), quote = F, col.names = F, row.names = F, sep = "\t")
