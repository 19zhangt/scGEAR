# ------------------------------------------------------------------------------
# File:        mws_cs_site_collapse.R
# Date:        2023-10-01
# Description: This script performs data collapse
# ------------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(pbmcapply)
library(seqinr)


# Set global options
argvs <- commandArgs(trailingOnly = T)
set.seed(123)


# Define file paths
data_dir <- argvs[1]
genome_file <- argvs[2]
f_threads <- argvs[3]


# DEFINE CUSTOM FUNCTIONS
para_extract <- function(chr_strand, input_data, site_dist, tmp_dir){
    tmp_chr <- strsplit(chr_strand, split = ":")[[1]][1]
    tmp_strand <- strsplit(chr_strand, split = ":")[[1]][2]
    sub_data <- subset(input_data, chr==tmp_chr&strand==tmp_strand)
    tmp_file <- sprintf("%s/tmp_%s.txt", tmp_dir, chr_strand)
    write.table(sub_data$name, tmp_file, quote = F, col.names = F, row.names = F)
    tmp_vec <- c()
    while(file.size(tmp_file) > 0){
        retain_names <- read.table(tmp_file)
        cat(nrow(retain_names), "\n")
        tmp_data <- sub_data %>% filter(name%in%retain_names$V1)
        if(tmp_strand=="+"){
            max_index <- which.max(tmp_data$count)[length(which.max(tmp_data$count))]
        }else{
            max_index <- which.max(tmp_data$count)[1]
        }
        max_coord <- tmp_data$end[max_index]
        max_name <- tmp_data$name[max_index]
        inner_index <- abs(tmp_data$end-max_coord) <= site_dist
        tmp_count <- sum(tmp_data$count[inner_index])
        tmp_vec <- c(tmp_vec, tmp_chr, max_coord-1, max_coord, max_name, tmp_count, tmp_strand, sum(inner_index))
        write.table(tmp_data$name[!inner_index], tmp_file, quote = F, col.names = F, row.names = F)
    }
    tmp_vec <- matrix(tmp_vec, byrow = T, ncol = 7) %>% as.data.frame()
    file.remove(tmp_file)
    tmp_vec
}

para_out <- function(input_data, site_dist, thread, tmp_dir){
    colnames(input_data) <- c("chr", "start", "end", "name", "count", "strand")
    for(i in c(2,3,5)){input_data[,i] <- as.numeric(input_data[,i])}
    input_data <- input_data %>% arrange(chr, start)
    chr_strand_set <- paste0(input_data$chr, ":", input_data$strand) %>% unique()
    out_clean_data <- pbmcapply::pbmclapply(chr_strand_set, para_extract, 
                                            input_data = input_data,
                                            site_dist = site_dist,
                                            tmp_dir = tmp_dir,
                                            mc.cores = thread, 
                                            mc.preschedule = F) %>% bind_rows()
    out_clean_data <- as.data.frame(out_clean_data)
    colnames(out_clean_data) <- c("chr", "start", "end", "name", "count", "strand", "support")
    for(i in c(2,3,5,7)){out_clean_data[,i] <- as.numeric(out_clean_data[,i])}
    out_clean_data <- out_clean_data %>% arrange(chr, start)
    out_clean_data
}


## cs from mws-seq
mws_site <- vroom::vroom(sprintf("%s/mws_cs_last_exon.bed", data_dir), delim = "\t", col_names = NA) %>% as.data.frame()
mws_site <- mws_site[, 1:6] %>% unique()
mws_site <- mws_site %>% group_by(X4) %>% mutate(total=sum(X5)) %>% mutate(CPM=X5/total*10^6)

selected_cells <- c("94-dendritic-cell", "93-myeloid-cell", "13-monocyte", "22-dendritic-cell", "24-t-cell", "2-macrophage", "31-monocyte", "37-b-cell",
                    "45-macrophage", "52-proliferating-t-cell", "47-m2-macrophage", "51-macrophage", "69-macrophage", "6-t-cell", "78-macrophage", "100-b-cell", "14-b-cell-plasmocyte",
                    "3-b-cell-plasmocyte")

input_data <- subset(mws_site, X4%in%selected_cells)
# hist(input_data$CPM, xlim=c(0,10), breaks = 250000)

clean_data <- input_data %>% group_by(X1,X3,X6) %>% summarise(CPM=max(CPM)) %>% as.data.frame()
clean_data <- clean_data %>% mutate(start=X3-1, name=paste0(X1, ":", X3, ":", X6)) %>% 
    select(chr=X1, start, end=X3, name, CPM, strand=X6)

head(clean_data)
# dist_summary <- abs(diff(clean_data$end))
# hist(dist_summary[dist_summary<1000], breaks = 1000, xlim = c(0,80))

clean_mws_res_info <- para_out(input_data = clean_data, site_dist = 200, 
                               thread = as.numeric(f_threads), tmp_dir = data_dir)  # 48932 --> 26341
nrow(clean_mws_res_info)

range(clean_mws_res_info$count)
clean_mws_res <- subset(clean_mws_res_info, count > 5)

## internal priming
ext_len <- 20
input_sites <- apply(clean_mws_res, 1, function(x){
    if(x[6]=="+"){
        site_val <- as.numeric(x[3])
        site_name <- x[4] #sprintf("%s:%s:%s", x[1], site_val, x[4])
        return(c(x[1], site_val-11, site_val + ext_len, site_name, ".", x[6]))
    }else{
        site_val <- as.numeric(x[3])
        site_name <- x[4] #sprintf("%s:%s:%s", x[1], site_val, x[4])
        return(c(x[1], site_val - ext_len - 1, site_val+10, site_name, ".", x[6]))
    }
}) %>% t() %>% as.data.frame()

write.table(input_sites, file = sprintf("%s/PAS_ext_%s.bed", data_dir, ext_len), 
            quote = F, sep = "\t", row.names = F, col.names = F)

system(command = sprintf("sort -k1,1 -k2,2n %s/PAS_ext_%s.bed | bedtools closest -s -D a -a - -b %s/../refSeq/refseq_tts.bed > %s/PAS_ext_%s_refseq.closest",
                         data_dir, ext_len, data_dir, data_dir, ext_len))
IPS_NF1 <- read.table(sprintf("%s/PAS_ext_%s_refseq.closest", data_dir, ext_len))
IPS_NF1 <- IPS_NF1[-10<=IPS_NF1$V13 & IPS_NF1$V13<=0, ]

system(command = sprintf("sort -k1,1 -k2,2n %s/PAS_ext_%s.bed | bedtools closest -s -D a -a - -b reference/UCSC_polyDB_hg38.bed > %s/PAS_ext_%s_polyadb.closest",
                         data_dir, ext_len, data_dir, ext_len))
IPS_NF2 <- read.table(sprintf("%s/PAS_ext_%s_polyadb.closest", data_dir, ext_len))
IPS_NF2$V13 <- IPS_NF2$V15
IPS_NF2 <- IPS_NF2[-10<=IPS_NF2$V13 & IPS_NF2$V13<=0, ]

system(command = sprintf("bedtools getfasta -s -name -fi %s -bed %s/PAS_ext_%s.bed -fo %s/PAS_ext_%s.fasta", 
                         genome_file, data_dir, ext_len, data_dir, ext_len))
system(command = sprintf("python src/internal_priming.py %s/PAS_ext_%s.fasta %s/PAS_ext_%s_internal_priming.out", 
                         data_dir, ext_len, data_dir, ext_len))

IPS_names <- read.table(paste0(data_dir, "/PAS_ext_", ext_len, "_internal_priming.out"), sep = "\t") 
IPS_names <- gsub("::.+", "", unique(IPS_names$V1))
IPS_names <- gsub("\\(\\+\\)", "", IPS_names)
IPS_names <- gsub("\\(-)", "", IPS_names)
sub_IPS_names <- setdiff(IPS_names, c(IPS_NF1$V4[-10<=IPS_NF1$V13 & IPS_NF1$V13<=0], IPS_NF2$V4[-10<=IPS_NF2$V15 & IPS_NF2$V15<=0]))

sub_clean_mws_res <- subset(clean_mws_res, !name%in%sub_IPS_names)
sub_clean_mws_res$end <- case_when(
    sub_clean_mws_res$strand=="+" ~ sub_clean_mws_res$end+1,
    sub_clean_mws_res$strand=="-" ~ sub_clean_mws_res$end-1
)
sub_clean_mws_res$start <- sub_clean_mws_res$end-1
rownames(sub_clean_mws_res) <- sub_clean_mws_res$name

write.table(sub_clean_mws_res, sprintf("%s/mws_cs_info.bed", data_dir), quote = F, row.names = F, col.names = F, sep = "\t") # 21943
