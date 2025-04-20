#!/usr/bin/env Rscript

##
para_extract <- function(chr_strand, input_data, site_dist){
    tmp_chr <- strsplit(chr_strand, split = ":")[[1]][1]
    tmp_strand <- strsplit(chr_strand, split = ":")[[1]][2]
    sub_data <- subset(input_data, chr==tmp_chr&strand==tmp_strand)
    tmp_file <- sprintf("tmp_%s.txt", chr_strand)
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

para_out <- function(input_data, site_dist, thread){
    colnames(input_data) <- c("chr", "start", "end", "name", "count", "strand")
    for(i in c(2,3,5)){input_data[,i] <- as.numeric(input_data[,i])}
    input_data <- input_data %>% arrange(chr, start)
    chr_strand_set <- paste0(input_data$chr, ":", input_data$strand) %>% unique()
    out_clean_data <- parallel::mclapply(chr_strand_set, para_extract, 
                                         input_data = input_data,
                                         site_dist = site_dist,
                                         mc.cores = thread, 
                                         mc.preschedule = F) %>% bind_rows()
    out_clean_data <- as.data.frame(out_clean_data)
    colnames(out_clean_data) <- c("chr", "start", "end", "name", "count", "strand", "support")
    for(i in c(2,3,5,7)){out_clean_data[,i] <- as.numeric(out_clean_data[,i])}
    out_clean_data <- out_clean_data %>% arrange(chr, start)
    out_clean_data
}


library(tidyverse)
argvs <- commandArgs(trailingOnly = T)

## cs from five studies
cs_res <- vroom::vroom(sprintf("%s/findtail_sites_in_last_exon.bed", argvs[1]), col_names = NA) %>% as.data.frame()
cs_res <- unique(cs_res[, 1:6]) %>% arrange(X1, as.numeric(X2))
nrow(cs_res)

clean_cs_res_info <- para_out(input_data = cs_res, site_dist = 200, thread = 20) # 357463 --> 48069
nrow(clean_cs_res_info)

# dist_summary <- abs(diff(as.numeric(clean_cs_res_info$start)))
# hist(dist_summary[dist_summary<1000], breaks = 1000, xlim = c(0,500))

clean_cs_res <- subset(clean_cs_res_info, count > 0.5)
clean_cs_res$start <- gsub(" ", "", format(clean_cs_res$start, scientific = F))
clean_cs_res$end <- gsub(" ", "", format(clean_cs_res$end, scientific = F))

## find internal priming
ext_len <- 20
input_sites <- apply(clean_cs_res, 1, function(x){
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

write.table(input_sites, file = sprintf("%s/PAS_ext_%s.bed", argvs[1], ext_len), 
            quote = F, sep = "\t", row.names = F, col.names = F)

system(command = sprintf("sort -k1,1 -k2,2n %s/PAS_ext_%s.bed |bedtools closest -s -D a -a - -b %s/../Refseq/refseq_tts.bed > %s/PAS_ext_%s_refseq.closest",
                         argvs[1], ext_len, argvs[1], argvs[1], ext_len))
IPS_NF1 <- read.table(sprintf("%s/PAS_ext_%s_refseq.closest", argvs[1], ext_len))
IPS_NF1 <- IPS_NF1[-10<=IPS_NF1$V13 & IPS_NF1$V13<=0, ]

system(command = sprintf("sort -k1,1 -k2,2n %s/PAS_ext_%s.bed |bedtools closest -s -D a -a - -b %s/../../src/UCSC_polyDB_hg38.bed > %s/PAS_ext_%s_polyadb.closest",
                         argvs[1], ext_len, argvs[1], argvs[1], ext_len))
IPS_NF2 <- read.table(sprintf("%s/PAS_ext_%s_polyadb.closest", argvs[1], ext_len))
IPS_NF2$V13 <- IPS_NF2$V15
IPS_NF2 <- IPS_NF2[-10<=IPS_NF2$V13 & IPS_NF2$V13<=0, ]


system(command = sprintf("bedtools getfasta -s -name -fi %s -bed %s/PAS_ext_%s.bed -fo %s/PAS_ext_%s.fasta", 
                         argvs[2], argvs[1], ext_len, argvs[1], ext_len))
system(command = sprintf("python %s/../../bin/internal_priming.py %s/PAS_ext_%s.fasta %s/PAS_ext_%s_internal_priming.out", 
                         argvs[1], argvs[1], ext_len, argvs[1], ext_len))
       
IPS_names <- read.table(paste0(argvs[1], "/PAS_ext_", ext_len, "_internal_priming.out"), sep = "\t") 
IPS_names <- gsub("::.+", "", unique(IPS_names$V1))
IPS_names <- gsub("\\(\\+\\)", "", IPS_names)
IPS_names <- gsub("\\(-)", "", IPS_names)
sub_IPS_names <- setdiff(IPS_names, c(IPS_NF1$V4, IPS_NF2$V4))

# ggvenn::ggvenn(list("raw"=unique(clean_cs_res$name), "IPS"=sub_IPS_names))

IPS_bed <- lapply(sub_IPS_names, function(x){
    ips_info <- strsplit(x, split = ":")[[1]]
    ips_loc <- as.numeric(ips_info[2])
    return(c(ips_info[1], ips_loc - 1, ips_loc, x, ".", ips_info[3]))
})
IPS_bed <- do.call('rbind', IPS_bed) %>% as.data.frame()

IPS_clean_cs_res <- subset(clean_cs_res, name%in%sub_IPS_names)
nrow(IPS_clean_cs_res)

IPS_clean_cs_res <- subset(IPS_clean_cs_res, count > 10)
write.table(IPS_clean_cs_res, sprintf("%s/high_abundance_IPS_sites.bed", argvs[1]), quote = F, row.names = F, col.names = F, sep = "\t")


## modify location of cs
refseq_polydb_location <- rbind(IPS_NF1, IPS_NF2[!IPS_NF2$V4%in%IPS_NF1$V4, 1:13]) %>% unique()
refseq_polydb_location$diff <- apply(refseq_polydb_location, 1, function(x){
    site_loc <- strsplit(x[4], split = ":")[[1]][2] %>% as.numeric()
    if(site_loc>=as.numeric(x[8])&site_loc<=as.numeric(x[9])){
        return(0)
    }else{
        min(abs(site_loc-as.numeric(x[8])), abs(site_loc-as.numeric(x[9])))
    }
})
refseq_polydb_location <- subset(refseq_polydb_location, !V4%in%refseq_polydb_location$V4[refseq_polydb_location$diff==0])
refseq_polydb_location <- refseq_polydb_location[!duplicated(refseq_polydb_location$V4), ] %>% as.data.frame()
refseq_polydb_modify <- apply(refseq_polydb_location , 1, function(x){
    site_start <- as.numeric(x[8])
    site_end <- as.numeric(x[9])
    if(x[6]=="+"){
        return(c(x[1], site_end-1, site_end, 
                 sprintf("%s:%s:%s", x[1], site_end, x[6]), ".", x[6]))
    }else{
        return(c(x[1], site_start, site_start+1, 
                 sprintf("%s:%s:%s", x[1], site_start+1, x[6]), ".", x[6]))
    }
}) %>% t() %>% as.data.frame()
rownames(refseq_polydb_modify) <- refseq_polydb_location$V4

sub_clean_cs_res <- subset(clean_cs_res, !name%in%sub_IPS_names)
rownames(sub_clean_cs_res) <- sub_clean_cs_res$name
sub_clean_cs_res[rownames(refseq_polydb_modify), 1:4] <- refseq_polydb_modify[, 1:4]

nrow(sub_clean_cs_res)
write.table(sub_clean_cs_res, sprintf("%s/findtail_PASs.bed", argvs[1]), quote = F, row.names = F, col.names = F, sep = "\t")

