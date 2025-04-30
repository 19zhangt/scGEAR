# ------------------------------------------------------------------------------
# File:        findAtail_collapse.R
# Date:        2023-10-01
# Description: This script performs data collapse
# ------------------------------------------------------------------------------

# Set global options
argvs <- commandArgs(trailingOnly = T)
set.seed(123)

# Load required libraries
suppressMessages(library(tidyverse))
suppressMessages(library(pbmcapply))
suppressMessages(library(seqinr))
suppressMessages(library(vroom))

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

## cs from studies
cs_res <- vroom::vroom(sprintf("%s/findtail_sites_in_last_exon.bed", data_dir), col_names = NA) %>% as.data.frame()
cs_res <- unique(cs_res[, 1:6]) %>% arrange(X1, as.numeric(X2))

nrow(cs_res)
clean_cs_res <- para_out(input_data = cs_res, site_dist = 100, 
                         thread = as.numeric(f_threads), tmp_dir = data_dir) # 389207 --> 67296
nrow(clean_cs_res)
clean_cs_res_info <- para_out(input_data = clean_cs_res[,1:6], site_dist = 200, 
                              thread = as.numeric(f_threads), tmp_dir = data_dir) ## 67296 --> 55069
nrow(clean_cs_res_info)

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

IPS_clean_cs_res <- subset(IPS_clean_cs_res, count > 5)
write.table(IPS_clean_cs_res, sprintf("%s/high_abundance_IPS_sites.bed", data_dir), quote = F, row.names = F, col.names = F, sep = "\t")


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
refseq_polydb_modify <- refseq_polydb_modify[rownames(refseq_polydb_modify)%in%clean_cs_res_info$name, ]


sub_clean_cs_res <- subset(clean_cs_res_info, !name%in%sub_IPS_names)
rownames(sub_clean_cs_res) <- sub_clean_cs_res$name
sub_clean_cs_res[rownames(refseq_polydb_modify), 1:4] <- refseq_polydb_modify[, 1:4]

# ggvenn::ggvenn(list("raw"=rownames(sub_clean_cs_res), "IPS"=rownames(refseq_polydb_modify)))

## PAS and CFi/II motif
test_motif <- list(c("up", 50, "A[AT]TAAA"),
                   c("up", 100, "TGTA"),
                   c("down", 50, "T[TG]T[TG]T[TG]"))
motif_summary <- lapply(test_motif, function(x){
    ext_strand <- x[1]
    ext_len <- as.numeric(x[2])
    use_motif <- x[3]
    cat(c(use_motif, ext_len, ext_strand), "\n")
    ##
    bed_file <- paste0(data_dir, "/PAS_ext_", ext_strand, "_", ext_len, ".bed")
    fasta_file <- paste0(data_dir, "/PAS_ext_", ext_strand, "_", ext_len, ".fasta")
    merge_input_sites <- apply(sub_clean_cs_res, 1, function(x){
        site_val <- as.numeric(x[3])
        site_name <- x[4]
        if(x[6]=="+"){
            if(ext_strand=="up"){
                return(c(x[1], site_val-ext_len-1, site_val, site_name, ".", x[6]))
            }else{
                return(c(x[1], site_val-1, site_val + ext_len, site_name, ".", x[6]))
            }
        }else{
            if(ext_strand=="up"){
                return(c(x[1], site_val-1, site_val + ext_len, site_name, ".", x[6]))
            }else{
                return(c(x[1], site_val-ext_len-1, site_val, site_name, ".", x[6]))
            }
        }
    }) %>% t() %>% as.data.frame()
    merge_input_sites$V2 <- gsub(" ", "", format(as.numeric(merge_input_sites$V2), scientific = F))
    merge_input_sites$V3 <- gsub(" ", "", format(as.numeric(merge_input_sites$V3), scientific = F))
    write.table(merge_input_sites, file = bed_file, quote = F, sep = "\t", row.names = F, col.names = F)
    system(command = paste0("bedtools getfasta -s -name -fi ",
                            genome_file,
                            " -bed ", bed_file, " -fo ", fasta_file))
    fasta_res <- seqinr::read.fasta(fasta_file, as.string = T)
    fasta_df <- data.frame("Name"=names(fasta_res), "seq"=toupper(unlist(fasta_res)))
    fasta_df$source <- sapply(strsplit(fasta_df$Name, "\\|"), function(x){x[1]})
    fasta_df$Name <- gsub("\\(\\+\\)", "", fasta_df$Name)
    fasta_df$Name <- gsub("\\(-)", "", fasta_df$Name)
    fasta_df$motif <- use_motif 
    fasta_df$have_motif <- sapply(fasta_df$seq, function(x){
        if(use_motif == "A[AT]TAAA"){
            grepl("AA[ACGT]AAA", x) | grepl("AATA[CGT]A", x) | grepl("AATGAA", x) | grepl("A[CGT]TAAA", x) | grepl("[CG]ATAAA", x) | grepl("T[AT]TAAA", x)
        }else{
            grepl(use_motif, x)
        }
        
    })
    fasta_df
})
motif_summary_df <- motif_summary %>% bind_rows()
motif_summary_res <- subset(motif_summary_df, motif=="A[AT]TAAA"&have_motif)
motif_summary_res$Name <- sapply(strsplit(motif_summary_res$Name, split = "::"), function(x){x[1]})

sub_clean_cs_res <- subset(sub_clean_cs_res, name%in%motif_summary_res$Name)
nrow(sub_clean_cs_res)
write.table(sub_clean_cs_res, sprintf("%s/findtail_PASs.bed", data_dir), quote = F, row.names = F, col.names = F, sep = "\t")
