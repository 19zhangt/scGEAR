
## Merge ####
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
        if(!all(is.na(tmp_data$count))){
            if(tmp_strand == "+"){
                max_index <- which.max(tmp_data$count)[length(which.max(tmp_data$count))]
            }else{
                max_index <- which.max(tmp_data$count)[1]
            }
            max_coord <- tmp_data$end[max_index]
            max_name <- tmp_data$name[max_index]
            inner_index <- abs(tmp_data$end-max_coord) <= site_dist
            tmp_count <- tmp_data$count[max_index]
            tmp_vec <- c(tmp_vec, tmp_chr, max_coord-1, max_coord, max_name, tmp_count, tmp_strand, sum(inner_index))
            write.table(tmp_data$name[!inner_index], tmp_file, quote = F, col.names = F, row.names = F)
        }else{
            if(tmp_strand == "+"){
                max_index <- nrow(tmp_data)
            }else{
                max_index <- 1
            }
            max_coord <- tmp_data$end[max_index]
            max_name <- tmp_data$name[max_index]
            inner_index <- abs(tmp_data$end-max_coord) <= site_dist
            tmp_count <- tmp_data$count[max_index]
            tmp_vec <- c(tmp_vec, tmp_chr, max_coord-1, max_coord, max_name, tmp_count, tmp_strand, sum(inner_index))
            write.table(tmp_data$name[!inner_index], tmp_file, quote = F, col.names = F, row.names = F)
        }
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
    out_clean_data <- pbmcapply::pbmclapply(chr_strand_set, para_extract, 
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
library(GenomicRanges)
library(rtracklayer)


## loading data ----------------------------------------------------------------
argvs <- commandArgs(trailingOnly = T)

known_cs1 <- read.table(sprintf("%s/data/Refseq/refseq_tts_in_last_exon.bed", argvs[1]), sep = "\t")
known_cs1 <- unique(known_cs1[, 1:6]) # 56208

known_cs2 <- read.table(sprintf("%s/data/polyA_DB/UCSC_polyDB_in_last_exon.bed", argvs[1]), sep = "\t")
known_cs2 <- unique(known_cs2[, 1:6]) # 33908

known_cs <- rbind(known_cs1[, 1:6], known_cs2[, 1:6]) %>% unique()
known_cs$V4 <- paste0(known_cs$V1, ":", known_cs$V3, ":", known_cs$V6)
known_cs$V5 <- 10

seq_data_cs1 <- read.table(sprintf("%s/data/findpolyAtail/findtail_PASs.bed", argvs[1]), sep = "\t") # 29865
summary(seq_data_cs1$V5)
seq_data_cs2 <- read.table(sprintf("%s/data/mws_cs_annotation/mws_cs_info.bed", argvs[1]), sep = "\t") # 21936
summary(seq_data_cs2$V5)
seq_data_cs2$V5 <- seq_data_cs2$V5/2
seq_data_cs <- rbind(seq_data_cs1[, 1:6], seq_data_cs2[, 1:6]) %>% unique()


## merge all results
all_sites <- rbind(known_cs, seq_data_cs) %>% 
    as.data.frame() %>% arrange(V1, as.numeric(V2)) %>% unique()
colnames(all_sites) <- c("chr", "start", "end", "name", "count", "strand")
nrow(all_sites)

tmp_dist <- abs(diff(all_sites$end))
hist(tmp_dist[tmp_dist < 1000 & tmp_dist > 10], breaks = 500, xlim = c(10,300))

nrow(all_sites)
clean_all_sites <- para_out(input_data = all_sites, site_dist = 200, thread = 22) ## 139682 --> 75043
nrow(clean_all_sites)


clean_all_sites$start <- gsub(" ", "", format(clean_all_sites$start, scientific = F))
clean_all_sites$end <- gsub(" ", "", format(clean_all_sites$end, scientific = F))

Annot_dir <- sprintf("%s/data/Annotation", argvs[1])
if(!dir.exists(Annot_dir)){dir.create(Annot_dir)}

write.table(clean_all_sites, sprintf("%s/Potential_PAS.bed", Annot_dir), 
            quote = F, row.names = F, col.names = F, sep = "\t")


#############################################
## PAS annotation ####
gene_annotation <- read.table(sprintf("%s/data/Refseq/refseq.annotation.txt", argvs[1]))
rownames(gene_annotation) <- gene_annotation$V5

system(sprintf("bedtools intersect -wo -s -a %s/Potential_PAS.bed -b %s/data/Refseq/last_exon_annotation.bed > %s/Potential_PAS_last_exon.bed", 
               Annot_dir, argvs[1], Annot_dir))
    

anno_res <- read.table(sprintf("%s/Potential_PAS_last_exon.bed", Annot_dir), sep = "\t")
anno_res <- anno_res[, -ncol(anno_res)]
colnames(anno_res) <- c("site_chr", "site_start", "site_end", "site_name", "site_count", "site_strand",
                        "site_cluster", "LE_chr", "LE_start", "LE_end", "LE_name", "LE_score", 
                        "LE_strand", "LE_transcript")

table(table(anno_res$site_name))
table(anno_res$site_name)[table(anno_res$site_name)==2]
anno_res$gene <- gsub(":I[0-9]+", "", anno_res$LE_name)
anno_res$annotation <- gene_annotation[anno_res$gene, "V6"]
head(anno_res)
clean_anno_res <- anno_res %>% group_by(site_name) %>% mutate(len=length(site_name)) %>% 
    filter(len==1 | annotation!="pseudogene")
nrow(clean_anno_res)
setdiff(anno_res$site_name, clean_anno_res$site_name)
table(table(clean_anno_res$site_name))
table(table(clean_anno_res$LE_name)==1)


select_anno_res <- clean_anno_res
select_anno_res$type <- "PAS"

## tandem APA and alternative last exon
rmNPS_anno_res <- select_anno_res %>% filter(type!="NPS")
gene_filter <- rmNPS_anno_res %>% group_by(gene) %>% summarise(len=length(unique(LE_name)))
gene_filter <- gene_filter$gene[gene_filter$len > 1]
length(unique(gene_filter))
gene_filter <- unique(rmNPS_anno_res$LE_name[rmNPS_anno_res$gene%in%gene_filter])
APA <- table(rmNPS_anno_res$LE_name)[table(rmNPS_anno_res$LE_name) >= 2]
APA <- names(APA)

select_anno_res$APA_type <- "SPAS"
select_anno_res$APA_type[select_anno_res$LE_name%in%gene_filter] <- "ALE"
select_anno_res$APA_type[select_anno_res$LE_name%in%APA] <- "TUTR"
table(select_anno_res$APA_type)
# ALE  SPAS  TUTR 
# 17629 24069 34531


PAS_num <- table(select_anno_res$LE_name[select_anno_res$APA_type=="TUTR"])
PAS_num[PAS_num > 10] <- 10
pie(table(PAS_num))
sum(PAS_num==2)/length(PAS_num)

select_anno_res <- select_anno_res %>%
    group_by(LE_name, type) %>%
    mutate(
        sort_key = if (site_strand[1] == "+") site_end else -site_end
    ) %>%
    arrange(sort_key, .by_group = TRUE) %>%
    mutate(order = row_number()) %>%
    ungroup() %>%
    select(-sort_key) %>% as.data.frame()

select_anno_res$rename <- paste0(select_anno_res$LE_name, ":", 
                                 select_anno_res$type, select_anno_res$order)

publication_data <- select_anno_res %>% as.data.frame()
publication_data <- publication_data %>% group_by(LE_name) %>% 
    mutate(mean_dist=mean(site_end[type=="PAS"])) %>% 
    mutate(PAS_type=case_when(site_strand[1] == "+" & site_end < mean_dist ~ "proximal",  
                              site_strand[1] == "+" & site_end > mean_dist ~ "distal",
                              site_strand[1] == "-" & site_end < mean_dist ~ "distal",  
                              site_strand[1] == "-" & site_end > mean_dist ~ "proximal"))
publication_data$PAS_type[publication_data$APA_type!="TUTR"] = "-"
publication_data <- publication_data %>% select(chr=site_chr, site_position=site_end, site_name=rename, 
                                                strand=site_strand, site_type=type, loc_type = PAS_type,
                                                LE_start, LE_end, LE_name, LE_type=APA_type, 
                                                gene_name=gene, gene_type=annotation)
publication_data <- as.data.frame(publication_data)

write.table(publication_data, sprintf("%s/KeepPAS_data_annotation.txt", Annot_dir), quote = F, sep = "\t", row.names = F)

table(publication_data[c("site_type","loc_type")])
length(unique(publication_data$LE_name[publication_data$LE_type=="TUTR"]))
length(unique(publication_data$LE_name[publication_data$LE_type=="TUTR"]))
length(unique(publication_data$gene_name[publication_data$LE_type=="TUTR"]))
length(unique(publication_data$gene_name[publication_data$LE_type=="ALE"]))


system(sprintf("bedtools closest -s -D a -a %s -b %s > %s",
               sprintf("%s/data/findpolyAtail/high_abundance_IPS_sites.bed", argvs[1]),
               sprintf("%s/Potential_PAS_last_exon.bed", Annot_dir),
               sprintf("%s/IPS_site_intersect.bed", Annot_dir)))

high_abundance_IPS <- read.table( sprintf("%s/IPS_site_intersect.bed", Annot_dir), header = F, sep = "\t")
high_abundance_around <- subset(high_abundance_IPS, abs(V23) < 200)
high_abundance_IPS <- high_abundance_IPS[!high_abundance_IPS$V11%in%high_abundance_around$V11, ]
high_abundance_IPS <- high_abundance_IPS[, c(1:4,18,6,16,17)]
colnames(high_abundance_IPS) <- c("site_chr", "site_start", "site_end", "site_name", "LE_name", "site_strand", "LE_start", "LE_end")
high_abundance_IPS$site_name <- paste0(high_abundance_IPS$site_name, ":NPS")

select_anno_res1 <- rbind(select_anno_res %>% select(site_chr, site_start, site_end, site_name=rename, LE_name, site_strand, LE_start, LE_end),
                         high_abundance_IPS)
select_anno_res1 <- select_anno_res1 %>% arrange(site_chr, site_start)
select_anno_res1 <- select_anno_res1 %>% group_by(LE_name) %>% 
    mutate(diff=case_when(length(site_strand)==1 & site_strand[1] == "+" ~ site_end-LE_start-1,
                          length(site_strand)==1 & site_strand[1] == "-" ~ LE_end-site_end,
                          site_strand[1] == "+" ~ c(site_end[1]-LE_start[1]-1, 
                                                    diff(site_end)), 
                          site_strand[1] == "-" ~ c(diff(site_end), 
                                                    LE_end[length(LE_end)]-site_end[length(LE_end)])),
           index=case_when(length(site_strand)==1 ~ 1,
                           site_strand[1] == "+" ~ c(1, rep(0, length(LE_end)-1)), 
                           site_strand[1] == "-" ~ c(rep(0, length(LE_end)-1), 1)))

select_anno_tmp1 <- subset(select_anno_res1, diff < 400 & index==1)
select_anno_tmp1 <- subset(select_anno_tmp1, !grepl("NPS", site_name))

# Convert data to GRanges objects
pas_granges <- GRanges(
    seqnames = select_anno_tmp1$site_chr,
    ranges = IRanges(
        start = select_anno_tmp1$site_start + 1L, 
        end = select_anno_tmp1$site_end
    ),
    strand = select_anno_tmp1$site_strand,
    site_name = select_anno_tmp1$site_name,
    LE_name = select_anno_tmp1$LE_name,
    LE_start = select_anno_tmp1$LE_start,
    LE_end = select_anno_tmp1$LE_end
)

# load GTF
gtf <- import(sprintf("%s/data/Refseq/GRCh38_rechr.gtf", argvs[1]))
exons <- gtf[gtf$type == "exon"]
exons_by_transcript <- split(exons, exons$transcript_id)


isoform_transcript <- read.table(sprintf("%s/data/Refseq/last_exon_annotation.bed", argvs[1]), sep = "\t")
isoform_transcript <- isoform_transcript %>% separate_rows(V7, sep = ",") %>% select(UTR=V4, transcript=V7) %>% unique() %>% as.data.frame()
rownames(isoform_transcript) <- isoform_transcript$transcript
isoform_transcript <- split(isoform_transcript, isoform_transcript$UTR)

exons_by_UTR <- pbmcapply::pbmclapply(names(isoform_transcript), function(x){
    if(nrow(isoform_transcript[[x]]) == 1){
        return(exons_by_transcript[[isoform_transcript[[x]][1,2]]][, 1])
    }else{
        return(reduce(unlist(exons_by_transcript[isoform_transcript[[x]][,2]])[, 1]))
    }
}, mc.cores = 10, mc.preschedule = TRUE)
names(exons_by_UTR) <- names(isoform_transcript)


# Optimized processing function
PAS_structure_info_optimized <- function(tmp_index) {
    # Directly obtain preprocessed GRanges data
    pas_site_info <- pas_granges[valid_indices[tmp_index]]
    UTR_id <- pas_site_info$LE_name
    gene_exons <- exons_by_UTR[[UTR_id]]
    
    # Sorting exons according to strand
    strand_dir <- as.character(strand(pas_site_info))
    transcript_exons <- if (strand_dir == "+") {
        sort(gene_exons, decreasing = TRUE)
    } else {
        sort(gene_exons)
    }
    
    # Handle overlaps
    overlaps <- findOverlaps(pas_site_info, transcript_exons)
    if (length(overlaps) == 0) return(NULL)
    exon_idx <- subjectHits(overlaps)[1]
    
    # Modify exon boundaries
    if (strand_dir == "+") {
        end(transcript_exons[exon_idx]) <- end(pas_site_info)
    } else {
        start(transcript_exons[exon_idx]) <- start(pas_site_info)
    }
    
    # Calculate cumulative width
    widths <- width(transcript_exons)
    cum_widths <- cumsum(widths)
    total_width <- cum_widths[length(cum_widths)]
    
    if (total_width >= 400) {
        cut_index <- which.max(cum_widths >= 400)
        transcript_exons <- transcript_exons[1:cut_index]
        remaining <- 400 - if (cut_index == 1) 0L else cum_widths[cut_index - 1L]
        
        if (strand_dir == "+") {
            start(transcript_exons[cut_index]) <- end(transcript_exons[cut_index]) - remaining + 1L
        } else {
            end(transcript_exons[cut_index]) <- start(transcript_exons[cut_index]) + remaining - 1L
        }
    }
    
    # Constructing the result
    data.frame(
        site_chr = as.character(seqnames(pas_site_info)),
        site_start = start(transcript_exons),
        site_end = end(transcript_exons),
        site_name = pas_site_info$site_name,
        LE_name = pas_site_info$LE_name,
        site_strand = strand_dir,
        LE_start = pas_site_info$LE_start,
        LE_end = pas_site_info$LE_end,
        stringsAsFactors = FALSE
    )
}

# Preprocessing step: extracting gene_id and filtering valid indexes
valid_indices <- which(select_anno_tmp1$LE_name %in% names(isoform_transcript))

# Parallel execution
select_anno_tmp1 <- pbmcapply::pbmclapply(
    X = seq_along(valid_indices),
    FUN = PAS_structure_info_optimized,
    mc.cores = 10,
    mc.preschedule = TRUE  
) %>% data.table::rbindlist() %>% dplyr::distinct()
select_anno_tmp1 <- as.data.frame(select_anno_tmp1)

## PAS within UTR
select_anno_tmp2 <- subset(select_anno_res1, !(diff < 400 & index==1))
select_anno_tmp2 <- subset(select_anno_tmp2, !grepl("NPS", site_name))

select_anno_tmp2 <- apply(select_anno_tmp2, 1, function(x){
    names(x) <- colnames(select_anno_tmp2)
    chr <- x["site_chr"]
    ##
    if(as.numeric(x["diff"]) >= 400){
        region_len <- 400
    }else{
        region_len <- as.numeric(x["diff"])
    }
    site <- as.numeric(x["site_end"])
    if(x["site_strand"]=="+"){
        start_v <- site - region_len + 1
        end_v <- site
    }else{
        start_v <- site
        end_v <- site + region_len - 1
    }
    return(c(chr, start_v-1, end_v, x["site_name"], x["LE_name"], x["site_strand"], x["LE_start"], x["LE_end"]))
}) %>% t() %>% as.data.frame()
colnames(select_anno_tmp2) <- colnames(select_anno_tmp1)


select_anno_tmp <- rbind(select_anno_tmp1, select_anno_tmp2)

len_num <- as.numeric(select_anno_tmp2$site_end) - as.numeric(select_anno_tmp2$site_start)
table(len_num < 200)

## GTF of PASs
select_anno_gr <- GRanges(seqnames = select_anno_tmp$site_chr,
                          ranges = IRanges(start = as.numeric(select_anno_tmp$site_start)+1,
                                           end = as.numeric(select_anno_tmp$site_end),
                                           names = select_anno_tmp$site_name),
                          strand = select_anno_tmp$site_strand)
rtracklayer::export(select_anno_gr, con = sprintf("%s/PAS_annotation.gtf", Annot_dir), format = "gtf")
system(command = sprintf("sed -i s\"/sequence_feature/exon/;s/ID /transcript_id /\" %s/PAS_annotation.gtf", Annot_dir))
