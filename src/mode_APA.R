
library(tidyverse)
library(pbmcapply)


run_apa_analysis <- function(bam_data_=bam_data, args_=args, log_=log_file) {
    #' Run complete RNA editing analysis pipeline
    #'
    #' @param bam_data_ BAM files and Cell type annotation data
    #' @param args_ parameters
    #' @param log_ log file
    
    bam_info_out <- bam_data_$bam_anno %>% select(pool,bam_dedup_path)
    write.table(bam_info_out, paste0(args_$outdir, "/bam_list.csv"), quote = F, row.names = F, sep = ",")
    
    system(command = sprintf("bash src/identify_polyA_sites.sh %s %s %s %s",
                             paste0(args_$outdir, "/bam_list.csv"),
                             args_$threads,
                             args_$genome_fa,
                             args_$outdir))
    
    
    bam_run_dir <- paste0(args_$outdir, "/BAMcount")
    dir.create(bam_run_dir)
    apply(bam_data_$bam_anno[, c("pool", "bam_dedup_path")], 1, function(x){
        dir.create(paste0(bam_run_dir, "/", x[1]))
        system(command = sprintf("featureCounts %s -a %s/Annotation/PAS_annotation.gtf -o %s/%s/assigned %s",
                                 sprintf("-M -O --largestOverlap -F GTF -t exon -g transcript_id -s 1 -T %s -R BAM", args_$threads),
                                 args_$outdir, bam_run_dir, x[1], x[2]))
        system(command = sprintf("samtools sort -@ %s %s/%s/%s.featureCounts.bam -o %s/%s/reassigned.bam",
                                 args_$threads, bam_run_dir, x[1], basename(x[2]), bam_run_dir, x[1]))
        system(command = sprintf("samtools index -@ %s %s/%s/reassigned.bam",
                                 args_$threads, bam_run_dir, x[1]))
        system(command = sprintf("rm %s/%s/%s.featureCounts.bam",
                                 bam_run_dir, x[1], basename(x[2])))
    })
    
    apply(bam_data_$bam_anno[, c("pool", "bam_dedup_path")], 1, function(x){
        system(command = sprintf("umi_tools count %s %s -I %s/%s/reassigned.bam -S %s/%s/UMItools.tsv.gz -E %s/%s/UMItools.error.log -L %s/%s/UMItools.logging.log",
                                 "--extract-umi-method=tag --umi-tag=UB --cell-tag=CB --per-gene --gene-tag=XT --assigned-status-tag=XS",
                                 "--per-cell --wide-format-cell-counts",
                                 bam_run_dir, x[1], bam_run_dir, x[1], bam_run_dir, x[1],
                                 bam_run_dir, x[1]))
    })
    
    sample_list <- paste0(bam_run_dir, "/", bam_data_$bam_anno$pool, "/UMItools.tsv.gz")
    lapply(sample_list, function(x){
        rawcount <- data.table::fread(x, header=T) %>% as.data.frame()
        rawcount <- Matrix::Matrix(rawcount %>% column_to_rownames('gene') %>% as.matrix(), sparse = TRUE)
        saveRDS(rawcount, gsub(".gz", ".rds", x))
    })
    
    
    PAS_anno <- read.table(sprintf("%s/Annotation/KeepPAS_data_annotation.txt", args_$outdir), header = T)
    meta_data <- bam_data$bam_cell_anno
    ct_out_dir <- sprintf("%s/PAS_selection", args_$outdir)
    dir.create(ct_out_dir)
    
    ind_id_list <- unique(meta_data$donor)
    cellquant_list <- pbmclapply(ind_id_list, function(ind_id){
        ind_info <- subset(meta_data, donor==ind_id) %>% as.data.frame()
        rownames(ind_info) <- paste0(ind_info$pool, "_", gsub("-1", "", ind_info$barcode))
        
        pool_list <- unique(ind_info$pool)
        txs_info_df <- lapply(pool_list, function(sample_id){
            tmp_file <- sprintf("%s/%s/UMItools.tsv.rds", 
                                bam_run_dir, sample_id)
            if(file.exists(tmp_file)){
                rawcount <- readRDS(tmp_file)
                colnames(rawcount) <- paste0(sample_id, "_", colnames(rawcount))
                if(sum(colnames(rawcount)%in%rownames(ind_info)) > 1){
                    sub_txs_sce <- rawcount[, colnames(rawcount)%in%rownames(ind_info)]
                    utrome_row_pos <- sub_txs_sce@i+1
                    utrome_col_pos <- findInterval(seq(sub_txs_sce@x)-1, sub_txs_sce@p[-1])+1
                    utrome_val <- sub_txs_sce@x
                    sub_txs_df <- data.frame("cell"=sub_txs_sce@Dimnames[[2]][utrome_col_pos],
                                             "PAS"=sub_txs_sce@Dimnames[[1]][utrome_row_pos],
                                             "PAS_rawcount"=utrome_val)
                    return(sub_txs_df)
                }
                
            }
        }) %>% bind_rows()
        
        ##
        fun_ct_ratio <- function(ct){
            sub_ind_info <- subset(ind_info, cell_type == ct)
            dir.create(paste0(ct_out_dir, "/", ct))
            cell_names <- rownames(sub_ind_info)
            # APA
            out_df <- subset(txs_info_df, cell%in%cell_names)
            multinames <- out_df$PAS[grepl(",", out_df$PAS)]
            rename_PAS <- sapply(multinames, function(x){
                if(grepl(",", x)){
                    x_list <- strsplit(x, split=",")[[1]]
                    sub_anno <- subset(PAS_anno, site_name%in%x_list)
                    if(sub_anno$strand[1] == "+"){
                        return(sub_anno$site_name[which.max(sub_anno$site_position)])
                    }else{
                        return(sub_anno$site_name[which.min(sub_anno$site_position)])
                    }
                }else{
                    return(x)
                }
            }) %>% unlist()
            out_df$PAS[grepl(",", out_df$PAS)] <- rename_PAS
            ## remove NPS
            out_df <- out_df[!grepl(":NPS", out_df$PAS), ]
            out_df <- out_df %>% group_by(cell, PAS) %>% summarise(PAS_count=sum(PAS_rawcount)) %>% as.data.frame()
            write.table(out_df, paste0(ct_out_dir, "/", ct, "/cellquant_", ind_id, ".txt"), row.names = F, quote = F, sep = "\t")
            
            out_df$UTR <- gsub(":PAS[0-9]+|:NPS[0-9]+", "", out_df$PAS)
            out_df$gene <- gsub(":I[0-9]+", "", out_df$UTR)
            
            ## pseudobulk
            pseudobulk_df <- out_df %>%
                group_by(PAS, UTR) %>%
                summarise(PAS_sum = sum(PAS_count), support=length(unique(cell))) %>%
                group_by(UTR) %>%
                mutate(UTR_sum = sum(PAS_sum)) %>%
                mutate(ratio = PAS_sum/UTR_sum, Ind=ind_id) %>%
                select(Ind, Name=PAS, support, PAS_sum, UTR, UTR_sum, ratio) %>% 
                unique()
            PAS_index <- pseudobulk_df$support>=3 & pseudobulk_df$PAS_sum>=5 & pseudobulk_df$UTR_sum >= 10
            table(PAS_index)
            pseudobulk_df <- pseudobulk_df[PAS_index, ]
            
            write.table(unique(pseudobulk_df$Name), paste0(ct_out_dir, "/", ct, "/confident_", ind_id, ".txt"),
                        row.names = F, col.names = F, quote = F, sep = "\t")
            return(paste0(ct_out_dir, "/", ct, "/cellquant_", ind_id, ".txt"))
        }
        
        cell_types <- unique(ind_info$cell_type)
        return_list <- mclapply(cell_types, fun_ct_ratio, mc.cores = 2)
        unlist(return_list)
    })
    
    ## find confident PASs
    system(command = sprintf("find %s -name \"confident*\" | xargs cat > %s/confident_PAS.txt", ct_out_dir, args_$outdir))
    system(command = sprintf("awk '{sum[$1]+=1}END{for(i in sum)print i\"\t\"sum[i]}' %s/confident_PAS.txt > %s/confident_PAS_summary.txt", 
                             args_$outdir, args_$outdir))
    
    
    ## PAS selection based on data
    PAS_info <- read.table(sprintf("%s/confident_PAS_summary.txt", args_$outdir), sep = "\t") %>% as.data.frame()
    confident_PAS <- PAS_info$V1[PAS_info$V2 >= min(50, max(PAS_info$V2)*0.05)] %>% unique()
    
    ## PASs
    raw_anno <- read.table(sprintf("%s/Annotation/KeepPAS_data_annotation.txt", args_$outdir), header = T) 
    confident_anno <- subset(raw_anno, site_name%in%confident_PAS)
    # write.table(confident_anno, "2023-10-27-quantification/molecular_phenotypes/gene_quantification/00confident_36177_PASs.txt", row.names = F, sep = "\t", quote = F)
    
    ## TUTR and ALE
    TUTR_isoform <- names(table(confident_anno$LE_name))[table(confident_anno$LE_name) > 1]
    TUTR_anno <- subset(confident_anno, LE_name%in%TUTR_isoform)
    table(TUTR_anno$LE_type)
    
    UTR_anno <- confident_anno %>% group_by(LE_name) %>% summarise(chr=unique(chr), 
                                                                   site_position=case_when(strand[1]=="+"~unique(LE_start),
                                                                                           strand[1]=="-"~unique(LE_end)),
                                                                   site_name=unique(LE_name), 
                                                                   strand=unique(strand), site_type="-", LE_start=min(LE_start),
                                                                   LE_end=max(LE_end), LE_name=unique(gene_name), 
                                                                   LE_type=unique(LE_type), gene_name=unique(gene_name), gene_type=unique(gene_type))
    UTR_anno <- UTR_anno[, intersect(colnames(confident_anno), colnames(UTR_anno))]
    
    ALE_gene <- names(table(UTR_anno$LE_name))[table(UTR_anno$LE_name) > 1]
    if(length(ALE_gene) >=1 ){
        ALE_anno <- subset(UTR_anno, LE_name%in%ALE_gene)
        ALE_anno$LE_type <- "ALE"
        
        ALE_PAS <- subset(confident_anno[, colnames(UTR_anno)], LE_name%in%ALE_anno$site_name)
        ALE_PAS$LE_type <- "ALE"
        ALE_PAS <- subset(ALE_PAS, !site_name%in%TUTR_anno$site_name)
        
        clean_anno <- rbind(TUTR_anno[, colnames(UTR_anno)], ALE_anno, ALE_PAS)
    }else{
        clean_anno <- rbind(TUTR_anno[, colnames(UTR_anno)])
    }
    
    clean_anno_addpd <- clean_anno %>% group_by(LE_name) %>% 
        mutate(mean_dist=mean(site_position)) %>% 
        mutate(pd_type=case_when("+"%in%strand & site_position <= mean_dist ~ "proximal",  
                                 "+"%in%strand & site_position > mean_dist ~ "distal",
                                 "-"%in%strand & site_position < mean_dist ~ "distal",  
                                 "-"%in%strand & site_position >= mean_dist ~ "proximal",
                                 .default = "-"))
    clean_anno_addpd$site_type <- clean_anno_addpd$pd_type
    clean_anno_addpd <- clean_anno_addpd[, !colnames(clean_anno_addpd)%in%c("mean_dist", "pd_type")]
    
    nrow(clean_anno_addpd)
    table(grepl("PAS", clean_anno_addpd$site_name))
    
    write.table(clean_anno_addpd, sprintf("%s/Annotation/Selected_PAS_sites.txt", args_$outdir), row.names = F, sep = "\t", quote = F)
    
    clean_anno_bed <- clean_anno_addpd %>% as.data.frame() %>% filter(grepl("PAS", site_name)) %>% 
        mutate(start=site_position-1, end=site_position) %>% 
        select(chr, start, end, site_name, gene_type, strand) %>% arrange(chr, start)
    write.table(clean_anno_bed, sprintf("%s/Annotation/Selected_PAS_sites.bed", args_$outdir), row.names = F, col.names = F, sep = "\t", quote = F)
    
    
    ## APA quantification
    pbmclapply(unlist(cellquant_list), function(x){
        x_info <- strsplit(x, "/")[[1]]
        ct <- x_info[4]
        ind_id <- gsub(".txt", "", gsub("cellquant_", "", x_info[5]))
        
        ##
        PAS_selection <- vroom::vroom(sprintf("%s/Annotation/Selected_PAS_sites.txt", args_$outdir)) %>% as.data.frame()
        PAS_selection <- PAS_selection %>% group_by(LE_name) %>% 
            mutate(order=case_when(strand=="+" ~ order(site_position), 
                                   strand=="-" ~ order(site_position, decreasing = T))) %>% 
            mutate(weight=(order-1)/(length(order)-1))
        PAS_selection <- as.data.frame(PAS_selection)
        rownames(PAS_selection) <- PAS_selection$site_name
        
        TUTR_PAS <- subset(PAS_selection, LE_type%in%"TUTR")
        multiPAS_name <- names(table(TUTR_PAS$LE_name))[table(TUTR_PAS$LE_name)>2]
        
        ALE_UTR <- subset(PAS_selection, LE_type%in%"ALE")
        multiALE_name <- names(table(ALE_UTR$LE_name))[table(ALE_UTR$LE_name)>2]
        
        
        ## load PAS quantificant in each cell
        raw_cellquant <- read.table(x, header = T)
        
        raw_cellquant$UTR <- gsub(":PAS[0-9]+|:NPS[0-9]+", "", raw_cellquant$PAS)
        raw_cellquant$gene <- gsub(":I[0-9]+", "", raw_cellquant$UTR)
        raw_cellquant <- subset(raw_cellquant, PAS%in%PAS_selection$site_name)
        
        
        ## write out cell type APA; cell subtyep APA; and pseudo-PBMC APA;
        ## write out variance APA
        each_cell_APA_quant <- function(fin_cellquant){
            out_dir_name <- gsub("selection", "quantification", x)
            dir.create(dirname(out_dir_name), recursive = T)
            ## each PAS of isoform with multi-PASs
            cell_merge_df1 <- fin_cellquant %>% filter(UTR%in%TUTR_PAS$LE_name) %>% 
                group_by(cell, PAS, UTR) %>%
                summarise(PAS_sum = sum(PAS_count), support=length(unique(cell))) %>% 
                group_by(cell, UTR) %>%
                mutate(UTR_sum = sum(PAS_sum)) %>%
                mutate(ratio = PAS_sum/UTR_sum, Ind=ind_id) %>% 
                select(Ind, Name=PAS, support, PAS_sum, UTR, UTR_sum, ratio) %>% 
                unique()
            
            ## PDUI
            cell_merge_df3 <- cell_merge_df1
            cell_merge_df3$weight <- PAS_selection[cell_merge_df3$Name, "weight"]
            cell_merge_df3$w_ratio <- cell_merge_df3$ratio * cell_merge_df3$weight
            cell_merge_df3 <- cell_merge_df3 %>% 
                group_by(cell, UTR, UTR_sum) %>% 
                summarise(ratio=sum(w_ratio), support=sum(support)) %>%
                mutate(Ind=ind_id, gene=UTR, gene_sum=UTR_sum) %>%
                select(Ind, Name=UTR, support, UTR_sum, gene, gene_sum, ratio) %>% unique()
            # cell_merge_df3 <- cell_merge_df3[cell_merge_df3$Name%in%confident_UTR, ]
            
            cell_merge_df3$Name <- paste0(cell_merge_df3$Name, ":PDUI")
            # plot(density(cell_merge_df3$ratio))
            
            ## PPUI
            cell_merge_df1 <- subset(cell_merge_df1, UTR%in%multiPAS_name)
            # cell_merge_df1 <- subset(cell_merge_df1, UTR%in%confident_UTR)
            cell_merge_df1$Name <- paste0(cell_merge_df1$Name, ":PPUI")
            
            ## PIUI
            cell_merge_df2 <- fin_cellquant %>% filter(gene%in%ALE_UTR$gene_name) %>% 
                group_by(cell, UTR, gene) %>%
                summarise(UTR_sum = sum(PAS_count), support=length(unique(cell))) %>%
                group_by(cell, gene) %>%
                mutate(gene_sum = sum(UTR_sum)) %>%
                mutate(ratio = UTR_sum/gene_sum, Ind=ind_id) %>%
                select(Ind, Name=UTR, support, UTR_sum, gene, gene_sum, ratio) %>% unique()
            
            ## PDIUI
            cell_merge_df4 <- cell_merge_df2
            cell_merge_df4$weight <- PAS_selection[cell_merge_df4$Name, "weight"]
            cell_merge_df4$w_ratio <- cell_merge_df4$ratio * cell_merge_df4$weight
            cell_merge_df4 <- cell_merge_df4 %>% 
                group_by(cell, gene, gene_sum) %>% 
                summarise(ratio=sum(w_ratio), support=sum(support)) %>%
                mutate(Ind=ind_id, Name=gene, UTR_sum=gene_sum) %>%
                select(Ind, Name, support, UTR_sum, gene, gene_sum, ratio) %>% unique()
            # cell_merge_df4 <- cell_merge_df4[cell_merge_df4$Name%in%confident_LE, ]
            cell_merge_df4$Name <- paste0(cell_merge_df4$Name, ":PDIUI")
            
            ## PIUI
            cell_merge_df2 <- subset(cell_merge_df2, gene%in%multiALE_name)
            # cell_merge_df2 <- subset(cell_merge_df2, gene%in%confident_LE)
            cell_merge_df2$Name <- paste0(cell_merge_df2$Name, ":PIUI")
            
            colnames(cell_merge_df2) <- colnames(cell_merge_df3) <- colnames(cell_merge_df4) <- colnames(cell_merge_df1)
            cell_merge_df <- rbind(cell_merge_df3, cell_merge_df2, cell_merge_df1, cell_merge_df4)
            
            write.table(cell_merge_df, out_dir_name, row.names = F, quote = F, sep = "\t")
        }
        
        each_cell_APA_quant(fin_cellquant = raw_cellquant)
        
        ## cell type
        celltype_APA_quant <- function(fin_cellquant){
            out_dir_name <- gsub("selection", "quantification", x)
            dir.create(dirname(out_dir_name), recursive = T)
            ## each PAS of isoform with multi-PASs
            bulk_merge_df1 <- fin_cellquant %>% dplyr::filter(UTR%in%TUTR_PAS$LE_name) %>% 
                group_by(PAS, UTR) %>%
                summarise(PAS_sum = sum(PAS_count), support=length(unique(cell))) %>% 
                group_by(UTR) %>%
                mutate(UTR_sum = sum(PAS_sum)) %>%
                mutate(ratio = PAS_sum/UTR_sum, Ind=ind_id) %>% 
                dplyr::select(Ind, Name=PAS, support, PAS_sum, UTR, UTR_sum, ratio) %>% 
                unique()
            
            confident_UTR <- subset(bulk_merge_df1, support >= 3 & PAS_sum >= 5 & UTR_sum >= 10) 
            confident_UTR <- unique(confident_UTR$UTR)
            
            ## PDUI
            bulk_merge_df3 <- bulk_merge_df1
            bulk_merge_df3$weight <- PAS_selection[bulk_merge_df3$Name, "weight"]
            bulk_merge_df3$w_ratio <- bulk_merge_df3$ratio * bulk_merge_df3$weight
            bulk_merge_df3 <- bulk_merge_df3 %>% 
                group_by(UTR, UTR_sum) %>%
                summarise(ratio=sum(w_ratio), support=sum(support)) %>%
                mutate(Ind=ind_id, gene=UTR, gene_sum=UTR_sum) %>%
                select(Ind, Name=UTR, support, UTR_sum, gene, gene_sum, ratio) %>% unique()
            bulk_merge_df3 <- bulk_merge_df3[bulk_merge_df3$Name%in%confident_UTR, ]
            bulk_merge_df3$Name <- paste0(bulk_merge_df3$Name, ":PDUI")
            # plot(density(bulk_merge_df3$ratio))
            
            ## PPUI
            bulk_merge_df1 <- subset(bulk_merge_df1, UTR%in%multiPAS_name)
            bulk_merge_df1 <- subset(bulk_merge_df1, UTR%in%confident_UTR)
            bulk_merge_df1$Name <- paste0(bulk_merge_df1$Name, ":PPUI")
            
            ## PIUI
            bulk_merge_df2 <- fin_cellquant %>% filter(gene%in%ALE_UTR$gene_name) %>% group_by(UTR, gene) %>%
                summarise(UTR_sum = sum(PAS_count), support=length(unique(cell))) %>%
                group_by(gene) %>%
                mutate(gene_sum = sum(UTR_sum)) %>%
                mutate(ratio = UTR_sum/gene_sum, Ind=ind_id) %>%
                select(Ind, Name=UTR, support, UTR_sum, gene, gene_sum, ratio) %>% unique()
            
            confident_LE <- subset(bulk_merge_df2, support >= 3 & UTR_sum >= 5 & gene_sum >= 10) 
            confident_LE <- unique(confident_LE$gene)
            
            ## PDIUI
            bulk_merge_df4 <- bulk_merge_df2
            bulk_merge_df4$weight <- PAS_selection[bulk_merge_df4$Name, "weight"]
            bulk_merge_df4$w_ratio <- bulk_merge_df4$ratio * bulk_merge_df4$weight
            bulk_merge_df4 <- bulk_merge_df4 %>% 
                group_by(gene, gene_sum) %>% 
                summarise(ratio=sum(w_ratio), support=sum(support)) %>%
                mutate(Ind=ind_id, Name=gene, UTR_sum=gene_sum) %>%
                select(Ind, Name, support, UTR_sum, gene, gene_sum, ratio) %>% unique()
            bulk_merge_df4 <- bulk_merge_df4[bulk_merge_df4$Name%in%confident_LE, ]
            bulk_merge_df4$Name <- paste0(bulk_merge_df4$Name, ":PDIUI")
            
            ## PIUI
            bulk_merge_df2 <- subset(bulk_merge_df2, gene%in%multiALE_name)
            bulk_merge_df2 <- subset(bulk_merge_df2, gene%in%confident_LE)
            bulk_merge_df2$Name <- paste0(bulk_merge_df2$Name, ":PIUI")
            
            colnames(bulk_merge_df2) <- colnames(bulk_merge_df3) <- colnames(bulk_merge_df4) <- colnames(bulk_merge_df1)
            bulk_merge_df <- rbind(bulk_merge_df3, bulk_merge_df2, bulk_merge_df1, bulk_merge_df4)
            
            write.table(bulk_merge_df, gsub("cellquant_", "pseudobulk_", out_dir_name), row.names = F, quote = F, sep = "\t")
        }
        
        celltype_APA_quant(fin_cellquant = raw_cellquant)
    })

    
    ## Generate rsd
    celltype_list <- list.files(sprintf("%s/PAS_quantification", args_$outdir), full.names = T)
    celltype_list <- celltype_list[sapply(celltype_list, dir.exists)]
    pbmcapply::pbmclapply(celltype_list, function(quant_input){
        this_cell_type <- basename(quant_input)
        
        ## generate ratio information
        ratio_prepare <- function(ind_id, celltype){
            ## load data
            apa_matrix_file <- sprintf("%s/pseudobulk_%s.txt", quant_input, ind_id)
            if(file.exists(apa_matrix_file)){
                apa_matrix <- read.table(apa_matrix_file, header = T)
                
                if(nrow(apa_matrix) > 0){
                    apa_df <- apa_matrix
                    ## Average ratio between replicates
                    tot_counts <- apa_df %>% dplyr::select(Name, Ind, ratio) %>% 
                        spread(Ind, ratio) %>% column_to_rownames('Name')
                    tot_counts <- apply(tot_counts, 1, function(x){mean(x[!is.na(x)])})
                    ## Average support between replicates
                    n_support <- apa_df %>% dplyr::select(Name, Ind, support) %>% 
                        spread(Ind, support) %>% column_to_rownames('Name')
                    n_support <- apply(n_support, 1, function(x){mean(x[!is.na(x)])})
                    ## Average coverage between replicates
                    n_coverage <- apa_df %>% dplyr::select(Name, Ind, PAS_sum) %>% 
                        spread(Ind, PAS_sum) %>% column_to_rownames('Name')
                    n_coverage <- apply(n_coverage, 1, function(x){mean(x[!is.na(x)])})
                    overname <- intersect(names(tot_counts), names(n_support))
                    ## New summary table
                    tmp_gene_name <- sapply(names(tot_counts[overname]), function(x){strsplit(x, split = ":")[[1]][1]})
                    cnt_df <- data.frame(counts=tot_counts[overname],
                                         n_support=n_support[overname],
                                         n_coverage=n_coverage[overname]) %>% 
                        rownames_to_column('ensembl') %>% 
                        mutate(cell_type = celltype,
                               individual_id = ind_id,
                               symbol = tmp_gene_name) %>%
                        dplyr::select(cell_type,individual_id,ensembl,symbol,counts,n_coverage,n_support)
                    return(cnt_df)
                }
            }
        }
        
        ind_info <- list.files(quant_input)
        ## pseudobulk
        pseudobulk_ind_info <- ind_info[grepl("pseudobulk_", ind_info)]
        pseudobulk_ind_info <- gsub(".txt", "", pseudobulk_ind_info)
        pseudobulk_ind_info <- gsub("pseudobulk_", "", pseudobulk_ind_info)
        ## pseudobulk of cells
        pseudobulk_apa <- pbmcapply::pbmclapply(pseudobulk_ind_info, ratio_prepare, celltype=this_cell_type,
                                   mc.cores = min(length(pseudobulk_ind_info),
                                                  args_$threads), mc.preschedule = FALSE) %>% bind_rows()
        saveRDS(pseudobulk_apa, sprintf('%s/pseudobulk_%s.rds', dirname(quant_input), this_cell_type))
    })
    
    ## Get indenpendent APA events
    data_dir <- "example/APA/PAS_quantification/"
    rds_list <- list.files(data_dir, pattern = ".rds")
    rds_list <- rds_list[grepl("pseudobulk_", rds_list)]
    
    loadRDS <- function(x){
        tmp_quant <- readRDS(paste0(data_dir, x))
        if(nrow(tmp_quant) > 0){
            src_name <- gsub(".rds", "", x)
            tmp_quant <- tmp_quant %>% mutate(source=src_name) %>% 
                mutate(source=paste0(source, "_", individual_id)) %>% 
                dplyr::select(ensembl, source, counts) %>% unique()
            tmp_quant
        }
    }
    PAS_quants <- mclapply(rds_list, loadRDS, mc.cores = 30, mc.preschedule = TRUE) %>% bind_rows()
    gc()
    ## 104,355,317 | 131,296,197 | 171,142,605
    
    PAS_quants_matrix_raw <- reshape2::dcast(PAS_quants, formula = ensembl~source, 
                                             value.var = "counts") %>% column_to_rownames('ensembl')
    gc()
    ## 26000 13089 | 35528 13646
    
    na1_lower_index <- apply(PAS_quants_matrix_raw, 1, function(x){sd(as.numeric(x), na.rm = T)!=0})
    table(na1_lower_index)
    rownames(PAS_quants_matrix_raw)[!na1_lower_index]
    PAS_quants_matrix <- PAS_quants_matrix_raw[which(na1_lower_index), ]
    gc()
    
    
    preimpute_na_summary <- apply(PAS_quants_matrix, 1, function(x){sum(!is.na(as.numeric(x)))})
    preimpute_na_index <- preimpute_na_summary >= min(50, 0.05*ncol(PAS_quants_matrix)) # 0.99*ncol(PAS_quants_matrix) # 0.05*ncol(PAS_quants_matrix)
    table(preimpute_na_index)
    PAS_quants_matrix <- PAS_quants_matrix[preimpute_na_index, ]
    
    preimpute_range_summary <- apply(PAS_quants_matrix, 1, function(x){
        quantile(x, probs = 0.75, na.rm = T)-quantile(x, probs = 0.25, na.rm = T)
    })
    preimpute_range_index <- preimpute_range_summary >= 0.02
    # table(preimpute_range_index)
    # rm_dd <- rownames(PAS_quants_matrix)[!preimpute_range_index]
    # rm_dd_type <- sapply(strsplit(rm_dd, ":"), function(x){x[length(x)]})
    # table(rm_dd_type)
    PAS_quants_matrix <- PAS_quants_matrix[preimpute_range_index, ]
    
    # APA_list <- read.table("2023-10-27-quantification/molecular_phenotypes/APA/APA_35528_events.txt", header = T)
    # ggvenn::ggvenn(list("A"=APA_list$site_name, "B"=rownames(PAS_quants_matrix)))
    
    ## filter 3: remove minor PAS
    no_cor_PAS_iter <- function(num_ind, cor_df) {
        # cat("Depth:", depth, "| num_ind:", num_ind, "\n")
        result <- c()
        remaining_ind <- num_ind
        while (length(remaining_ind) > 0) {
            current <- remaining_ind[1]
            related_cols <- which(cor_df[current, ])
            result <- c(result, current)
            remaining_ind <- remaining_ind[!remaining_ind %in% c(current, related_cols)]
        }
        return(result)
    }
    # no_cor_PAS_iter(num_ind = tmp_index, cor_df = cor_df)
    
    no_cor_PAS_list <- function(num_ind, cor_df) {
        # cat("Depth:", depth, "| num_ind:", num_ind, "\n")
        result <- c()
        remaining_ind <- num_ind
        while (length(remaining_ind) > 0) {
            current <- remaining_ind[1]
            related_cols <- which(cor_df[current, ])
            result <- c(result, paste0(rownames(cor_df)[remaining_ind[remaining_ind %in% c(current, related_cols)]], collapse = ","))
            remaining_ind <- remaining_ind[!remaining_ind %in% c(current, related_cols)]
        }
        return(paste0(result, collapse = "|"))
    }
    
    library(Hmisc)
    APA_PAS_select <- function(APA_index){
        APA_df <- PAS_data[ratio_matrix_group==APA_index,]
        # APA_df <- abs(APA_df)
        abu_order <- apply(APA_df, 1, function(x){mean(abs(x), na.rm=T)})
        type_order <- factor(sapply(strsplit(rownames(APA_df), ":"), function(x){x[length(x)]}), levels = rev(c("PDUI", "PPUI", "PDIUI", "PIUI")))
        abu_order <- order(type_order, abu_order, decreasing = T)
        APA_df <- APA_df[abu_order, ]
        if(nrow(APA_df) >= 2){
            cor_df_data <- rcorr(t(APA_df), type = "spearman")
            cor_df <- abs(cor_df_data$r) >= 0.3 & cor_df_data$P < 0.05
            cor_df[is.na(cor_df)] <- FALSE
            ##
            tmp_index <- 1:nrow(cor_df)
            names(tmp_index) <- rownames(cor_df)
            out_index <- no_cor_PAS_iter(num_ind = tmp_index, cor_df = cor_df)
            out_relation <- no_cor_PAS_list(num_ind = tmp_index, cor_df = cor_df)
            write.table(out_relation, "example/APA/Annotation/00PAS_cor.txt", 
                        quote = F, sep = "\t", append = T, row.names = F, col.names = F)
            return(names(out_index))
        }else{
            return(rownames(APA_df))
        }
    }
    
    # PAS_index <- grepl("PDUI|PPUI", rownames(PAS_quants_matrix))
    # table(PAS_index)
    PAS_data <- PAS_quants_matrix
    ratio_matrix_group <- sapply(strsplit(rownames(PAS_data), ":"), function(x){x[1]})
    use_PAS <-  unique(ratio_matrix_group)
    
    PAS_clean <- mclapply(use_PAS, APA_PAS_select, mc.cores = 20)
    PAS_clean <- unlist(PAS_clean) %>% unique()
    
    out_file <- "example/APA/Annotation/Independent_PAS_set.txt"
    write.table(PAS_clean, out_file, row.names = F, col.names = F, quote = F)
    
    dim(PAS_data)
    write.table(rownames(PAS_data), "example/APA/Annotation/Independent_used_pas_utr.txt", 
                row.names = F, col.names = F, quote = F)
    
    ## Generate unified annotation
    raw_anno <- read.table(sprintf("example/APA/Annotation/Selected_PAS_sites.txt"), header = T) 
    length(unique(raw_anno$gene_name))
    length(unique(raw_anno$gene_name[raw_anno$LE_type=="TUTR"]))
    
    TUTR_PAS <- subset(raw_anno, LE_type%in%"TUTR")
    multiPAS_name <- names(table(TUTR_PAS$LE_name))[table(TUTR_PAS$LE_name)>2]
    
    ALE_UTR <- subset(raw_anno, LE_type%in%"ALE")
    multiALE_name <- names(table(ALE_UTR$LE_name))[table(ALE_UTR$LE_name)>2]
    
    ##
    PAS_UTR_anno <- raw_anno %>% filter(grepl(":PAS", site_name)&LE_type=="TUTR") %>% 
        group_by(LE_name) %>% summarise(chr=unique(chr), 
                                        site_position=case_when(strand[1]=="+"~unique(LE_start),
                                                                strand[1]=="-"~unique(LE_end)),
                                        site_name=unique(LE_name), 
                                        strand=unique(strand), site_type="-", LE_start=min(LE_start),
                                        LE_end=max(LE_end), LE_name=unique(gene_name), 
                                        LE_type=unique(LE_type), gene_name=unique(gene_name), gene_type=unique(gene_type))
    PAS_UTR_anno <- PAS_UTR_anno[, colnames(raw_anno)]
    
    ##
    multi_PAS_anno <- raw_anno[grepl(":PAS", raw_anno$site_name)&raw_anno$LE_type=="TUTR"&raw_anno$LE_name%in%multiPAS_name, ]
    
    ##
    # PAS_UTR_anno <- subset(PAS_UTR_anno, !site_name%in%raw_anno$site_name)
    # LE_UTR_anno <- raw_anno[!grepl(":PAS", raw_anno$site_name), ]
    # UTR_anno <- rbind(PAS_UTR_anno, LE_UTR_anno)
    
    LE_UTR_anno <- raw_anno[(!grepl(":PAS", raw_anno$site_name))&raw_anno$LE_type=="ALE", ]
    # ALE_genes <- unique(LE_UTR_anno$LE_name)
    # LE_PAS_UTR_anno <- subset(PAS_UTR_anno, LE_name%in%ALE_genes)
    # LE_PAS_UTR_anno <- subset(LE_PAS_UTR_anno, !site_name%in%LE_UTR_anno$site_name)
    
    ## 
    multi_PAS_anno$site_name <- paste0(multi_PAS_anno$site_name, ":PPUI")
    PAS_UTR_anno$site_name <- paste0(PAS_UTR_anno$site_name, ":PDUI")
   
    if(nrow(LE_UTR_anno) > 0){
        #
        ALE_anno <- LE_UTR_anno %>% group_by(LE_name) %>% summarise(chr=unique(chr),
                                                                    site_position=case_when(strand[1]=="+"~min(LE_start),
                                                                                            strand[1]=="-"~max(LE_end)),
                                                                    site_name=unique(LE_name),
                                                                    strand=unique(strand), site_type="-", LE_start=min(LE_start),
                                                                    LE_end=max(LE_end), LE_name=unique(gene_name),
                                                                    LE_type="-", gene_name=unique(gene_name),
                                                                    gene_type=unique(gene_type))
        ALE_anno <- ALE_anno[, colnames(raw_anno)]
        ALE_anno$site_name <- paste0(ALE_anno$site_name, ":PDIUI")
        
        table(table(LE_UTR_anno$LE_name))
        multi_LE_UTR_anno <- subset(LE_UTR_anno, LE_name%in%multiALE_name)
        multi_LE_UTR_anno$site_name <- paste0(multi_LE_UTR_anno$site_name, ":PIUI")
        
        merge_anno <- rbind(multi_PAS_anno, PAS_UTR_anno, multi_LE_UTR_anno, ALE_anno)
        
    } else{
        merge_anno <- rbind(multi_PAS_anno, PAS_UTR_anno) 
    }
    #
    clean_anno <- unique(merge_anno)
    
    nrow(clean_anno)
    
    #
    filtered_anno <- read.table("example/APA/Annotation/Independent_PAS_set.txt")
    
    # ggvenn::ggvenn(list("A"=filtered_anno$V1, "B"=clean_anno$site_name))
    out_clean_anno <- subset(clean_anno, site_name%in%filtered_anno$V1)
    
    table(table(out_clean_anno$site_name))
    write.table(out_clean_anno, "example/APA/Annotation/Independent_APA_events.txt", row.names = F, sep = "\t", quote = F)
    
    ## Generate matrix 
    rds_list <- list.files("example/APA/PAS_quantification/", pattern = "rds", full.names = T)
    dir.create("example/APA/cell_type_quant")
    pbmcapply::pbmclapply(rds_list, function(rds_file){
        clean_anno <- out_clean_anno
        ## PAS annotation and genomic location
        phenotype_using <- apply(clean_anno, 1, function(x){
            return(c(x[1], as.numeric(x[2])-1, x[2], x[3], x[8], x[4]))
        }) %>% t() %>% as.data.frame()
        colnames(phenotype_using) <- c("seqnames", "start", "end", "Name", "GeneName", "strand")
        phenotype_using <- phenotype_using %>% arrange(seqnames, as.numeric(start))
        
        ## generate ratio matrix
        sum_apa_df <- readRDS(rds_file)
        if(nrow(sum_apa_df)>0){
            ratio_matrix_raw <- sum_apa_df %>% filter(ensembl%in%phenotype_using$Name)
            
            ratio_matrix_raw <- ratio_matrix_raw %>% 
                dplyr::select(individual_id, ensembl, counts) %>% 
                unique() %>% 
                spread(individual_id, counts) %>% 
                column_to_rownames('ensembl')
            
            out_table <- inner_join(phenotype_using, ratio_matrix_raw %>% rownames_to_column("Name"), by="Name")
            write.table(out_table, file = sprintf("example/APA/cell_type_quant/%s.txt", 
                                                  gsub("pseudobulk_|.rds", "", basename(rds_file))), row.names = F, sep = "\t", quote = F) 
        }
    })
}
