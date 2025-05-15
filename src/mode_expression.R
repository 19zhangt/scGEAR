
# Load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(pbmcapply))
suppressMessages(library(Seurat))

#
run_expression_analysis <- function(bam_data_=bam_data, args_=args, log_=log_file) {
    #' Run complete expression analysis pipeline
    #'
    #' @param bam_data_ BAM files and Cell type annotation data
    #' @param args_ parameters
    #' @param log_ log file
    matrix_index <- split(bam_data_$bam_anno, bam_data_$bam_anno$pool)
    ##
    submetadata <- as.data.frame(bam_data_$bam_cell_anno)
    rownames(submetadata) <- paste0(submetadata$pool, "_", gsub("-[0-9]+", "", submetadata$barcode))

    data.list <- pbmcapply::pbmclapply(matrix_index, function(f_matrix_index){
        f_path <- sprintf(sprintf("%s/filtered_feature_bc_matrix", dirname(f_matrix_index$bam_dedup_path)))
        pool_data <- Read10X(f_path)
        colnames(pool_data) <- gsub("-[0-9]+", "", colnames(pool_data))
        colnames(pool_data) <- paste0(f_matrix_index$pool, "_", colnames(pool_data))
        col_index <- colnames(pool_data)%in%rownames(submetadata)
        table(col_index)
        pool_data <- pool_data[, col_index]
        return_data <- CreateSeuratObject(counts = pool_data)
        return(return_data)
    }, mc.cores = min(as.numeric(args_$threads), length(matrix_index)))

    ##
    merged.object <- base::merge(data.list[[1]], data.list[2:length(data.list)])
    merged.object <- JoinLayers(merged.object)
    tmp_pool_sm <- merged.object[['RNA']]$counts
    
    table(colnames(tmp_pool_sm)%in%rownames(submetadata))
    tmp_metadata <- submetadata[colnames(tmp_pool_sm), ]
    tmp_seurat <- CreateSeuratObject(tmp_pool_sm, meta.data = tmp_metadata)
    
    ## split seurat by cell type
    ct_names <- unique(tmp_seurat@meta.data$cell_type)
    tmp_seurat <- NormalizeData(tmp_seurat)
    
    ## Gene annotation and genomic location
    clean_anno <- read.table("reference/gencode.annotation.txt", header = F, sep = "\t")
    clean_anno <- subset(clean_anno, !V5%in%names(table(clean_anno$V5))[table(clean_anno$V5)==2])
    clean_anno$len <- clean_anno$V3-clean_anno$V2+1
    rownames(clean_anno) <- clean_anno$V5
    
    ##
    phenotype_using <- apply(clean_anno, 1, function(x){
        if(x[4]=="+"){
            return(c(x[1], as.numeric(x[2])-1, x[2], x[5], x[6], x[4]))
        }else{
            return(c(x[1], as.numeric(x[3])-1, x[3], x[5], x[6], x[4]))
        }
    }) %>% t() %>% as.data.frame()
    colnames(phenotype_using) <- c("seqnames", "start", "end", "Name", "GeneName", "strand")
    
    dir.create(sprintf("%s/cell_type_quant", args_$outdir), showWarnings = F)
    
    ## reshape expression data
    pbmcapply::pbmclapply(ct_names, function(f_ct_name){
        ct_seurat <- subset(tmp_seurat, cell_type==f_ct_name)
        ind_count <- unlist(table(ct_seurat@meta.data$donor))
        keep_list <- names(which(ind_count >= 10))
        if(length(keep_list) > 0){
            ct_seurat <- subset(x = ct_seurat, subset = donor %in% keep_list)
            print(paste0(length(keep_list)," individuals remain after excluding those with < 10 cells."))
            
            ct_seurat <- ScaleData(ct_seurat, vars.to.regress = 'pool', do.scale = FALSE, do.center = FALSE)
            expressed_cells <- rowSums(ct_seurat@assays$RNA$counts>0)
            expressed_genes <- names(expressed_cells)[expressed_cells/ncol(ct_seurat) > 0.01]
            sub_ct_seurat <- subset(ct_seurat, features = expressed_genes)
            
            bulk_inds <- AggregateExpression(sub_ct_seurat, group.by = "donor", return.seurat = TRUE)
            bulk_inds <- as(bulk_inds@assays$RNA$data, "dgCMatrix")
            
            utrome_row_pos <- bulk_inds@i+1
            utrome_col_pos <- findInterval(seq(bulk_inds@x)-1, bulk_inds@p[-1])+1
            utrome_val <- bulk_inds@x
            sub_txs_df <- data.frame("ind"=bulk_inds@Dimnames[[2]][utrome_col_pos],
                                     "gene"=bulk_inds@Dimnames[[1]][utrome_row_pos],
                                     "value"=utrome_val)
            
            sub_txs_df <- subset(sub_txs_df, gene%in%phenotype_using$Name)
            
            ratio_matrix_raw <- sub_txs_df %>% 
                dplyr::select(ind, gene, value) %>% 
                unique() %>% 
                spread(ind, value) %>% 
                column_to_rownames('gene')
            ratio_matrix_raw[is.na(ratio_matrix_raw)] <- 0
            
            ratio_matrix <- ratio_matrix_raw %>% as.data.frame() %>%
                rownames_to_column('Name') %>%
                inner_join(., phenotype_using, by='Name')
           
            ##
            ratio_matrix <- ratio_matrix %>%
                dplyr::select(-GeneName,-strand) %>%
                dplyr::select(seqnames,start,end,pid=Name,everything()) %>%
                arrange(seqnames, start)
            ratio_matrix[,1] <- gsub("chr", "", ratio_matrix[,1])
            ratio_matrix <- ratio_matrix %>% arrange(seqnames, as.numeric(start))
            colnames(ratio_matrix)[1:4] <- c('#chr','start','end','phenotype_id')
            
            write.table(ratio_matrix, sprintf("%s/cell_type_quant/%s.txt", args_$outdir, f_ct_name), row.names = F, sep = "\t", quote = F)
        }
    }, mc.cores = min(length(ct_names), args_$threads), mc.preschedule = F)
}
