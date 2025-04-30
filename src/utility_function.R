
load_bam_cell_anno3 <- function(bam_anno_file_="example/cell_annotation.csv",
                                bam_file_ = "example/bam_list.csv",
                                check_bam=FALSE){
  
  # pooled bam-cell-donor- details:
  bam_cell_anno <- fread(bam_anno_file_);
  # names(bam_cell_anno)[which(names(bam_cell_anno)=="POOL")] <- "pool"
  uniq_ID <- unique(bam_cell_anno[,.(cell_type, donor, pool)])
  uniq_ID <- cbind(data.table(split_id = paste0("splited_",1:nrow(uniq_ID))), uniq_ID)
  bam_cell_anno <- merge(bam_cell_anno, uniq_ID, by =c("cell_type","donor", "pool"))
  # check whether cells from the same cell type belong to one individuals are split into different bams:
  # print(unique(bam_cell_anno[,.(cell_type, donor, pool)])[,.(num_bam = nrow(.SD)), by =c("cell_type", "donor")][order(num_bam)])
  
  # pooled bam file path:
  bam_anno <- fread(bam_file_);
  if(check_bam){
    bam_anno <- bam_anno[pool %in% unique(bam_cell_anno$pool)]
    if(nrow(bam_anno) == 0){
        stop(" No overlap between cell annotation and bam file !")
    }
  }
  return(list(bam_cell_anno=bam_cell_anno, bam_anno=bam_anno))
}


remove_dups_pools2 <- function(bam_anno_, 
                               parThreads_=30, out_path ="",
                               log_=log_file){
  start_time_i <- Sys.time()
  message("==> total: ",nrow(bam_anno_))
  bam_anno_[,c("bam_dedup_path"):=.( str_replace_all(bamfile, "possorted_genome_bam.bam", "possorted_genome_dedup.bam") )]
  bam_anno_$dedup_exist <- unlist(lapply(1:nrow(bam_anno_),function(i){
    bam_dedup_path_i <- bam_anno_[i,]$bam_dedup_path
    if( !file.exists(bam_dedup_path_i) | !file.exists(paste0(bam_dedup_path_i,".bai")) | ifelse(is.na(file.size(bam_dedup_path_i)), 0, file.size(bam_dedup_path_i)) == 0  ){
      return(0)
    }else{
      return(1)
    }
  }))
  if(sum(bam_anno_$dedup_exist==0) > 0){
    do_bam_anno_ <- bam_anno_[dedup_exist==0,]
    message("==> left: ", nrow(bam_anno_))
    
    pbmcapply::pbmclapply(1:nrow(do_bam_anno_), function(i){
      message("==> ", i, " | ", do_bam_anno_[i,]$bamfile, " | ",do_bam_anno_[i,]$bam_dedup_path, " | ",date())
      bam_path_i <- do_bam_anno_[i,]$bamfile
      bam_dedup_path_i <- do_bam_anno_[i,]$bam_dedup_path
      if(!dir.exists(dirname(bam_dedup_path_i))){ dir.create(dirname(bam_dedup_path_i), recursive = TRUE) }
      if(!dir.exists(dirname(bam_dedup_path_i))){stop("outPath can not be created");}
      
      # start dedup:
      if( !file.exists(bam_dedup_path_i) | ifelse(is.na(file.size(bam_dedup_path_i)), 0, file.size(bam_dedup_path_i)) == 0  ){
        if(file.exists(bam_dedup_path_i)){ file.remove(bam_dedup_path_i) }
        system(paste0("umi_tools dedup -I ", bam_path_i, " -S ", bam_dedup_path_i, " --per-cell --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB "))
      }
      # create index:
      if(file.exists(bam_dedup_path_i) & !file.exists(paste0(bam_dedup_path_i,".bai"))){
        system(paste0("samtools index -@ ",4," ", bam_dedup_path_i, " "))
      }
    }, mc.cores = min(nrow(do_bam_anno_), parThreads_))
    end_time_i <- Sys.time()
    running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
    message("remove_dups_pools2 running time: ", running_time)
    fwrite(data.table(type = "remove_dups_pools2", time=running_time), log_, sep="\t", append = FALSE)
  }
  return(bam_anno_[,.(pool, bamfile, bam_dedup_path)])
}