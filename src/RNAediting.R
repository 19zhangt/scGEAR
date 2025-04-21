
split_dedup_of_pool_bam <- function(bam_data_, out_dir_ = "", 
                                    parThreads_=6, subThreads_=10, 
                                    sinto_exe_="sinto" ){
  if(!dir.exists(out_dir_)){ stop("==> |",out_dir_," | can not be created")}
 
  bam_anno_ <- bam_data_$bam_anno[file.exists(bam_dedup_path)]
  if(nrow(bam_anno_)==0){return()}
  bam_anno_$split_exist <- unlist(lapply(1:nrow(bam_anno_),function(i){
    bam_cell_anno_i <- bam_data_$bam_cell_anno[pool == bam_anno_[i,]$pool,]
    out_dir_i <- paste0(out_dir_, "/", bam_anno_[i,]$pool)
    if(length(setdiff(paste0(unique(bam_cell_anno_i$split_id), ".bam"), list.files(out_dir_i))) == 0){
      message(" ==> ",i," | All bams has been extracted! next!")
      return(1)
    }else{
      return(0)
    }
  }))
  
  if(sum(bam_anno_$split_exist==0) > 0){
    do_bam_anno_ <- bam_anno_[split_exist==0,]
    message("==> left: ",nrow(do_bam_anno_))
    if(nrow(do_bam_anno_)==0){return(0)}
    
    # 6 batches, each 12 threads:
    start_time_i <- Sys.time()
    run_time <- pbmcapply::pbmclapply( 1:nrow(do_bam_anno_), function(i){
      message("==> | ", i, " | ", do_bam_anno_[i,]$pool, " | ",date())
      out_dir_i <- paste0(out_dir_, "/", do_bam_anno_[i,]$pool)
      cell_id_path_i <- paste0(out_dir_, "/", do_bam_anno_[i,]$pool, "/bam_cell_anno.txt")
      if(!dir.exists(out_dir_i)){ dir.create(out_dir_i, recursive = TRUE) }
      bam_cell_anno_i <- bam_data_$bam_cell_anno[pool == do_bam_anno_[i,]$pool,]
      # write out
      fwrite(bam_cell_anno_i[,.(barcode, split_id)], cell_id_path_i, sep="\t", col.names = FALSE)
      # chekc output files exist:
      if(length(setdiff(paste0(unique(bam_cell_anno_i$split_id), ".bam"), list.files(out_dir_i))) == 0){
        message("All bams has been extracted! next!")
        return()
      }
      bam_dedup_path_i <- do_bam_anno_[i,]$bam_dedup_path
      if(!file.exists(bam_dedup_path_i) | !file.exists(paste0(bam_dedup_path_i,".bai")) | ifelse(is.na(file.size(bam_dedup_path_i)), 0, file.size(bam_dedup_path_i)) == 0 ){ return(0) }
      # split bam into separate files:
      system(paste0("sinto filterbarcodes -p ", subThreads_," -b ", bam_dedup_path_i, " -c ", cell_id_path_i, " --outdir ", out_dir_i))
      message("==> ",i, " | ", do_bam_anno_[i,]$pool, " | done! ",date())
      return()
    }, mc.cores = min(nrow(do_bam_anno_), parThreads_))
    end_time_i <- Sys.time()
    running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
    return(unlist(running_time))
  }
}


combine_bam_ind <- function(bam_cell_anno_=bam_cell_anno, 
                            bam_splited_dir_ ="", 
                            out_dir_ = "/", 
                            parThreads_=5, subThreads_=8){
  if(!dir.exists(out_dir_)){ dir.create(out_dir_, recursive = TRUE)}
  if(!dir.exists(out_dir_)){stop("==> ",out_dir_, " | can not be created")}
  inds <- unique(bam_cell_anno[,.(donor)])
  
  start_time_i <- Sys.time()
  tmp_out <- pbmcapply::pbmclapply(1:nrow(inds), function(i){
    message("==> i:",i," | ",date())
    ind_i <- inds[i,]$donor
    bam_cell_anno_i <- unique(bam_cell_anno_[donor==ind_i,.(cell_type, donor, pool, split_id)])
    bam_cell_anno_i$bam_splited_path <- paste0(bam_splited_dir_, "/", bam_cell_anno_i$pool ,"/",bam_cell_anno_i$split_id,".bam" )
    # if(!all(file.exists(bam_cell_anno_i$bam_splited_path))){stop("file does not exist!")}
    bam_cell_anno_i <- bam_cell_anno_i[file.exists(bam_cell_anno_i$bam_splited_path)]
    if(nrow(bam_cell_anno_i)==0){return()}
    bam_ind_out_path_i <-  paste0(out_dir_, "/", ind_i,".bam")
    bam_cell_anno_i$cluster_info <- paste0(bam_splited_dir_, "/", bam_cell_anno_i$pool, "/bam_cell_anno.txt" )
    cluster_info_i <- rbindlist(lapply(1:nrow(bam_cell_anno_i), function(k){ 
      a <- fread(bam_cell_anno_i[k]$cluster_info, header = FALSE);
      names(a)<-c("barcode", "split_id");
      a$pool <- bam_cell_anno_i[k]$pool;
      return(a[split_id==bam_cell_anno_i[k]$split_id]) 
    }))
    cluster_info_out_paht_i <- paste0(out_dir_, "/", inds[i,]$donor, "_cell_info.txt")
    fwrite(cluster_info_i, cluster_info_out_paht_i, sep="\t")
    
    # merge bams:
    if( !file.exists(bam_ind_out_path_i) | ifelse(is.na(file.size(bam_ind_out_path_i)), 0, file.size(bam_ind_out_path_i)) == 0 ){
      message(" ==> File does not exist: ",ind_i, " | start merge: ")
      if(nrow(bam_cell_anno_i)>1){
        system(paste0("samtools merge -f -@ ",subThreads_," -O bam ", bam_ind_out_path_i, " ",
                      paste0(bam_cell_anno_i$bam_splited_path, collapse = " ")))
      }else{
        system(paste0("cp ",bam_cell_anno_i[1,]$bam_splited_path, " ", bam_ind_out_path_i))
      }
    }
    # make index for merged bam file:
    if( file.exists(bam_ind_out_path_i) & !file.exists( paste0(bam_ind_out_path_i, ".bai") )){
      message( " ==> ", ind_i, " | create index:")
      system(paste0("samtools index -@ ",subThreads_," ", bam_ind_out_path_i, " &"))
    }
  }, mc.cores = min(nrow(inds), parThreads_))
  end_time_i <- Sys.time()
  running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
  return(running_time)
}


combine_bam_celltype_ind <- function(bam_cell_anno_=bam_cell_anno, 
                                     bam_splited_dir_ = "", 
                                     out_dir_ = "", 
                                     parThreads_=15, 
                                     subThreads_=4){
  start_time_i <- Sys.time()
  if(!dir.exists(out_dir_)){dir.create(out_dir_, recursive = TRUE)}
  if(!dir.exists(out_dir_)){ stop("==> ",out_dir_, " | can not be created.." ) }
  cell_types <- unique(bam_cell_anno_[,.(cell_type)])
  for(i in 1:nrow(cell_types)){
    message("==> i:",i," | ",date())
    bam_cell_anno_i <- unique(bam_cell_anno_[cell_type==cell_types[i,]$cell_type,.(cell_type, donor, pool, split_id)])
    inds_i <- unique(bam_cell_anno_i[,.(donor)])
    # parallel:
    tmp_out <- pbmcapply::pbmclapply(1:nrow(inds_i), function(j){
      bam_ind_out_path_j <-  paste0(out_dir_,"/", cell_types[i,]$cell_type,"/",inds_i[j,]$donor,".bam")
      cluster_info_out_paht_j <- paste0(out_dir_,"/", cell_types[i,]$cell_type,"/",inds_i[j,]$donor,"_cell_info.txt")
      bam_ind_dedup_out_path_j <-  paste0(out_dir_,"/", cell_types[i,]$cell_type,"/",inds_i[j,]$donor,"_dedup.bam")
      if(!dir.exists(dirname(bam_ind_out_path_j))){ dir.create(dirname(bam_ind_out_path_j), recursive = TRUE) }
      bam_cell_anno_j <- unique(bam_cell_anno_i[donor==inds_i[j,]$donor,.(donor, cell_type, pool, split_id)])
      bam_cell_anno_j$bam_splited_path <- paste0(bam_splited_dir_,"/", bam_cell_anno_j$pool, "/", bam_cell_anno_j$split_id, ".bam")
      bam_cell_anno_j$cluster_info <- paste0(bam_splited_dir_,"/", bam_cell_anno_j$pool, "/bam_cell_anno.txt" )
      # if(!all(file.exists(bam_cell_anno_j$bam_splited_path))){stop("bam file does not exist!")}
      # if(!all(file.exists(bam_cell_anno_j$cluster_info))){stop("splited info file does not exist!")}
      bam_cell_anno_j <- bam_cell_anno_j[file.exists(bam_cell_anno_j$bam_splited_path)]
      if(nrow(bam_cell_anno_j)==0){ return() }
      cluster_info_j <- rbindlist(lapply(1:nrow(bam_cell_anno_j), function(k){ 
        a <- fread(bam_cell_anno_j[k]$cluster_info, header = FALSE);
        names(a)<-c("barcode", "split_id");
        a$pool <- bam_cell_anno_j[k]$pool;
        return(a[split_id==bam_cell_anno_j[k]$split_id]) 
      }))
      fwrite(cluster_info_j, cluster_info_out_paht_j, sep="\t")
      # merge bams:
      # edit by Lihai Gong, add TRUE in if sequence to make sure new file will be re-generated
      # if(TRUE | !file.exists(bam_ind_out_path_j) | ifelse(is.na(file.size(bam_ind_out_path_j)), 0, file.size(bam_ind_out_path_j)) == 0 ){
      if(!file.exists(bam_ind_out_path_j) | ifelse(is.na(file.size(bam_ind_out_path_j)), 0, file.size(bam_ind_out_path_j)) == 0 ){
        message(" ==> File does not exist: ", cell_types[i,]$cell_type, " | "," | ", inds_i[j,]$donor, " | start merge: ")
        if(nrow(bam_cell_anno_j)>1){
          system(paste0("samtools merge -f -@ ",subThreads_," -O bam ", bam_ind_out_path_j, " ",
                        paste0(bam_cell_anno_j$bam_splited_path, collapse = " ")))
        }else{
          system(paste0("ln -s ", normalizePath(bam_cell_anno_j[1,]$bam_splited_path), " ", bam_ind_out_path_j))
        }
      }
      # make index for merged bam file:
      if( file.exists(bam_ind_out_path_j) & !file.exists( paste0(bam_ind_out_path_j, ".bai") )){
        message( " ==>", " - ",cell_types[i,]$cell_type,"  - ",inds_i[j,]$donor, " | create index:")
        system(paste0("samtools index -@ ",subThreads_," ", bam_ind_out_path_j, " &"))
      }
      return(j)
    }, mc.cores = min(nrow(inds_i), parThreads_))
  }
  end_time_i <- Sys.time()
  running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
  return(running_time)
}

# -t threads
# -c Min. read coverage [10]		
# -q Min. quality score [25]
# -m Min. mapping quality score [25]
# -O Min. homoplymeric length [5]
# -s 1 read1 as RNA,read2 not as RNA;
# -E Exclude positions with multiple changes
# -v Minimum number of reads supporting the variation [3]
# -l List of known RNA editing events
# -n Minimum editing frequency [0]
detect_editing_by_reditools_13m_noSNP_for_ind <- function(bam_cell_anno_= bam_cell_anno, 
                                                          bam_ind_dir_ = "", 
                                                          RE_out_dir_ = "", 
                                                          REDIportal_tab_ = "ref/REDIportal_noRSid.tab.gz",
                                                          genome_fa_ = "genome.fa",
                                                          min_cov_=3,
                                                          read_strand_=1,
                                                          min_edited_reads=2,
                                                          min_base_quality_=30,
                                                          min_mapping_quality = 30,
                                                          parThreads_=3, 
                                                          subThreads_=20, 
                                                          python2.7_exe_="/usr/bin/python2.7" ){
  start_time_i <- Sys.time()
  if(!dir.exists(RE_out_dir_)){
    dir.create(RE_out_dir_, showWarnings = TRUE, recursive = TRUE)
  }
  if(!dir.exists(RE_out_dir_)){
    if(!dir.exists(RE_out_dir_)){ stop("==> |",RE_out_dir_," | can not be created")}
  }
  
  # 
  bam_ind <- unique(bam_cell_anno_[,.(donor)])
  bam_ind$bam_path <- paste0(bam_ind_dir_,"/",bam_ind$donor,".bam")
  bam_ind$RE_out_path <- paste0(RE_out_dir_, "/", bam_ind$donor)
  bam_ind$out_exist <- unlist(lapply(bam_ind$RE_out_path, function(x){dir.exists(x)}))
  message("==> in total: ", nrow(bam_ind))
  message("==> left: ", nrow( bam_ind[out_exist==FALSE,]))
  bam_ind <- bam_ind[out_exist==FALSE,]
  Sys.sleep(1)
  # 
  tmp_out <- mclapply(1:nrow(bam_ind), function(i){
    ind_i <- bam_ind[i,]$donor
    bam_i <- bam_ind[i]$bam_path
    if(!file.exists(bam_i)){ warnings("  ===> i: ",i," | bam file doesn't exist!");return(i) }
    bam_i <- normalizePath(bam_i)
    Re_out_i <- bam_ind[i]$RE_out_path
    
    # Add -S, only consider the reads of inferred strand
    system(paste0(python2.7_exe_, " src/REDItoolKnown.py -i ",
                  bam_i," -f ",genome_fa_ ," -l ",REDIportal_tab_, 
                  " -S -r 5 -t ",subThreads_," -c ",min_cov_," -q ",min_base_quality_, " -m ", min_mapping_quality ," -s ",read_strand," -v ",min_edited_reads," -n 0 -o ",Re_out_i))
    
    # refine:
    refined_ed_path_i <- paste0(Re_out_i, "_refine.txt")
    #
    ed_exists <- data.table(RE_out_path = list.files(Re_out_i, pattern = "outTable_", recursive = TRUE, full.names = TRUE))
    if(nrow(ed_exists)==0){return()}
    red_i <- fread(ed_exists$RE_out_path, header = TRUE)
    if(nrow(red_i)==0){return(0)}
    # remove multiple subs at one site (overwrite)
    red_i <- red_i[!str_detect(AllSubs, " ")]
    if(nrow(red_i)==0){return(0)}
    # only keep chrom 1-22,X,Y:
    if(str_detect(red_i[1,]$Region, "chr")){
      red_i <- red_i[Region %in% paste0("chr", c(1:22, "X","Y"))]
    }else{
      red_i <- red_i[Region %in% paste0("", c(1:22, "X","Y"))]
    }
    if(nrow(red_i)==0){return(0)}
    # rename coverage quality name:
    names(red_i)[str_detect(names(red_i),"Coverage-q")] <- "Cov_over_q"
    
    # filter based on dynamic cov cutoff:
    red_i <- red_i[Cov_over_q >= min_cov_]
    
    # reformat:
    red_i_sub_1 <- red_i[AllSubs!="-"]; 
    red_i_sub_2 <- red_i[AllSubs=="-"];
    if(nrow(red_i_sub_1)>0){
      red_i_sub_1$`BaseCount[A,C,G,T]`<- str_remove_all(red_i_sub_1$`BaseCount[A,C,G,T]`, regex("\\[|\\]"))
      red_i_sub_1 <- as.data.table(tidyr::separate(red_i_sub_1, col = "BaseCount[A,C,G,T]", sep = ", ", into = c("A","C","G","T")))
      red_i_sub_1 <- as.data.table(tidyr::separate(red_i_sub_1, col = "AllSubs", sep = "", into = c("nn","Ref", "Alt")));red_i_sub_1$nn <- NULL
      red_i_sub_1 <- red_i_sub_1[,.(chrom=Region, position=Position, Ref, Alt, strand=ifelse(Strand==1,"+","-"), coverage=Cov_over_q, editlevel=Frequency, A, C, G, T )]
      red_i_sub_1 <- melt(red_i_sub_1, id.vars = c("chrom", "position", "Ref", "Alt", "strand", "coverage", "editlevel"), variable.name = "seq_base", value.name = "reads_num")
      red_i_sub_1 <- red_i_sub_1[Alt==seq_base,.(chrom, position, Ref, Alt, strand, coverage, editedreads=as.numeric(reads_num), editlevel)]
    }else{
      red_i_sub_1 <- data.table(chrom=character(0), position=numeric(0), Ref=character(0), Alt=character(0), strand=character(0), coverage=numeric(0), editedreads=numeric(0), editlevel=numeric(0))
    }
    if(nrow(red_i_sub_2)>0){
      red_i_sub_2 <- red_i_sub_2[,.(chrom=Region, position=Position, Ref=Reference, Alt=Reference, strand=ifelse(Strand==1,"+","-"), coverage=Cov_over_q, editedreads=0, editlevel=Frequency)]
    }else{
      red_i_sub_2 <- data.table(chrom=character(0), position=numeric(0), Ref=character(0), Alt=character(0), strand=character(0), coverage=numeric(0), editedreads=numeric(0), editlevel=numeric(0))
    }
    
    #
    red_i <- rbind(red_i_sub_1, red_i_sub_2)[order(chrom, position)]
    red_i$editlevel <- red_i$editedreads / red_i$coverage
    fwrite(red_i, refined_ed_path_i, sep="\t")
    
    return(NULL)
  }, mc.cores = parThreads_)
  end_time_i <- Sys.time()
  running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
  return(running_time)
}


summary_ed_cov_cellNum_for_cell_type_ind_pool <- function(bam_cell_anno_,
                                                          RE_out_reditools_ ="reditools_13m_ind/",
                                                          bam_split_by_celltype_ind_path = "bam_split_by_celltype_ind/",
                                                          ref_fa_ = "genome.fa",
                                                          out_dir_= "reditools_13m_ind_celltype/",
                                                          min_mapping_quality=0,
                                                          bed_splited_lines=1000,
                                                          parThreads_=40
){
  start_time_i <- Sys.time()
  bam_ind <- unique(bam_cell_anno_[,.(donor)])
  bam_ind$RE_out_path <- paste0(RE_out_reditools_, "/", bam_ind$donor,"_refine.txt")
  bam_ind <- bam_ind[file.exists(bam_ind$RE_out_path),]
  inds <- unique(bam_ind[,.(donor)])
  tmp_out <- lapply(1:nrow(inds), function(i){
    ind_i <- inds[i,]$donor
    bed_path_i <- str_replace(bam_ind[donor==ind_i]$RE_out_path, "_refine.txt", ".bed")
    # if(file.exists(bed_path_i)){message("==> ",i," | output exists! next");return()}
    
    # editing for donor:
    re_i <- fread(bam_ind[donor==ind_i]$RE_out_path, header=TRUE, sep="\t")
    re_i <- re_i[editlevel>0][order(chrom, position)]
    if(nrow(re_i)==0){ return() }
    fwrite( re_i[,.(chrom, as.character(as.integer(position-1)), as.character(as.integer(position)), 
                    paste(chrom, position, Ref, Alt, strand, sep="_"), editedreads/coverage, strand )], bed_path_i, sep="\t", col.names = FALSE)
    
    # bam file of cell type for donor:
    # bam_cell_anno_i <- bam_cell_anno_[donor==ind_i,]
    cell_types_i <- unique(bam_cell_anno_[donor==ind_i,][,.(cell_type)])
    cell_types_i$bam_path <- paste0(bam_split_by_celltype_ind_path, "/", cell_types_i$cell_type,"/", ind_i,".bam")
    cell_types_i <- cell_types_i[file.exists(bam_path)]
    if(nrow(cell_types_i)==0){return()}
    
    # extract reads info for each ed site for each cell type:
    tmp_out <- lapply(1:nrow(cell_types_i), function(j){
      message("==> ",i," | ",j)
      cell_type_j <- cell_types_i[j]$cell_type
      bam_path_j <- cell_types_i[j]$bam_path
      cell_info_out_dir_j <- paste0(out_dir_, "/", cell_type_j,"/",ind_i,"/")
      pileup_path_j = paste0(cell_info_out_dir_j,"pileup_all.txt")
      cell_reads_path_j = paste0(cell_info_out_dir_j,"cell_reads.txt")
      mmp_path_j <- paste0(cell_info_out_dir_j,"misMatch_pop_each_cell.txt")
      # if(dir.exists(cell_info_out_dir_j)){ message("==> ",i," | ",j," exists, next");return() }
      
      # extract pileup file
      if(!file.exists(pileup_path_j)){
        tmp_out2 <- extract_cell_infomation_poolID(pos_bed_ =bed_path_i, 
                                            bam_file_ = bam_path_j,
                                            cell_info_out_dir_ = cell_info_out_dir_j,
                                            ref_fa_= ref_fa_,
                                            min_mapping_quality=min_mapping_quality,
                                            bed_splited_lines=bed_splited_lines,
                                            parThreads_ = parThreads_)
      }
      if(!file.exists(pileup_path_j)){ return() }
      
      # extract cell_reads file
      tmp_out3 <- refine_summary_pileup_poolID_Nathan(bam_cell_anno_ = bam_cell_anno_,
                                               pileup_path_ = pileup_path_j,
                                               out_cell_reads_path = cell_reads_path_j,
                                               base_quality_score_cutoff_=0,
                                               type="ed",
                                               parThreads_= parThreads_)
      if(!is.null(tmp_out3)){
        # extract mmp file
        tmp_out4 <- cell_reads_summary(cell_reads_path = cell_reads_path_j,
                                       position_bed = bed_path_i,
                                       out_misMatch_pop_each_cell_ = mmp_path_j,
                                       type="ed",   
                                       parThreads = parThreads_)
      }
      
      message("=================")
      return()
    })
    
  })
  end_time_i <- Sys.time()
  running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
  return(running_time)
}


extract_cell_infomation_poolID <- function(pos_bed_ ="ed_refined.bed",
                                           bam_file_ = "HMN171215.bam",
                                           cell_info_out_dir_ = "/HMN171215/",
                                           ref_fa_ = "genome.fa",
                                           min_mapping_quality=0,
                                           bed_splited_lines=1000,
                                           parThreads_=20) {
  
  if(pos_bed_==""){
    if(file.exists(paste0(cell_info_out_dir_,"/pileup_all.txt"))){ message(" ==> output file exist! next"); message("[",paste0(cell_info_out_dir_,"/pileup_all.txt"),"]");return(NA) }
    if(!dir.exists(cell_info_out_dir_)){dir.create(cell_info_out_dir_, recursive = T)}
    bam_splited_dir <- paste0(cell_info_out_dir_,"/bam_splited/")
    if(!dir.exists(bam_splited_dir)){dir.create(bam_splited_dir, recursive = T)}
    pileup_splited_dir <- paste0(cell_info_out_dir_,"/pileup_splited/")
    if(!dir.exists(pileup_splited_dir)){dir.create(pileup_splited_dir, recursive = T)}
    
    # split chrom of bam file:
    mclapply(paste0("chr", c(1:22, "X", "Y") ), function(c){
      print(c)
      system(paste0( "samtools view -b ",bam_file_, " ",c," >",bam_splited_dir,c,".bam" ))
      return()
    }, mc.cores = parThreads_)
    splited_bams <- list.files(bam_splited_dir, regex(".bam"), full.names = TRUE)
    
    # start samtools mpileup, parallel:
    message("  ==> In total of ",length(splited_bams), " bam files plited by chromosome")
    tmp_out <- mclapply(1:length(splited_bams), function(i){
      message("   ==> running splited files: ",i,"/",length(splited_bams)," | ", date())
      bam_i <- splited_bams[i]
      pileup_out_i <- paste0(pileup_splited_dir, str_remove(basename(bam_i), regex(".bam$")),".txt")
      suppressMessages(system(paste0("samtools mpileup --output-extra CB -B -d 100000 -q ",min_mapping_quality," ", bam_i, " -f ", ref_fa_, " | awk '{if($4>=10 && $5 ~\"[ATGCatgc]\" )print}' >", pileup_out_i)))
      return(NULL)
    }, mc.cores = length(splited_bams))
    # 
    message(" ==> Combine splited pileup files:")
    system(paste0("cd ", pileup_splited_dir,"; cat ./* >../pileup_all.txt"))
    message(" ==> Remove splited pileup files:")
    # system(paste0("rm  ", pileup_splited_dir,"/*"))
    system(paste0("rm  ", bam_splited_dir,"/*"))
    return(NA)
    
  }else{
    if(file.exists(paste0(cell_info_out_dir_, "/pileup_all.txt"))){ 
      message(" ==> output file exist! next"); 
      message("[",paste0(cell_info_out_dir_, "/pileup_all.txt"),"]");return(NA) 
    }
    if(!dir.exists(cell_info_out_dir_)){dir.create(cell_info_out_dir_, recursive = T)}
    bed_splited_dir <- paste0(cell_info_out_dir_,"/bed_splited/")
    if(!dir.exists(bed_splited_dir)){dir.create(bed_splited_dir, recursive = T)}
    pileup_splited_dir <- paste0(cell_info_out_dir_,"/pileup_splited/")
    if(!dir.exists(pileup_splited_dir)){dir.create(pileup_splited_dir, recursive = T)}
    
    # split bed file:
    system(paste0("cd ", bed_splited_dir,"; split -l ",bed_splited_lines," ", normalizePath(pos_bed_)))
    splited_beds <- list.files(bed_splited_dir)
    
    # parallel:
    message("  ==> In total of ",length(splited_beds), "files with splited lines: ", bed_splited_lines)
    a <- mclapply(1:length(splited_beds), function(i){
      # message("   ==> running splited files: ",i,"/",length(splited_beds)," | ", date())
      bed_i <- paste0(bed_splited_dir, splited_beds[i])
      pileup_out_i <- paste0(pileup_splited_dir, splited_beds[i])
      suppressMessages(system(paste0("samtools mpileup --output-extra CB,RG -A -B -d 100000 -q ",min_mapping_quality," ", bam_file_, " -f ", ref_fa_, " -l ", bed_i," -o ", pileup_out_i)))
      return(NULL)
    }, mc.cores = parThreads_)
    # 
    message(" ==> Combine splited pileup files:")
    system(paste0("cd ", pileup_splited_dir,"; cat ./* >../pileup_all.txt"))
    if(file.exists(paste0(pileup_splited_dir,"../pileup_all.txt" ))){
      num_lines <- as.numeric(system(paste0("cd ", pileup_splited_dir,"; wc -l ../pileup_all.txt|awk '{print $1}' "), intern = TRUE))
      if(num_lines ==0){
        file.remove(paste0(pileup_splited_dir,"../pileup_all.txt" ))
      }
    }
    message(" ==> Remove splited pileup files:")
    system(paste0("rm  ", pileup_splited_dir,"/*"))
    system(paste0("rm  ", bed_splited_dir,"/*"))
    system(paste0("rmdir  ", pileup_splited_dir))
    system(paste0("rmdir  ", bed_splited_dir))
    return(NA)
  }
  
}


refine_summary_pileup_poolID_Nathan <- function(bam_cell_anno_,
                                                pileup_path_ = "pileup_all.txt",
                                                out_cell_reads_path = "cell_reads.txt",
                                                type="ed",
                                                splite_input=FALSE,
                                                base_quality_score_cutoff_=30,
                                                parThreads_= 20){
  
  if(type=="ed"){
    message("  ==> Start loading file: ")
    ind_x <- basename(dirname(pileup_path_))
    cluster_x <- basename(dirname(dirname(pileup_path_)))
    bam_cell_anno_x <- bam_cell_anno_[cell_type==cluster_x & donor == ind_x,.(cell_type, donor, pool, barcode)]
    setindex(bam_cell_anno_x, barcode)
    
    pileup <- fread(pileup_path_, sep="\t", header=FALSE)
    if(nrow(pileup)==0){return()}
    names(pileup) <- c("chrom", "position", "ref", "readCount", "bases","baseQuality", "barcode", "ReadGroup")
    # remove "^*" and "$", note , first remove ^*, then remove $, due to ^$
    message("  ==> Start remove read start and end sign, total of [", nrow(pileup[str_detect(bases, regex("[$\\^]"))]), "] records")
    pileup[str_detect(bases, regex("[$\\^]")),"bases"] <- unlist(mclapply(pileup[str_detect(bases, regex("[$\\^]"))]$bases , function(x){ str_remove_all(str_remove_all(x,regex("[\\^].{1}")), fixed("$")) }, mc.cores = round(parThreads_/2,0) ))
    # remove indels:
    message("  ==> Start remove read indels, total of [", nrow( pileup[str_detect(bases, regex('[+-]')),] ), "] records")
    
    indel_bases <- pileup[str_detect(bases, regex('[+-]'))]$bases
    if(length(indel_bases)>0){
      pileup[str_detect(bases, regex('[+-]')),"bases"] <- unlist(mclapply(1:length(indel_bases), function(x){ 
        x <- indel_bases[x]
        x_char = str_split_1(x,""); 
        num_pos = as.numeric(str_match_all(x, regex("[0-9]{1,10}"))[[1]])
        if(length(num_pos)==0){return(x)}
        xx=x
        for(n in unique(sort(num_pos, decreasing = TRUE)) ){
          n_pos <- as.data.table(str_locate_all(xx, as.character(n))[[1]])
          # 
          for( nn in 1:nrow(n_pos) ){
            x_char[(n_pos[nn]$start-1): (n_pos[nn]$end+n)] <- "@"
            xx=paste0(x_char, collapse = "")
          }
        }
        xx = str_remove_all(xx, "@")
        return(xx)
      }, mc.cores = round(parThreads_/2,0) ))
    }
    #
    # if(length(which(nchar(pileup$bases) != nchar(pileup$baseQuality)))>0){
    #   warning("  ==> Not all positions are cleaned, ", length(which(nchar(pileup$bases) != nchar(pileup$baseQuality)))," records ,has been removed.")
    #   pileup <- pileup[ which(nchar(pileup$bases) == nchar(pileup$baseQuality)), ]
    #   return()
    # }else{
    #   message("  ==> Pileup file has been cleaned!")
    # }
    # message(" ==> removed ", nrow(pileup[bases==""]), "/", nrow(pileup))
    pileup <- pileup[bases!=""]
    if(nrow(pileup)==0){return()}
    # The remaining symbol types are as follows："<" ">" "c" "," "a" "g" "t" "." "G" "C" "A" "T" "*"
    
    message("  ==> Start remove pos with skip(<>*), which due to junction:")
    # dim(pileup[str_detect(bases, "[ATCGatcg]")])
    # dim(pileup[str_detect(bases, "[Nn]")])
    no_junction_pos <- str_locate_all(pileup$bases, "[*><]")
    junction_pos <- lapply(no_junction_pos, invert_match)
    #
    pileup$bases <- unlist(lapply(str_sub_all(pileup$bases, junction_pos ), function(x){ paste0(x, collapse = "") }))
    pileup$baseQuality <- unlist(lapply(str_sub_all(pileup$baseQuality, junction_pos), function(x){ paste0(x, collapse = "") }))
    pileup$barcode <- unlist(mclapply(1:nrow(pileup), function(i){
      cell_id_splited <- str_split_1(pileup[i,]$barcode, ",")
      return( paste0(cell_id_splited[ setdiff(1:length(cell_id_splited), no_junction_pos[[i]][,1] ) ] , collapse=",") )
    }, mc.cores = parThreads_))
    # pileup$ReadGroup <- unlist(mclapply(1:nrow(pileup), function(i){
    #   rd_i_splited <- unlist(lapply(str_split(pileup[i,]$ReadGroup,",")[[1]], function(x){ str_split(x,":")[[1]][1] }))
    #   return( paste0(rd_i_splited[ setdiff(1:length(rd_i_splited), no_junction_pos[[i]][,1] ) ] , collapse=",") )
    # }, mc.cores = parThreads_))
    
    # message(" ==> removed ", nrow(pileup[bases==""]), "/", nrow(pileup))
    pileup <- pileup[bases!=""]
    if(nrow(pileup)==0){return()}
    
    message(" ==> filter reads with low quality with cutoff: ",base_quality_score_cutoff_ )
    pileup <- rbindlist(mclapply(1:nrow(pileup), function(i){
      quality_score_i <- DescTools::CharToAsc(str_split_1(pileup[i]$baseQuality,""))
      if(all(quality_score_i>=base_quality_score_cutoff_)){
        return(cbind(pileup[i], data.table(removed_reads=0)))
      }else if(length(which(quality_score_i>=base_quality_score_cutoff_))==0){
        return(NULL)
      }else{
        bases_good_i <- which(quality_score_i>=base_quality_score_cutoff_)
        return(pileup[i,][,.(chrom, position, ref, readCount=length(bases_good_i), bases=paste0(str_split_1(bases, "")[bases_good_i],collapse = ""), baseQuality = paste0(str_split_1(pileup[i]$baseQuality,"")[bases_good_i], collapse = ""), barcode = paste0(str_split_1(barcode, ",")[bases_good_i], collapse = ""), removed_reads=length(quality_score_i)-length(which(quality_score_i>=base_quality_score_cutoff_)) )])
      }
    }, mc.cores = parThreads_))
    # message(" ==> removed ", nrow(pileup[bases==""]), "/", nrow(pileup))
    pileup <- pileup[bases!=""]
    if(nrow(pileup)==0){return()}
    message("  ==> Start summary: ")
    # If there are two strands of mutations at each site, two rows of results are retained and counted.
    # foward reads by cell
    
    # procell ReadGroup
    # combine barcode with ReadGroup:
    pileup$rg_cell_id <- unlist(mclapply(1:nrow(pileup), function(i){
      cell_id_i <- unlist(str_split(pileup[i,]$barcode,",")[[1]])
      cell_id_splited_dt <- bam_cell_anno_x[cell_id_i, on="barcode"]
      cell_id_splited_dt[is.na(pool), "pool"] <- na.omit(cell_id_splited_dt[1,])$pool
      cell_id_splited_dt <- cell_id_splited_dt[1:length(cell_id_i),]
      return(paste0(paste0(cell_id_splited_dt$pool, "#", cell_id_splited_dt$barcode), collapse = ","))
    }, mc.cores = parThreads_))
    #
    cell_reads <- rbindlist(mclapply(1:nrow(pileup), function(i){
      pileup_i <- pileup[i,]
      cell_reads_i <- as.data.table(table(pileup_i[,.(bases=str_split_1(pileup_i$bases, ""), barcode = str_split_1(pileup_i$rg_cell_id, ","))]))[N>0]
      names(cell_reads_i) <- c("base_pileup", "barcode", "num_reads")
      # message(" ==> i: ",i," | ",ncol(cell_reads_i))
      return(cbind(pileup_i[,.(chrom, position, ref)], cell_reads_i))
    }, mc.cores = parThreads_))
    message("  ==> Finished.. total of ",nrow(cell_reads), " records.. writing... to..")
    message("  ==> [",out_cell_reads_path,"]")
    
    # forward strand: "." "G" "C" "A" "T"
    # if(str_detect(pileup_i$bases, regex("[ATCG.]")))
    fwrite(cell_reads, out_cell_reads_path, sep="\t")
    if(nrow(cell_reads) > 0){
      return("Done")
    }else{
      return()
    }
  }
}

cell_reads_summary <- function(cell_reads_path = "cell_reads.txt",
                               position_bed ="ed_refined.bed",
                               out_misMatch_pop_each_cell_ = "misMatch_pop_each_cell.txt",
                               cov_reads_cutoff = 3,
                               mismatch_reads_cutoff = 1,
                               min_cov_cells_per_site = 3,
                               type="ed",   
                               parThreads = 30
){
  if(type=="deNovo"){
    misMatch_reads_per_site_strand_cutoff <- 3
    cov_reads_per_site_cutoff <- 10
    #
    cell_reads <- fread(cell_reads_path, header = TRUE, sep="\t")
    setindex(cell_reads, chrom, position)
    #
    # determine Alt and strand:
    cell_reads[,c("strand"):=.( ifelse(base_pileup %in% c(",","a","t", "c","g"), "-","+") )]
    cell_reads[base_pileup %in% c(",", "."), "base_pileup"] <- cell_reads[base_pileup %in% c(",", "."), ]$ref
    cell_reads$base_pileup <- str_to_upper(cell_reads$base_pileup)
    # summary mismatch reads:
    cell_reads_pheno_summary <- cell_reads[ref != base_pileup][,.(cell_reads=sum(num_reads)),by= c("chrom", "position","ref","base_pileup", "strand")]
    # filter out site with mismatch reads <=2:
    cell_reads_pheno_summary <- cell_reads_pheno_summary[cell_reads>= misMatch_reads_per_site_strand_cutoff]
    cell_reads_pheno_summary_dup <- merge(cell_reads_pheno_summary, unique(cell_reads_pheno_summary[which(duplicated(cell_reads_pheno_summary[,1:3])),1:3]), by=c("chrom", "position", "ref"), sort=FALSE)
    cell_reads_pheno_summary_dedup <- fsetdiff(cell_reads_pheno_summary, cell_reads_pheno_summary_dup)
    # only retation the base with the largest num at the same position and same strand:
    cell_reads_pheno_summary_dup <- cell_reads_pheno_summary_dup[order(-cell_reads)][,.SD[1,],by=c("chrom", "position","strand")]
    # retaion both strand:
    # cell_reads_pheno_summary_dup2 <- merge(cell_reads_pheno_summary_dup, unique(cell_reads_pheno_summary_dup[which(duplicated(cell_reads_pheno_summary_dup[,.(chrom, position, ref)])),.(chrom, position, ref)]),by=c("chrom", "position", "ref"))
    cell_reads_pheno_summary <- rbind(cell_reads_pheno_summary_dup, cell_reads_pheno_summary_dedup); rm(cell_reads_pheno_summary_dedup,cell_reads_pheno_summary_dup)
    cell_reads_pheno_summary[,c("phenotype_id"):=.(paste0(chrom,"_",position,"_",ref,"_", base_pileup,"_",strand))]
    cell_reads <- merge(cell_reads, cell_reads_pheno_summary[,.(chrom, position, ref, strand, phenotype_id)], by=c("chrom", "position","ref","strand"))
    
    # only keeep cov reads >= 10
    cell_reads <- merge(cell_reads, cell_reads[,.(cov_reads=sum(num_reads)),by= c("chrom", "position","ref")][cov_reads >=cov_reads_per_site_cutoff][,.(chrom, position, ref)], by= c("chrom", "position", "ref"))
    #
    phenotypes <- unique(cell_reads[,.(phenotype_id, chrom, position, strand)])
    phenotypes$batch <- rep(1:parThreads,each =ceiling(nrow(phenotypes)/parThreads))[1:nrow(phenotypes)]
    misMatch_pop_each_cell <- rbindlist(mclapply(1:parThreads, function(i){
      cell_reads[phenotypes[batch==i,]$phenotype_id, on="phenotype_id"][,.(misMatch_reads=sum(.SD[base_pileup!=ref]$num_reads), cov_reads=sum(.SD$num_reads) ),by=c("barcode", "phenotype_id")]
    }, mc.cores = parThreads))
    misMatch_pop_each_cell[,c("misMatch_prop"):=.(misMatch_reads/cov_reads )]
    misMatch_pop_each_cell <- merge(misMatch_pop_each_cell, phenotypes[,-c("batch")], by= "phenotype_id", sort=FALSE )[,.(phenotype_id, chrom, position,strand, barcode, misMatch_reads,cov_reads, misMatch_prop)]
    
    # 
    # misMatch_pop_each_cell <- rbindlist(mclapply(1:nrow(phenotypes), function(i){
    #   cell_reads_i <- na.omit(cell_reads[phenotypes[i,]$phenotype_id, on=c("phenotype_id")])
    #   misMatch_pop_each_cell_i <- cell_reads_i[,.(misMatch_reads=sum(.SD[base_pileup!=ref]$num_reads), cov_reads=sum(.SD$num_reads) ),by="barcode"]
    #   misMatch_pop_each_cell_i[,c("misMatch_prop"):=.(misMatch_reads/cov_reads )]
    #   return(cbind( unique(cell_reads_i[,.(chrom, position, strand, phenotype_id)]), misMatch_pop_each_cell_i)[,.(phenotype_id, chrom, position,strand, barcode, misMatch_reads,cov_reads, misMatch_prop)] )
    # }, mc.cores = parThreads))
    fwrite(misMatch_pop_each_cell, out_misMatch_pop_each_cell_, sep="\t", col.names = TRUE )
    # misMatch_cell_prop <- misMatch_pop_each_cell[,.(cell_prop =nrow(.SD[misMatch_prop>0.5]) / nrow(.SD)),by=c("chrom", "position", "strand", "phenotype_id")]
    return()
  }
  if(type=="ed"){
    # for ed:
    ed <- fread(position_bed, header=FALSE)
    names(ed) <- c("chrom", "start","position", "phenotype_id", "coverage", "strand")
    cell_reads <- fread(cell_reads_path, header = TRUE, sep="\t")
    setindex(cell_reads, chrom, position)
    ed <- merge(ed, unique(cell_reads[,.(chrom, position)]), by=c("chrom", "position"),sort=FALSE )
    misMatch_pop_each_cell <- rbindlist(mclapply(1:nrow(ed), function(i){
      ed_i <- ed[i,]
      cell_reads_i <- na.omit(cell_reads[ed_i, on=c("chrom", "position")])
      if(nrow(cell_reads_i)==0){return(NULL)}
      if( ed_i$strand=="-" ){
        cell_reads_i <- cell_reads_i[ base_pileup %in% c(",","a","t", "c","g"), ]
        if(nrow(cell_reads_i)==0){return(NULL)}
        # mismatch proportion in each cell:
        misMatch_pop_each_cell_i <- cell_reads_i[,.(misMatch_reads=sum(.SD[base_pileup!=","]$num_reads), cov_reads=sum(.SD$num_reads) ),by="barcode"][order(misMatch_reads)]
        misMatch_pop_each_cell_i
      }else{
        cell_reads_i <- cell_reads_i[ base_pileup %in% c(".","A","T", "C","G"), ]
        if(nrow(cell_reads_i)==0){return(NULL)}
        # mismatch proportion in each cell:
        misMatch_pop_each_cell_i <- cell_reads_i[,.(misMatch_reads=sum(.SD[base_pileup!="."]$num_reads), cov_reads=sum(.SD$num_reads) ),by="barcode"][order(misMatch_reads)]
        misMatch_pop_each_cell_i
      }
      misMatch_pop_each_cell_i$misMatch_prop <- misMatch_pop_each_cell_i$misMatch_reads / misMatch_pop_each_cell_i$cov_reads
      return(cbind(cell_reads_i[,.(chrom, position, strand, phenotype_id)], misMatch_pop_each_cell_i))
    }, mc.cores = parThreads))
    fwrite(misMatch_pop_each_cell, out_misMatch_pop_each_cell_, sep="\t", col.names = TRUE )
    # misMatch_cell_prop <- misMatch_pop_each_cell[,.(cell_prop =nrow(.SD[misMatch_prop>0.5]) / nrow(.SD)),by=c("chrom", "position", "strand", "phenotype_id")]
    # p <- boxplot(misMatch_cell_prop$cell_prop)
    # print(p)
    # quantile(misMatch_cell_prop$cell_prop, seq(0,1,0.01))
    # message(" ==> ",round(nrow(misMatch_cell_prop[cell_prop<0.5])/nrow(misMatch_cell_prop)*100,2), "% editing sites with < 50% mismatched cells")
    return()
  }else{
    # for variant:
    gt <- fread(position_bed, header=FALSE)
    names(gt) <- c("chrom", "start","position", "phenotype_id", "gt", "NA")
    cell_reads <- fread(cell_reads_path)
    # ed <- fread("/media/bora_A/dingruofan/2023-10-16_RNA_editing_project/outputs/Randolph2021Science/reditools_13m_ind/HMN171215.bed", header=FALSE)
    # names(ed) <- c("chrom", "start","position", "phenotype_id", "coverage", "strand")
    gt <- gt[position %in% unique(cell_reads$position),-c("NA")]
    gt <- cbind(gt, as.data.table(tidyr::separate(gt[,.(phenotype_id)], col="phenotype_id", into=c("chrom", "pos", "gt_Ref","gt_Alt"), remove =TRUE, convert = TRUE))[,.(gt_Ref, gt_Alt)])
    gt_cell_reads <- cell_reads[gt, on=c("chrom", "position")]
    # correct gt_ref/gt_alt due to plink error:
    need_correct = FALSE
    if(need_correct){
      gt_cell_reads_corrected <- cbind(gt_cell_reads[ref!=gt_Ref,-c("gt_Ref", "gt_Alt")], gt_cell_reads[ref!=gt_Ref][,.(gt_Ref=gt_Alt, gt_Alt=gt_Ref)])
      gt_cell_reads_corrected$gt_old <- gt_cell_reads_corrected$gt; gt_cell_reads_corrected[gt_old=="1/1", "gt"]<-"0/0";gt_cell_reads_corrected[gt_old=="0/0", "gt"]<-"1/1"
      gt_cell_reads_correct <- gt_cell_reads[ref==gt_Ref,]
      gt_cell_reads <- rbind(gt_cell_reads_correct, gt_cell_reads_corrected[,-c("gt_old")])
      gt_cell_reads <- gt_cell_reads[gt!="0/0"]
    }
    # only keep ref=gt_ref:
    gt_cell_reads <- gt_cell_reads[ref==gt_Ref]
    # correct the alt base in negtive strand based on gt_alt
    # gt_cell_reads <- rbind(gt_cell_reads[!(base_pileup %in% c("a", "t","g","c"))], 
    #                        merge(data.table(base_pileup_old=c("a", "t","g", "c"),base_pileup=c("T","A","C","G")), gt_cell_reads,  by.x="base_pileup_old",by.y="base_pileup", sort=FALSE)[,-c("base_pileup_old")])
    # 
    # only keep pileup Alt == gt alt:()
    gt_cell_reads <- rbind(gt_cell_reads[base_pileup %in% c(".", ",")], gt_cell_reads[!base_pileup %in% c(".", ",") & str_to_upper(base_pileup) == gt_Alt])
    
    # Retention site cov reads num>=3 from all cells
    gt_cell_reads <- merge(gt_cell_reads, gt_cell_reads[,.(cell_reads=sum(num_reads)),by=c("barcode", "phenotype_id")][cell_reads >= cov_reads_cutoff,-c("cell_reads")], by=c("barcode", "phenotype_id"), sort=FALSE)
    # There is at least one mismatch in all cells at a locus so that it can be shown to be a mutation
    gt_cell_reads <- merge(gt_cell_reads,  gt_cell_reads[base_pileup!="." & base_pileup!="," ][,.(cell_reads=sum(num_reads)),by= c("barcode", "phenotype_id")][cell_reads >= mismatch_reads_cutoff], by=c("barcode", "phenotype_id"), sort=FALSE )
    # At least 1 cell at a locus
    gt_vars <- gt_cell_reads[,.(num_cells = length(unique(barcode)) ), by="phenotype_id"][num_cells >= min_cov_cells_per_site]
    gt_cell_reads <- gt_cell_reads[phenotype_id %in% gt_vars$phenotype_id]
    
    phenotypes <- unique(gt_cell_reads[,.(phenotype_id, chrom, position, gt)])
    phenotypes$batch <- rep(1:parThreads,each =ceiling(nrow(phenotypes)/parThreads))[1:nrow(phenotypes)]
    mutated_pop_each_cell <- rbindlist(mclapply(1:parThreads, function(i){
      gt_cell_reads[phenotypes[batch==i,]$phenotype_id, on="phenotype_id"][,.(misMatch_reads=sum(.SD[base_pileup!=ref]$num_reads), cov_reads=sum(.SD$num_reads) ),by=c("barcode", "phenotype_id")]
    }, mc.cores = parThreads))
    mutated_pop_each_cell[,c("misMatch_prop"):=.(misMatch_reads/cov_reads )]
    mutated_pop_each_cell <- merge(mutated_pop_each_cell, phenotypes[,-c("batch")], by= "phenotype_id", sort=FALSE )[,.(phenotype_id, chrom, position, barcode, misMatch_reads,cov_reads, misMatch_prop, gt)]
    
    
    # mutated_pop_each_cell <- rbindlist(mclapply(1:nrow(gt), function(i){
    #   gt_vars_i <- gt_vars[i,]
    #   gt_cell_reads_i <- na.omit(gt_cell_reads[gt_vars_i, on=c("phenotype_id")])
    #   if(nrow(gt_cell_reads_i)==0){return(NULL)}
    #   # Statistics the proportion of cells with at least one mutated reads
    #   misMatch_pop_each_cell_i <- gt_cell_reads_i[,.(misMatch_reads=sum(.SD[base_pileup !="." & base_pileup !=","]$num_reads), cov_reads=sum(.SD$num_reads) ),by="barcode"][order(misMatch_reads)]
    #   misMatch_pop_each_cell_i$misMatch_prop <- misMatch_pop_each_cell_i$misMatch_reads / misMatch_pop_each_cell_i$cov_reads
    #   return( cbind(unique(gt_cell_reads_i[,.(chrom, position, phenotype_id)]), misMatch_pop_each_cell_i) )
    # }, mc.cores = parThreads))
    
    fwrite(mutated_pop_each_cell, out_misMatch_pop_each_cell_, sep="\t", col.names = TRUE )
    # quantile(mutated_pop_each_cell$misMatch_prop, seq(0,1,0.01))
    message(" ==> ",round(nrow(mutated_pop_each_cell[misMatch_prop>0.5])/nrow(mutated_pop_each_cell)*100,2), "% SNPs with > 50% mismatched cells")
    return(mutated_pop_each_cell)
  }
}


# filter edsites for each cell type each donor (>=3 cells per site)
filter_ed_by_cell_and_cov_ind <- function(bam_cell_anno_=bam_cell_anno,
                                          RE_out_reditools_ ="data/ed/reditools_13m_ind/",
                                          mmp_out_dir= "data/ed/reditools_13m_ind_celltype/",
                                          mmp_file_name="misMatch_pop_each_cell.txt",
                                          min_num_covered_cells = 3,
                                          min_cov_reads = 3,
                                          min_num_covered_cells_0 =2,
                                          min_cov_reads_0 = 3,
                                          parThreads=20
){
  start_time_i <- Sys.time()
  bam_ind <- unique(bam_cell_anno_[,.(donor)])
  bam_ind$RE_out_path <- paste0(RE_out_reditools_, "/", bam_ind$donor,"_refine.txt")
  bam_ind <- bam_ind[file.exists(bam_ind$RE_out_path),]
  inds <- unique(bam_ind[,.(donor)])
  clusters <- unique(bam_cell_anno_[,.(cell_type)])
  lapply(1:nrow(clusters), function(i){
    cluster_i <- clusters[i]$cell_type
    message("==> ",i," | ",cluster_i, " | mmp file: ",mmp_file_name," | ",date())
    a <- mclapply(1:nrow(inds), function(j){
      ind_j <- inds[j]$donor
      mmp_path_j <- paste0(mmp_out_dir, "/", cluster_i,"/",ind_j,"/",mmp_file_name)
      ed_filter_cell_out_path_j <- paste0(mmp_out_dir, "/", cluster_i,"/",ind_j,"/ed_filtered_cells.txt")
      if(!file.exists(mmp_path_j)){ message("mmp does not exist!");return() }
      if(file.exists(ed_filter_cell_out_path_j)){ file.remove(ed_filter_cell_out_path_j) }
      mmp_j <- fread(mmp_path_j, header = TRUE, sep="\t")
      phenotypes_j <- unique(mmp_j[,.(chrom, position, strand,phenotype_id)])
      phenotypes_j$ed_type <- unlist(mclapply(phenotypes_j$phenotype_id, function(x){ x=str_split_fixed(x, "_", 5); return(paste0(x[3],x[4])) }, mc.cores = parThreads))
      
      # if ed level>0, ed of covered cells >= 3:
      ed_j_1 <- mmp_j[phenotype_id %in% mmp_j[phenotype_id %in% mmp_j[misMatch_prop>0]$phenotype_id,][,.(cell_num=length(barcode)), by="phenotype_id"][cell_num>= min_num_covered_cells, ]$phenotype_id ]
      # calculate editing level:
      ed_j_1 <- ed_j_1[,.(cell_num=length(barcode),cov_reads=sum(cov_reads), 
                          edited_reads = sum(misMatch_reads),
                          ed_level = sum(misMatch_reads) / sum(cov_reads)), 
                       by="phenotype_id"]
      ed_j_1 <- ed_j_1[cov_reads >= min_cov_reads, ]
      ed_j_1 <- merge(phenotypes_j,ed_j_1, by="phenotype_id")
      
      ####### characterization and comparison --- start
      # # editing level negtivel correlated with cell num:
      # ed_j_1[,.(cor(cell_num, ed_level)), by="ed_type"]
      # # editing level are different among ed types:
      # ggplot(ed_j_1)+
      #   geom_boxplot(aes(x=ed_type, y=ed_level))
      # # compare with last result derived from psuedobulk of cell type + donor:
      # ed_raw <- fread(paste0("/media/bora_A/dingruofan/2023-10-16_RNA_editing_project/outputs/Randolph2021Science/reditools_13m/",cluster_i,"#",ind_j,"_refine_filterGT.txt"))
      # ed_raw$editlevel <- ed_raw$editedreads / ed_raw$coverage
      # ed_raw <- ed_raw[editlevel>0]
      # ed_raw$phenotype_id <- paste(ed_raw$chrom,ed_raw$position, ed_raw$Ref, ed_raw$Alt,ed_raw$strand, sep="_")
      # ed_merged <- merge(ed_raw[,.(phenotype_id, editlevel, coverage)], ed_j_1[,.(phenotype_id, ed_level, cov_reads,edited_reads, ed_type)], by="phenotype_id")
      # ed_raw[phenotype_id %in% setdiff(ed_raw$phenotype_id, ed_j_1$phenotype_id)]
      # ed_j_1[phenotype_id %in% setdiff(ed_j_1$phenotype_id, ed_raw$phenotype_id)][edited_reads>2 & str_detect(phenotype_id,"chr1_")]
      ####### characterization and comparison --- end
      
      # if ed level ==0, just keep
      # if ed level>0, ed of covered cells >= 2:
      ed_j_0 <- mmp_j[phenotype_id %in% mmp_j[phenotype_id %in% mmp_j[,.( ed_level=sum(misMatch_reads) / sum(cov_reads) ), by="phenotype_id"][ed_level==0]$phenotype_id,][,.(cell_num=length(barcode)), by="phenotype_id"][cell_num>= min_num_covered_cells_0, ]$phenotype_id ]
      # calculate editing level:
      ed_j_0 <- ed_j_0[,.(cell_num=length(barcode),cov_reads=sum(cov_reads), edited_reads = sum(misMatch_reads),ed_level = sum(misMatch_reads) / sum(cov_reads)), by="phenotype_id"]
      ed_j_0 <- ed_j_0[cov_reads >= min_cov_reads, ]
      ed_j_0 <- merge(phenotypes_j,ed_j_0, by="phenotype_id")
      
      #
      ed_filtered_cells <- rbind(ed_j_1, ed_j_0)
      fwrite(ed_filtered_cells[,.(phenotype_id, chrom, position=as.character(as.integer(position)), 
                                  strand, ed_type, cell_num, cov_reads, edited_reads, ed_level )], 
             ed_filter_cell_out_path_j, sep="\t", col.names = TRUE)
      return()
    }, mc.cores = parThreads)
  })
  end_time_i <- Sys.time()
  running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
  return(running_time)
}

# raw name: merge_eds_ind_filtered
merge_eds_ind_filtered <- function(bam_cell_anno_= bam_cell_anno, 
                                    mmp_out_dir= "reditools_13m_ind_celltype/",
                                    RE_dcast_out_dir_ = "reditools_13m_ind_merged/",
                                    prop_of_ed_shared_inds = 0.3,
                                    filter_ed_gt_0=FALSE,
                                    parThreads_merge_ = 5 ){
      start_time_i <- Sys.time()
      if(!dir.exists(RE_dcast_out_dir_)){dir.create(RE_dcast_out_dir_)}
      if(!dir.exists(RE_dcast_out_dir_)){ stop("==> |",RE_dcast_out_dir_," | can not be created")}
      
      bam_cluster <- unique(bam_cell_anno_[,.(donor, cell_type)])
      bam_cluster$ed_path <- paste0(mmp_out_dir, "/", bam_cluster$cell_type,"/",bam_cluster$donor,"/ed_filtered_cells.txt")
      bam_cluster <- bam_cluster[file.exists(bam_cluster$ed_path)]
      #
      inds <- unique(bam_cluster[,.(donor)])
      clusters <- unique(bam_cell_anno_[,.(cell_type)])
      #
      a <- lapply(1:nrow(clusters), function(i){
        cluster_i <- clusters[i]$cell_type
        out_bed_i <- paste0(RE_dcast_out_dir_, "/", cluster_i,"_ed.txt")
        out_bed_AG_i <- paste0(RE_dcast_out_dir_, "/", cluster_i,"_ed_AG.txt")
        
        if(file.exists(out_bed_i)){ file.remove(out_bed_i)}
        if(file.exists(out_bed_AG_i)){ file.remove(out_bed_AG_i) }
        
        message("==> ",i,"/",nrow(clusters)," | ",cluster_i, " | ",date())
        bam_cluster_i <- bam_cluster[cell_type== cluster_i]
        # Note that each cell_type person here may be different! So it has to be calculated separately!
        min_num_of_ind_i <- round(prop_of_ed_shared_inds *nrow(bam_cluster_i))
        eds_i <- rbindlist(mclapply(1:nrow(bam_cluster_i), function(j){
          ed_j <- fread(bam_cluster_i[j]$ed_path); ed_j$donor <- bam_cluster_i[j]$donor
          # if(only_AG){ed_j <- ed_j[ed_type=="AG"]}
          ed_j <- ed_j[chrom %in% paste0("chr", c(1:22, "X", "Y"))]
          return(ed_j)
        }, mc.cores = parThreads_merge_))
        if(nrow(eds_i)==0){ message(" ==> no data!"); return()}
        # filter ed by donor num:
        if(filter_ed_gt_0){
          eds_i <- eds_i[phenotype_id %in% eds_i[ed_level>0,][,.(num_ind= length(donor)), by="phenotype_id"][num_ind>=min_num_of_ind_i]$phenotype_id,]
        }else{
          eds_i <- eds_i[phenotype_id %in% eds_i[,.(num_ind= length(donor)), by="phenotype_id"][num_ind>=min_num_of_ind_i]$phenotype_id,]
        }
        # filter ed by donor num( ed level>0):
        if(nrow(eds_i)==0){ message(" ==> after filter, no data!");return()}
        
        # dcast by chrom for numerous records:
        chroms_i <- sort(unique(eds_i$chrom))
        # construct bed output format:
        a <- data.table( matrix(inds$donor, nrow=1) )
        names(a) <- unique( inds$donor )
        a <- a[0,]
        eds_dcast_i <- cbind(eds_i[0,.(`#chr`=chrom, start=position, end=position+1, phenotype_id)], a)
        # create the initiaed output file:
        fwrite(eds_dcast_i, out_bed_i, sep="\t", col.names = TRUE)
        fwrite(eds_dcast_i, out_bed_AG_i, sep="\t", col.names = TRUE)
        rm(a)
        # dcast by chrom for numerous records:
        a <- mclapply(1:length(chroms_i), function(j){
          chrom_j <- chroms_i[j]
          eds_j <- eds_i[chrom %in% chrom_j]
          eds_j$end=eds_j$position+1
          eds_dcast_j <- dcast(eds_j, chrom+position+end+phenotype_id~donor, value.var = "ed_level", fun.aggregate = function(x) if(length(x) == 0) NA_real_ else sum(x, na.rm = TRUE) )
          names(eds_dcast_j)[1:4] <- c("#chr", "start", "end", "phenotype_id")
          eds_dcast_j <- rbind(eds_dcast_i, eds_dcast_j, fill=TRUE)
          # remove editing site all is 0 or NA:
          # remove editing site all is 0 or NA:
          a <- copy(eds_dcast_j[,-1:-4])
          a[is.na(a)]<-0
          a <- as.data.table(lapply(a, as.numeric))
          a$rowSum_value <- rowSums(a)
          eds_dcast_j <- eds_dcast_j[a$rowSum_value>0,]
          eds_dcast_AG_j <- eds_dcast_j[str_detect(phenotype_id, "_A_G_")]
          if(nrow(eds_dcast_j)>0){
            fwrite(eds_dcast_j[order(`#chr`, start)], out_bed_i, sep="\t", append=TRUE, col.names=FALSE)
          }
          if(nrow(eds_dcast_AG_j)>0){
            fwrite(eds_dcast_AG_j[order(`#chr`, start)], out_bed_AG_i, sep="\t", append=TRUE, col.names=FALSE)
          }
          message("==> ",j," | ",chrom_j," | ", nrow(eds_dcast_j)," | AG: ",nrow(eds_dcast_AG_j))
        }, mc.cores = parThreads_merge_)
        return()
      })
      end_time_i <- Sys.time()
      running_time <- as.numeric(difftime(end_time_i, start_time_i,units = "mins"))
      return(running_time)
}


run_editing_analysis <- function(bam_data_=bam_data, args_=args, log_) {
  #' Run complete RNA editing analysis pipeline
  #'
  #' @param bam_data_ BAM files and Cell type annotation data
  #' @param args_ parameters
  #' @param log_ log file
  
  bam_anno <- bam_data_$bam_anno
  bam_cell_anno <- bam_data_$bam_cell_anno
  out_dir <- args$outdir
  parThreads <- args$threads
  py2_path <- args$py2_path
  genome_fa <- args$genome_fa
  read_strand <- args$read_strand
  REDIportal_tab <- "ref/REDIportal_noRSid.tab.gz"
  
  # Create output directories if they don't exist
  dir.create(file.path(out_dir, "bam_split_by_celltype_ind_split"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "bam_split_by_ind"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "bam_split_by_celltype_ind"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "reditools_13m_ind"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "reditools_13m_ind_celltype"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir, "reditools_13m_ind_merged"), showWarnings = FALSE, recursive = TRUE)
  
  # Split BAM files by cell type and individual
  running_time <- split_dedup_of_pool_bam(
    bam_data_ = bam_data_,
    out_dir_ = file.path(out_dir, "bam_split_by_celltype_ind_split"),
    parThreads_ = ifelse(round(parThreads/10) >= 1, round(parThreads/10), 1),
    subThreads_ = ifelse(parThreads >= 10, 10, parThreads)
  )
  message("remove_dups_pools2 running time: ", running_time)
  fwrite(data.table(type = "split_dedup_of_pool_bam", time = running_time), 
         file.path(out_dir, "log.txt"), sep = "\t", append = TRUE)
  
  # Merge BAM files by individual
  running_time <- combine_bam_ind(
    bam_cell_anno_ = bam_data_$bam_cell_anno,
    bam_splited_dir_ = file.path(out_dir, "bam_split_by_celltype_ind_split"),
    out_dir_ = file.path(out_dir, "bam_split_by_ind"),
    parThreads_ = parThreads,
    subThreads_ = 4
  )
  fwrite(data.table(type = "combine_bam_ind", time = running_time), 
         file.path(out_dir, "log.txt"), sep = "\t", append = TRUE)
  message("combine_bam_ind running time: ", running_time)
  
  # Merge BAM files by cell type and individual
  running_time <- combine_bam_celltype_ind(
    bam_cell_anno_ = bam_data_$bam_cell_anno,
    bam_splited_dir_ = file.path(out_dir, "bam_split_by_celltype_ind_split"),
    out_dir_ = file.path(out_dir, "bam_split_by_celltype_ind"),
    parThreads_ = parThreads,
    subThreads_ = 4
  )
  message("combine_bam_celltype_ind running time: ", running_time)
  fwrite(data.table(type = "combine_bam_celltype_ind", time = running_time), 
         file.path(out_dir, "log.txt"), sep = "\t", append = TRUE)
  
  
  # Identify RNA editing sites using REDItools
  running_time <- detect_editing_by_reditools_13m_noSNP_for_ind(
    bam_cell_anno_ = bam_data_$bam_cell_anno,
    bam_ind_dir_ = file.path(out_dir, "bam_split_by_ind"),
    RE_out_dir_ = file.path(out_dir, "reditools_13m_ind"),
    REDIportal_tab_ = REDIportal_tab,
    genome_fa_ = genome_fa,
    min_cov_ = 5,
    read_strand_ = read_strand,
    min_edited_reads = 2,
    min_base_quality_ = 30,
    min_mapping_quality = 30,
    parThreads_ = ifelse(round(parThreads/10) >= 1, round(parThreads/10), 1),
    subThreads_ = ifelse(parThreads >= 10, 10, parThreads),
    python2.7_exe_ = py2_path
  )
  message("detect_editing_by_reditools_13m_noSNP_for_ind running time: ", running_time)
  fwrite(data.table(type = "detect_editing", time = running_time), 
         file.path(out_dir, "log.txt"), sep = "\t", append = TRUE)
  
  
  # Calculate coverage statistics per cell type
  running_time <- summary_ed_cov_cellNum_for_cell_type_ind_pool(
    bam_cell_anno_ = bam_data_$bam_cell_anno,
    RE_out_reditools_ = file.path(out_dir, "reditools_13m_ind"),
    bam_split_by_celltype_ind_path = file.path(out_dir, "bam_split_by_celltype_ind"),
    ref_fa_ = genome_fa,
    out_dir_ = file.path(out_dir, "reditools_13m_ind_celltype"),
    bed_splited_lines = 1000,
    parThreads_ = parThreads
  )
  message("summary_ed_cov running time: ", running_time)
  fwrite(data.table(type = "summary_ed_cov", time = running_time), 
         file.path(out_dir, "log.txt"), sep = "\t", append = TRUE)
  
  
  # Apply coverage-based filtering
  running_time <- filter_ed_by_cell_and_cov_ind(
    bam_cell_anno_ = bam_data_$bam_cell_anno,
    RE_out_reditools_ = file.path(out_dir, "reditools_13m_ind"),
    mmp_out_dir = file.path(out_dir, "reditools_13m_ind_celltype"),
    min_num_covered_cells = 1,
    min_cov_reads = 1,
    min_num_covered_cells_0 = 1,
    min_cov_reads_0 = 1,
    parThreads = parThreads
  )
  message("filter_ed_by_cell running time: ", running_time)
  
  # Combine results across individuals
  running_time <- merge_eds_ind_filtered(
    bam_cell_anno_ = bam_data_$bam_cell_anno,
    mmp_out_dir = file.path(out_dir, "reditools_13m_ind_celltype"),
    RE_dcast_out_dir_ = file.path(out_dir, "reditools_13m_ind_merged"),
    prop_of_ed_shared_inds = 0.1,
    parThreads_merge_ = ifelse(parThreads > 10, 10, parThreads)
  )
  message("merge_eds_ind running time: ", running_time)
  
  # Return final output directory path
  return(out_dir)
}