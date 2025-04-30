#!/usr/bin/env Rscript
# SERA Main Processing Pipeline
# Author: Ting Zhang and Ruofan Ding
# Date: 2025-04-20

# Load required packages with error handling
required_packages <- c("BiocManager", "devtools", "data.table", "pbmcapply", 
                       "vroom", "Hmisc", "reshape2", "seqinr", "stringr", 
                       "DescTools", "argparse", "Matrix", "Seurat")
suppressPackageStartupMessages({
    for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg)
        }
    }
    if (!require("bedtoolsr", character.only = TRUE)) {
        devtools::install_github("PhanstielLab/bedtoolsr")
    }
    if (!require("GenomicRanges", character.only = TRUE)) {
        BiocManager::install("GenomicRanges")
    }
    if (!require("rtracklayer", character.only = TRUE)) {
        BiocManager::install("rtracklayer")
    }
    library(argparse)
    library(data.table)
    library(stringr)
})


# Defining the conversion function
convert_logical <- function(x) {
    x_upper <- toupper(x)
    if (x_upper %in% c("TRUE", "T")) {
        return(TRUE)
    } else if (x_upper %in% c("FALSE", "F")) {
        return(FALSE)
    } else {
        stop("The parameter must be TRUE/T or FALSE/F")
    }
}


#' Initialize and configure command line arguments
#' @return ArgumentParser object with configured parameters
initialize_parser <- function() {
    parser <- ArgumentParser(description = "SERA: a computational framework for analyzing multiple RNA-based phenotypes in scRNA-seq data")
    
    parser$add_argument("--mode", type = "character", choices = c("expression", "APA", "editing", "all"), default = "expression",
                        help = "Processing mode")
    
    # Required parameters
    parser$add_argument("--cellinfo", type = "character", required = TRUE,
                        help = "Unique cell barcodes assigned to each sample within each cluster.")
    parser$add_argument("--baminfo", type = "character", required = TRUE,
                        help = "Directory containing Cell Ranger outputs for each pool.")
    parser$add_argument("--remove_duplicates", type=convert_logical, default = FALSE,
                        help = "Logical parameters, accept TRUE/T or FALSE/F (case insensitive).")
    
    # Reference file
    parser$add_argument("--genome_fa", type = "character", required = TRUE,
                        help = "Path to genome FASTA file.")
    
    # Processing parameters
    parser$add_argument("--REDIportal_anno", type = "character", required = TRUE,
                        help = "Path to REDIportal annotation file (tab.gz format).")
    parser$add_argument("--read_strand", type = "integer", choices = c(0, 1), default = 1,
                        help = "Strand orientation (0=unstranded, 1=stranded)")
    
    parser$add_argument("--threads", type = "integer", default = 4,
                        help = "Number of parallel processing threads")
    parser$add_argument("--outdir", type = "character", required = TRUE,
                        help = "Output directory path")
    return(parser)
}


#' Validate input file existence
#' @param path Path to validate
#' @param name Human-readable name for error messages
validate_input_file <- function(path, name) {
    if (!file.exists(path)) {
        stop(name, " not found at: ", path)
    }
}

#' Validate input directory existence
#' @param path Path to validate
#' @param name Human-readable name for error messages
validate_input_dir <- function(path, name) {
    if (!dir.exists(path)) {
        stop(name, " not found at: ", path)
    }
}


#' Execute processing based on selected mode
#' @param gear Processing mode
#' @param bam_data List containing BAM data
#' @param args Parsed command line arguments
process_mode <- function(gear, bam_data_, args_) {
    log_file <- file.path(args_$outdir, "processing_log.txt")
    
    # Common preprocessing steps
    if (gear %in% c("editing", "all", "APA")) {
        if (system2("umi_tools", args = "--version", stdout = NULL, stderr = NULL) != 0) {
            stop("umi_tools is not installed! Please install via pip: pip install umi_tools")
        }
        if (system2("samtools", args = "--version", stdout = NULL, stderr = NULL) != 0) {
            stop("samtools is not installed! Please install via Conda: conda install bioconda::samtools")
        }
    }

    if(args_$remove_duplicates){
        ## remove duplicates for each raw pooled bam file
        bam_data_$bam_anno <- remove_dups_pools2(bam_anno_ = bam_data_$bam_anno, 
                                                parThreads_ = args_$threads, 
                                                out_path = args_$outdir, log_ = log_file)
    }else{
        colnames(bam_data_$bam_anno)[2] <- "bam_dedup_path"
    }
    
    # Mode-specific processing
    switch(gear,
           "editing" = {
               args_$outdir <- paste0(args_$outdir, "/editing")
               # Load processing functions
               source("src/mode_editing.R")
               run_editing_analysis(bam_data_ = bam_data_, args_ = args_, log_ = log_file)
           },
           "APA" = {
               args_$outdir <- paste0(args_$outdir, "/", gear)
               source("src/mode_APA.R")
               run_apa_analysis(bam_data_ = bam_data_, args_ = args_, log_ = log_file)
           },
           "expression" = {
               source("src/mode_expression.R")
               args_$outdir <- paste0(args_$outdir, "/", gear)
               run_expression_analysis(bam_data_ = bam_data_, args_ = args_, log_ = log_file)
           },
           "all" = {
               dir_output <- args_$outdir
               args_$outdir <- paste0(dir_output, "/expression")
               source("src/mode_expression.R")
               run_expression_analysis(bam_data_, args_, log_file)
               
               args_$outdir <- paste0(dir_output, "/APA")
               source("src/mode_APA.R")
               run_apa_analysis(bam_data_, args_, log_file)
               
               args_$outdir <- paste0(dir_output, "/editing")
               source("src/mode_editing.R")
               run_editing_analysis(bam_data_, args_, log_file)
           }
    )
}


args <- list(
    mode = "all",
    cellinfo = "example/cell_annotation.csv",
    baminfo = "example/bam_list.csv",
    remove_duplicates = TRUE,
    genome_fa = "/media/iceland/share/Index/Genome_index/Human_hg38/GRCh38.primary_assembly.genome.fa",
    REDIportal_anno = "reference/REDIportal_noRSid.tab.gz",
    py2_path = "python2",
    read_strand = as.integer("1"),
    threads = as.integer("4"),
    outdir = "example"
)


# Main execution flow ---------------------------------------------------------
main <- function() {
    # Parse command line arguments
    parser <- initialize_parser()
    args <- parser$parse_args()
    
    # Parameter calibration
    if (!args$mode %in% c("expression", "APA", "editing", "all")) {
        stop("Usage: Rscript SERA.R --mode <mode>\n  mode options: expression / APA / editing / all")
    }
    
    # Define the folder
    target_dirs <- switch(args$mode,
                          "expression" = paste0(args$outdir, "/expression"),
                          "APA" = paste0(args$outdir, "/APA"),
                          "editing" = paste0(args$outdir, "/editing"),
                          "all" = paste0(args$outdir, c("/expression", "/APA", "/editing"))
    )
    
    # Interactive deletion function
    confirm_removal <- function(dir_path) {
        response <- tolower(readline(prompt = paste0(
            dir_path, " already exists. Delete existing directory? (y/n)"
        )))
        
        if (response == "y") {
            unlink(dir_path, recursive = TRUE)
            message("Deleted directories: ", dir_path)
            return(TRUE)
        } else if (response == "n") {
            message("Keep the existing directory and exit the script")
            return(FALSE)
        } else {
            stop("Invalid input, please enter y or n")
        }
    }
    
    # Main process control
    for (dir in target_dirs) {
        if (dir.exists(dir)) {
            message("\n>>> Conflict Detection <<<")
            success <- confirm_removal(dir)
            if (!success) quit(save = "no", status = 1)
        }
        # Create output directory
        dir.create(dir, recursive = TRUE)
        if (dir.exists(dir)) {
            message("Creating output directory: ", dir)
        } else {
            stop("Directory creation failed: ", dir)
        }
    }
    
    # Load and validate input data
    validate_input_file(args$cellinfo, "Cell information file")
    validate_input_file(args$baminfo, "Cell Ranger bam files")
    
    # Process input data
    source("src/utility_function.R")
    bam_data <- load_bam_cell_anno3(bam_anno_file_ = args$cellinfo, 
                                    bam_file_ = args$baminfo, 
                                    check_bam = TRUE)
    
    # Main processing based on mode
    process_mode(gear = args$mode, bam_data_ = bam_data, args_ = args)
}

# Execute main function
if (!interactive()) {
    main()
}