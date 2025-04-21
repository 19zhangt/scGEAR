#!/usr/bin/env Rscript
# scGEAR Main Processing Pipeline
# Author: Ting Zhang and Ruofan Ding
# Date: 2025-04-20


# Load required packages with error handling
required_packages <- c("data.table", "stringr", "parallel", "argparse")
suppressPackageStartupMessages({
    for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE)) {
            stop("Package ", pkg, " not found. Please install first.", call. = FALSE)
        }
    }
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
    parser <- ArgumentParser(description = "scGEAR: a computational framework for analyzing multiple RNA-based phenotypes in scRNA-seq data")
    
    parser$add_argument("--gear", type = "character", choices = c("expression", "apa", "editing", "all"), default = "expression",
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


args <- list(
    gear = "apa",
    cellinfo = "example/cell_annotation.csv",
    baminfo = "example/bam_list.csv",
    remove_duplicates = TRUE,
    genome_fa = "/media/iceland/share/Index/Genome_index/Human_hg38/GRCh38.primary_assembly.genome.fa",
    REDIportal_anno = "ref/REDIportal_noRSid.tab.gz",
    py2_path = "python2",
    read_strand = as.integer("1"),
    threads = as.integer("4"),
    outdir = "example/output"
)


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
process_mode <- function(gear, bam_data, args) {
    log_file <- file.path(args$outdir, "processing_log.txt")
    
    # Common preprocessing steps
    if (gear %in% c("editing", "all", "apa")) {
        if (system2("umi_tools", args = "--version", stdout = NULL, stderr = NULL) != 0) {
            stop("umi_tools is not installed! Please install via pip: pip install umi_tools")
        }
        if (system2("samtools", args = "--version", stdout = NULL, stderr = NULL) != 0) {
            stop("samtools is not installed! Please install via Conda: conda install bioconda::samtools")
        }
        
        if(args$remove_duplicates){
            ## remove duplicates for each raw pooled bam file
            bam_data$bam_anno <- remove_dups_pools2(bam_anno_ = bam_data$bam_anno, 
                                                    parThreads_ = args$threads, 
                                                    out_path = args$outdir, log_ = log_file)
        }else{
            colnames(bam_data$bam_anno)[2] <- "bam_dedup_path"
        }
    }

    # Mode-specific processing
    switch(gear,
           "editing" = {
               args$outdir <- paste0(args$outdir, "/", args$gear)
               # Load processing functions
               functions_path <- file.path("src", "RNAediting.R")
               if (!file.exists(functions_path)) {
                   stop("Functions script not found at: ", functions_path)
               }
               source(functions_path)
               run_editing_analysis(bam_data, args, log_file)
           },
           "apa" = {
               # Load processing functions
               functions_path <- file.path("src", "apa.R")
               if (!file.exists(functions_path)) {
                   stop("Functions script not found at: ", functions_path)
               }
               source(functions_path)
               
               args$outdir <- paste0(args$outdir, "/", args$gear)
               run_apa_analysis(bam_data, args, log_file)
           },
           "exp" = {
               args$outdir <- paste0(args$outdir, "/", args$gear)
               run_expression_analysis(bam_data, args, log_file)
           },
           "all" = {
               # Load processing functions
               functions_path <- file.path("src", "functions.R")
               if (!file.exists(functions_path)) {
                   stop("Functions script not found at: ", functions_path)
               }
               source(functions_path)
               
               dir_output <- args$outdir
               args$outdir <- paste0(dir_output, "/expression")
               run_expression_analysis(bam_data, args, log_file)
               
               args$outdir <- paste0(dir_output, "/apa")
               run_apa_analysis(bam_data, args, log_file)
               
               args$outdir <- paste0(dir_output, "/editing")
               run_editing_analysis(bam_data, args, log_file)
           }
    )
}

# Main execution flow ---------------------------------------------------------
main <- function() {
    # Parse command line arguments
    parser <- initialize_parser()
    args <- parser$parse_args()
    
    # Create output directory
    if (!dir.exists(args$outdir)) {
        message("Creating output directory: ", args$outdir)
        dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Load and validate input data
    validate_input_file(args$cellinfo, "Cell information file")
    validate_input_file(args$baminfo, "Cell Ranger bam files")
    
    # Process input data
    source("src/base.R")
    bam_data <- load_bam_cell_anno3(args$cellinfo, args$baminfo, check_bam = TRUE)
    
    # Main processing based on mode
    process_mode(args$gear, bam_data, args)
}

# Execute main function
if (!interactive()) {
    main()
}