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

#' Initialize and configure command line arguments
#' @return ArgumentParser object with configured parameters
initialize_parser <- function() {
    parser <- ArgumentParser(description = "scGEAR: a population-scale computational framework for identifying multiple 
                             molecular phenotypes from polyA-enriched single-cell transcriptome data")
    
    # Required parameters
    parser$add_argument("--cellinfo", type = "character", required = TRUE,
                        help = "Unique cell barcodes assigned to each sample within each cluster.")
    parser$add_argument("--cellrangerDir", type = "character", required = TRUE,
                        help = "Directory containing Cell Ranger outputs for each pool.")
    
    # # Tool paths with validation
    # parser$add_argument("--umitools_path", type = "character", default = "umi_tools",
    #                     help = "Path to umi_tools executable.")
    # parser$add_argument("--samtools_path", type = "character", default = "/usr/bin/samtools",
    #                     help = "Path to samtools executable.")
    # parser$add_argument("--sinto_path", type = "character", default = "sinto",
    #                     help = "Path to sinto executable.")
    
    # Reference files
    parser$add_argument("--genome_fa", type = "character", required = TRUE,
                        help = "Path to genome FASTA file.")
    parser$add_argument("--REDIportal_tab", type = "character", required = TRUE,
                        help = "Path to REDIportal annotation file (tab.gz format).")
    
    # Processing parameters
    parser$add_argument("--read_strand", type = "integer", choices = c(0, 1), default = 1,
                        help = "Strand orientation (0=unstranded, 1=stranded)")
    parser$add_argument("--mode", type = "character", choices = c("exp", "apa", "ed", "all"), default = "exp",
                        help = "Processing mode: exp (expression), apa (alternative polyadenylation), ed (RNA editing), all")
    parser$add_argument("--threads", type = "integer", default = 4,
                        help = "Number of parallel processing threads")
    parser$add_argument("--outdir", type = "character", required = TRUE,
                        help = "Output directory path")
    
    return(parser)
}

#' Initialize and configure command line arguments
#' @return ArgumentParser object with configured parameters
initialize_parser <- function() {
    parser <- ArgumentParser(description = "scGEAR: A Population-scale Computational Framework for identifying multiple transcriptional phenotypes from polyA-enriched Single-cell Transcriptome Data")
    
    # Required parameters
    parser$add_argument("--cellinfo", type = "character", required = TRUE,
                        help = "Unique cell barcodes assigned to each sample within each cluster.")
    parser$add_argument("--cellrangerDir", type = "character", required = TRUE,
                        help = "Directory containing Cell Ranger outputs for each pool.")
    
    # Tool paths with validation
    parser$add_argument("--umitools_path", type = "character", default = "umi_tools",
                        help = "Path to umi_tools executable.")
    parser$add_argument("--samtools_path", type = "character", default = "/usr/bin/samtools",
                        help = "Path to samtools executable.")
    parser$add_argument("--sinto_path", type = "character", default = "sinto",
                        help = "Path to sinto executable.")
    
    # Reference files
    parser$add_argument("--genome_fa", type = "character", required = TRUE,
                        help = "Path to genome FASTA file.")
    parser$add_argument("--REDIportal_tab", type = "character", required = TRUE,
                        help = "Path to REDIportal annotation file (tab.gz format).")
    
    # Processing parameters
    parser$add_argument("--read_strand", type = "integer", choices = c(0, 1), default = 1,
                        help = "Strand orientation (0=unstranded, 1=stranded)")
    parser$add_argument("--mode", type = "character", choices = c("exp", "apa", "ed", "all"), default = "exp",
                        help = "Processing mode: exp (expression), apa (alternative polyA), ed (RNA editing), all")
    parser$add_argument("--threads", type = "integer", default = 4,
                        help = "Number of parallel processing threads")
    parser$add_argument("--outdir", type = "character", required = TRUE,
                        help = "Output directory path")
    
    return(parser)
}

# Main execution flow ---------------------------------------------------------
main <- function() {
    # Parse command line arguments
    parser <- initialize_parser()
    args <- parser$parse_args()
    
    
    args <- list(
        cellinfo = normalizePath(cellinfo),
        cellrangerDir = normalizePath(cellrangerDir),
        umitools_path = umitools_path,
        samtools_path = samtools_path,
        py2_path = py2_path,
        sinto_path = sinto_path,
        genome_fa = genome_fa,
        REDIportal_tab = REDIportal_tab,
        read_strand = as.integer(read_strand),
        mode = mode,
        threads = as.integer(threads),
        outdir = outdir
    )
    
    
    # # Validate tool paths
    # validate_executable <- function(path, name) {
    #     if (!file.exists(path)) {
    #         stop("Path for ", name, " not found: ", path)
    #     }
    #     if (file.access(path, 1) == -1) {
    #         stop("No execute permission for ", name, " at: ", path)
    #     }
    # }
    # 
    # validate_executable(args$samtools_path, "samtools")
    # validate_executable(args$umitools_path, "umi_tools")
    # validate_executable(args$sinto_path, "sinto")
    
    # Create output directory
    if (!dir.exists(args$outdir)) {
        message("Creating output directory: ", args$outdir)
        dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Load processing functions
    script_dir <- get_script_directory()
    functions_path <- file.path(script_dir, "src", "functions.R")
    if (!file.exists(functions_path)) {
        stop("Functions script not found at: ", functions_path)
    }
    source(functions_path)
    
    # Load and validate input data
    validate_input_file(args$cellinfo, "Cell information file")
    validate_input_dir(args$cellrangerDir, "Cell Ranger directory")
    
    # Process input data
    bam_data <- load_bam_data(args$cellinfo, args$cellrangerDir)
    
    # Main processing based on mode
    process_mode(args$mode, bam_data, args)
}

#' Retrieve script directory path
#' @return Absolute path to script directory
get_script_directory <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    script_path <- sub("^--file=", "", cmd_args[grepl("^--file=", cmd_args)])
    if (length(script_path) == 0) return(".")
    return(normalizePath(dirname(script_path)))
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

#' Load and validate BAM data
#' @param cellinfo_path Path to cell information file
#' @param cellranger_dir Path to Cell Ranger directory
#' @return List containing validated BAM data
load_bam_data <- function(cellinfo_path, cellranger_dir) {
    message("Loading cell information from: ", cellinfo_path)
    bam_anno <- fread(cellinfo_path)
    
    if (nrow(bam_anno) == 0) {
        stop("Empty cell information file: ", cellinfo_path)
    }
    
    message("Validating BAM files in: ", cellranger_dir)
    # Add BAM file validation logic here
    
    return(list(
        bam_anno = bam_anno,
        cell_count = nrow(bam_anno)
    ))
}

#' Execute processing based on selected mode
#' @param mode Processing mode
#' @param bam_data List containing BAM data
#' @param args Parsed command line arguments
process_mode <- function(mode, bam_data, args) {
    log_file <- file.path(args$outdir, "processing_log.txt")
    
    # Common preprocessing steps
    if (mode %in% c("ed", "all")) {
        run_duplicate_removal(bam_data, args, log_file)
        run_bam_processing(bam_data, args, log_file)
    }
    
    # Mode-specific processing
    switch(mode,
           "ed" = {
               run_editing_analysis(bam_data, args, log_file)
           },
           "apa" = {
               run_apa_analysis(bam_data, args, log_file)
           },
           "all" = {
               run_full_analysis(bam_data, args, log_file)
           },
           "exp" = {
               run_expression_analysis(bam_data, args, log_file)
           }
    )
}

# Execute main function
if (!interactive()) {
    main()
}