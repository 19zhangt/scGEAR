#!/bin/bash

# PAS (PolyA Site) Identification Pipeline
# Description: A workflow for identifying polyadenylation sites using multiple data sources
# Usage: ./pas_pipeline.sh <CURRENT_DIR> <BAMFILELIST> <THREAD_NUM> <REFGENOME>

# --------------------------
# Initialize Parameters
# --------------------------
# Check input parameters
if [ $# -ne 4 ]; then
    echo "Usage: $0 <CURRENT_DIR> <BAMFILELIST> <THREAD_NUM> <REFGENOME>"
    exit 1
fi


# Set global variables
CURRENT_DIR="/media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/scGEAR/ref/PAS_sites"
BAMFILELIST=${CURRENT_DIR}/src/bamfiles.csv
THREAD_NUM=2
REFGENOME="/media/iceland/share/Index/Genome_index/Human_hg38/GRCh38.primary_assembly.genome.fa"

# CURRENT_DIR=$1
# BAMFILELIST=$2
# thread_num=$3
# REFGENOME=$4
RSDIR="${CURRENT_DIR}/data/Refseq"
PATDIR="${CURRENT_DIR}/data/findpolyAtail"

# Create directory structure
mkdir -p "${CURRENT_DIR}/data/"{Refseq,polyA_DB,mws_cs_annotation,findpolyAtail}

# --------------------------
# 1. RefSeq Processing
# --------------------------
process_refseq() {
    # Download and prepare RefSeq data
    wget -q -P "${RSDIR}" \
        https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
        https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
    
    gunzip -f "${RSDIR}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
    dos2unix -q "${RSDIR}/GCF_000001405.40_GRCh38.p14_assembly_report.txt"

    # Process chromosome names
    awk -F"\t" '
    NR==FNR{if($0!~/#/&&$7~/NC_/){a[$7]=$NF};next}
    {OFS="\t";if($0~/#/){print $0}else if(a[$1]){$1=a[$1];print $0}}' \
    "${RSDIR}/GCF_000001405.40_GRCh38.p14_assembly_report.txt" "${RSDIR}/GCF_000001405.40_GRCh38.p14_genomic.gtf" > "${RSDIR}/GRCh38_rechr.gtf"

    # Install required Python packages
    declare -A PKGS=(
        ["cgat"]="git+https://github.com/CGATOxford/cgat.git"
        ["cgatcore"]="git+https://github.com/CGATOxford/cgat-core.git"
    )
    
    for pkg in "${!PKGS[@]}"; do
        python3 -m pip show "${pkg}" &>/dev/null || \
        python3 -m pip install --user "${PKGS[$pkg]}"
    done

    # Extract and process last exons
    python "${CURRENT_DIR}/bin/extract_last_exon.py" "${RSDIR}/GRCh38_rechr.gtf" > "${RSDIR}/last_exon_each_transcript.gtf"
    
    awk -F"\t" '{
        OFS="\t";
        split($9,a," ");
        print $1, $4-1, $5, a[2], ".", $7, a[4]
    }' "${RSDIR}/last_exon_each_transcript.gtf" | sed "s/[\";]//g" | sort -k1,1 -k2,2n | uniq > "${RSDIR}/last_exon_each_transcript.bed"

    # Additional processing steps...
    ## extract last exon for each gene
    awk '{tid = $4;
        if (chain[tid] == "") {
            chain[tid] = $6;
        }
        if ($6 == "+") {
            if (last_exon_end[tid] < $3) {
                last_exon_end[tid] = $3;
                last_exon_line[tid] = $0;
            }
        } else {
            if (last_exon_start[tid] == "" || $2 < last_exon_start[tid]) {
                last_exon_start[tid] = $2;
                last_exon_line[tid] = $0;
            }
        }
    }
    END {
        for (tid in last_exon_line) {
            print last_exon_line[tid];
        }
    }' ${RSDIR}/last_exon_each_transcript.bed > ${RSDIR}/last_exon_each_gene.bed

    awk -F "\t" '$0!~/#/&&$3=="gene"{OFS="\t";split($9,a," ");
        print $1, $4-1, $5, a[2], 0, $7
    }' ${RSDIR}/GRCh38_rechr.gtf | sed "s/\"//g;s/;//" | sort -k1,1 -k2,2n > ${RSDIR}/genes.bed

    extend_length="2000"
    bedtools slop -i ${RSDIR}/last_exon_each_gene.bed -g ${REFGENOME}.fai -l 0 -r ${extend_length} -s > ${RSDIR}/last_exon_last_expanded.bed
    bedtools intersect -wao -a ${RSDIR}/last_exon_last_expanded.bed -b ${RSDIR}/genes.bed | awk '{if($4!=$11){print $0}}' > ${RSDIR}/last_exon_last_expanded_overlap.txt

    awk 'NR==FNR{a[$4]=$4;next}{if(!a[$4]){print $0}}' ${RSDIR}/last_exon_last_expanded_overlap.txt ${RSDIR}/last_exon_last_expanded.bed > ${RSDIR}/last_exon_last_expanded_nooverlap.bed

    cat ${RSDIR}/last_exon_each_transcript.bed ${RSDIR}/last_exon_last_expanded_nooverlap.bed | sort -k1,1 -k2,2n >${RSDIR}/last_exon.bed
    rm ${RSDIR}/last_exon_*

    ## last exon within genes
    Rscript ${CURRENT_DIR}/bin/refseq.annotation.R ${RSDIR}/GRCh38_rechr.gtf ${RSDIR}/refseq.annotation.txt
    awk -F"\t" 'NR==FNR{a[$5]=$5;next}{if(a[$4]){print $0}}' ${RSDIR}/refseq.annotation.txt ${RSDIR}/last_exon.bed > ${RSDIR}/last_exon.genes.bed

    ## merge overlapped last exons
    sort -k 1,1 -k 2,2n ${RSDIR}/last_exon.genes.bed | bedtools merge -s -c 4,5,6,7 -o distinct,distinct,distinct,distinct -i - > ${RSDIR}/last_exon.merge.bed
    Rscript ${CURRENT_DIR}/bin/last_exon_rename.R ${RSDIR}/last_exon.merge.bed ${RSDIR}/last_exon.genes.bed ${RSDIR}/last_exon_annotation.bed

    rm ${RSDIR}/last_exon.*
}


# --------------------------
# 2. RefSeq_polyDB_mws Annotation
# --------------------------
RefSeq_polyDB_mws() {
    # Process RefSeq TTS
    awk -F"\t" '$3=="transcript" {
        OFS="\t";
        if($7=="+") {print $1, $5-1, $5, ".", ".", $7}
        else {print $1, $4-1, $4, ".", ".", $7}
    }' "${RSDIR}/GRCh38_rechr.gtf" | sort -k1,1 -k2,2n > "${RSDIR}/refseq_tts.bed"
    bedtools intersect -wo -s -a ${RSDIR}/refseq_tts.bed -b ${RSDIR}/last_exon_annotation.bed > ${RSDIR}/refseq_tts_in_last_exon.bed

    # Process external databases
    bedtools intersect -wo -s -a "${CURRENT_DIR}/src/UCSC_polyDB_hg38.bed" \
        -b "${RSDIR}/last_exon_annotation.bed" > "${CURRENT_DIR}/data/polyA_DB/UCSC_polyDB_in_last_exon.bed"
    
    # Process MWS annotations
    awk -F"\t" '{
        OFS="\t";
        $3=$3+1;
        print "chr"$0
    }' "${CURRENT_DIR}/src/mws_cs_celltypes_hg38.bed" | \
    bedtools intersect -wo -s -a - -b "${RSDIR}/last_exon_annotation.bed" | \
    cut -f 1-6 > "${CURRENT_DIR}/data/mws_cs_annotation/mws_cs_last_exon.bed"

    Rscript ${CURRENT_DIR}/bin/mws_cs_site_aggr.R ${CURRENT_DIR}/data/mws_cs_annotation ${REFGENOME}
}

# --------------------------
# 3. PolyA Tail Detection
# --------------------------
detect_polya_tails() {
    # mkdir -p "${PATDIR}/site_report"

    # # Run with GNU Parallel
    # parallel --colsep ',' --jobs "${THREAD_NUM}" --progress \
    #     python "${CURRENT_DIR}/bin/findpolyAtail.py" {1} {2} "${PATDIR}" \
    #     < "${BAMFILELIST}"

    # # Post-processing
    # Rscript "${CURRENT_DIR}/bin/findtail_for_annotation.R" "${PATDIR}"

    bedtools intersect -wo -s -a "${PATDIR}/findtail_sites.bed" -b "${RSDIR}/last_exon_annotation.bed" \
        > "${PATDIR}/findtail_sites_in_last_exon.bed"
    
    Rscript "${CURRENT_DIR}/bin/findtail_merge.R" "${PATDIR}" "${REFGENOME}"
}


integrate_data() {
    # Final integration
    Rscript "${CURRENT_DIR}/bin/PAS_identification.R" "${CURRENT_DIR}"
}

# --------------------------
# Main Execution Flow
# --------------------------
main() {
    # process_refseq
    # RefSeq_polyDB_mws
    detect_polya_tails
    # integrate_data
}

# Execute main function with timing
start_time=$(date +%s)
main
end_time=$(date +%s)
echo "Total execution time: $((end_time - start_time)) seconds"
