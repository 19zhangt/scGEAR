#!/src/bash

# PAS (PolyA Site) Identification Pipeline
# Description: A workflow for identifying polyadenylation sites using multiple data PAS_sources
# Usage: ./pas_pipeline.sh <BAMFILELIST> <BAMFILELIST> <THREAD_NUM> <REFGENOME> <OUTDIR>


# Set global variables

# BAMFILELIST=example/APA/bam_list.csv
# THREAD_NUM=4
# REFGENOME="/media/iceland/share/Index/Genome_index/Human_hg38/GRCh38.primary_assembly.genome.fa"
# OUTDIR=example/APA

BAMFILELIST=$1
THREAD_NUM=$2
REFGENOME=$3
OUTDIR=$4

# Create directory structure
mkdir -p "${OUTDIR}/PAS_source/"{refSeq,polyA_DB,mws_cs_annotation,findpolyAtail}
RSDIR="${OUTDIR}/PAS_source/refSeq"
MWSDIR="${OUTDIR}/PAS_source/mws_cs_annotation"
POLYADBDIR="${OUTDIR}/PAS_source/polyA_DB"
TAILDIR="${OUTDIR}/PAS_source/findpolyAtail"

# --------------------------
# 1. RefSeq Processing
# --------------------------
process_refseq() {
    # Download and prepare RefSeq data
    if [[ ! -f ${RSDIR}/GCF_000001405.40_GRCh38.p14_genomic.gtf || ! -f ${RSDIR}/GCF_000001405.40_GRCh38.p14_assembly_report.txt ]];then
        wget -q -P "${RSDIR}" \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt

        gunzip -f "${RSDIR}/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
        dos2unix -q "${RSDIR}/GCF_000001405.40_GRCh38.p14_assembly_report.txt"
    fi

    # Process chromosome names
    awk -F"\t" 'NR==FNR{if($0!~/#/&&$7~/NC_/){a[$7]=$NF};next}
    {OFS="\t";if($0~/#/){print $0}else if(a[$1]){$1=a[$1];print $0}}' \
    "${RSDIR}/GCF_000001405.40_GRCh38.p14_assembly_report.txt" "${RSDIR}/GCF_000001405.40_GRCh38.p14_genomic.gtf" > "${RSDIR}/GRCh38_rechr.gtf"

    # Extract and process last exons
    python "src/extract_last_exon.py" "${RSDIR}/GRCh38_rechr.gtf" "${RSDIR}/last_exon_each_transcript.gtf"
    
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
    Rscript src/refseq.annotation.R ${RSDIR}/GRCh38_rechr.gtf ${RSDIR}/refseq.annotation.txt
    awk -F"\t" 'NR==FNR{a[$5]=$5;next}{if(a[$4]){print $0}}' ${RSDIR}/refseq.annotation.txt ${RSDIR}/last_exon.bed > ${RSDIR}/last_exon.genes.bed
    rm ${RSDIR}/last_exon.bed

    Rscript src/last_exon_rename.R ${RSDIR}/last_exon.genes.bed ${RSDIR}/last_exon_annotation.bed
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
    }' "${RSDIR}/GRCh38_rechr.gtf" | sort -k1,1 -k2,2n | uniq > "${RSDIR}/refseq_tts.bed"
    bedtools intersect -wo -s -a ${RSDIR}/refseq_tts.bed -b ${RSDIR}/last_exon_annotation.bed > ${RSDIR}/refseq_tts_in_last_exon.bed

    # Process external databases
    bedtools intersect -wo -s -a reference/UCSC_polyDB_hg38.bed \
        -b "${RSDIR}/last_exon_annotation.bed" > "${POLYADBDIR}/UCSC_polyDB_in_last_exon.bed"
    
    # Process MWS annotations
    zcat "reference/mws_cs_celltypes_hg38.bed.gz" | awk -F"\t" '{
        OFS="\t";
        $3=$3+1;
        print "chr"$0
    }' | \
    bedtools intersect -wo -s -a - -b "${RSDIR}/last_exon_annotation.bed" | \
    cut -f 1-6 > "${MWSDIR}/mws_cs_last_exon.bed"

    Rscript "src/mws_cs_site_collapse.R" "${MWSDIR}" "${REFGENOME}" "${THREAD_NUM}"
}

# --------------------------
# 3. PolyA Tail Detection
# --------------------------
detect_polya_tails() {
    mkdir -p "${TAILDIR}/site_report"

    # Run with GNU Parallel
    tail -n+2 "${BAMFILELIST}" | parallel --colsep ',' --jobs "${THREAD_NUM}" --progress \
        python "src/findAtail_from_bam.py" {1} {2} "${TAILDIR}"

    # Post-processing
    Rscript "src/findAtail_collection.R" "${TAILDIR}"

    bedtools intersect -wo -s -a "${TAILDIR}/findtail_sites.bed" -b "${RSDIR}/last_exon_annotation.bed" \
        > "${TAILDIR}/findtail_sites_in_last_exon.bed"
    
    Rscript "src/findAtail_collapse.R" "${TAILDIR}" "${REFGENOME}" "${THREAD_NUM}"
}


integrate_data() {
    # Final integration
    Rscript "src/PAS_identification.R" "${OUTDIR}" "${THREAD_NUM}"
}

# --------------------------
# Main Execution Flow
# --------------------------
main() {
    # process_refseq
    # RefSeq_polyDB_mws
    detect_polya_tails
    integrate_data
    
}

# Execute main function with timing
start_time=$(date +%s)
main
end_time=$(date +%s)
echo "Total execution time: $((end_time - start_time)) seconds"
