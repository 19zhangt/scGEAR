import pybedtools

# load exons and genes
genes = pybedtools.BedTool('/media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/scGEAR/ref/PAS_sites/Refseq_last_exon/genes.bed')
exons = pybedtools.BedTool('/media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/scGEAR/ref/PAS_sites/Refseq_last_exon/last_exons_last_expanded.bed')

for exon in exons:
    chr_exon = exon.chrom
    start_exon = exon.start
    end_exon = exon.end
    strand_exon = exon.strand
    exon_name = exon.name

    if strand_exon == '+':
        ext_start = end_exon
        ext_end = end_exon + 2000
    else:
        ext_start = start_exon - 2000
        ext_end = start_exon

    extension = pybedtools.BedTool(
        f"{chr_exon}\t{ext_start}\t{ext_end}\tExtension\t0\t{strand_exon}",
        from_string=True
    )

    overlapping = genes.intersect(extension, wa=True).filter(
        lambda x: x.name != exon_name
    )

    min_distance = None

    for gene in overlapping:
        gene_start = gene.start
        gene_end = gene.end
        gene_strand = gene.strand

        if gene_strand == '+':
            downstream_point = gene_start
        else:
            downstream_point = gene_end

        if strand_exon == '+':
            if downstream_point < end_exon:
                continue
            distance = downstream_point - end_exon
        else:
            if downstream_point > start_exon:
                continue
            distance = start_exon - downstream_point

        if distance >= 0 and (min_distance is None or distance < min_distance):
            min_distance = distance

    if min_distance is not None:
        if strand_exon == '+':
            adjusted_end = end_exon + (min_distance // 2)
            new_start, new_end = end_exon, adjusted_end
        else:
            adjusted_start = start_exon - (min_distance // 2)
            new_start, new_end = adjusted_start, start_exon
    else:
        if strand_exon == '+':
            new_start, new_end = end_exon, end_exon + 2000
        else:
            new_start, new_end = ext_start, ext_end


    if new_start > new_end:
        new_start, new_end = new_end, new_start

    print(f"{chr_exon}\t{new_start}\t{new_end}\t{exon_name}\t0\t{strand_exon}")
