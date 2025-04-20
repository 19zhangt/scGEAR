
np_polt <- function(dna.fasta){
    s <- read_lines(dna.fasta)
    s <- s[!str_detect(s,"^>")]
    d <- t(str_split_fixed(s,"",nchar(s[1])))
    colnames(d) <- str_c("X",1:ncol(d))
    d <- as_tibble(d)
    d1 <- tibble(A = apply(d,1,\(x) sum(x=="A")),U = apply(d,1,\(x) sum(x=="T")),C = apply(d,1,\(x) sum(x=="C")),G=apply(d,1,\(x) sum(x=="G")))
    d <- as_tibble(d1/length(s)*100)
    d <- mutate(d,pos = 1:nchar(s[1]) - (nchar(s[1])/2 + 0.5)) %>% 
        pivot_longer(!pos,names_to = "nucleotide",values_to = "frq")
    d$nucleotide <- factor(d$nucleotide, levels = c("A","U","C", "G"))
    print(ggplot(d,aes(pos,frq,color = nucleotide))+
              geom_line() +
              theme_bw(base_size = 12) + 
              scale_color_manual(values = c("#ff6363", "#8d63ad", 
                                            "#338660", "#f0d277")) +
              xlab("Sequence position") +
              ylab("Nucleotide frequency") +
              theme(axis.text = element_text(color="black"))+
              scale_x_continuous(expand = c(0, 0)) +
              annotate("text", x=-70, y=20, label= sprintf("N=%s", length(s))))
}


test_bed <- as.data.frame(sub_clean_cs_res)
test_dir <- "/media/bora_A/zhangt/2023-09-01-sc-aQTL-Project/scGEAR/ref/PAS_sites/polyAtail//"
colnames(test_bed)[1:6] <- c("chr", "start", "end", "name", "score", "strand")

np_bed <- data.frame("chr"=test_bed$chr,
                     "start"=as.numeric(test_bed$end)-101,
                     "end"=as.numeric(test_bed$end)+100,
                     "name"=test_bed$name,
                     "score"=".",
                     "strand"=test_bed$strand)

if(!dir.exists(sprintf("%s/tmp", test_dir))){dir.create(sprintf("%s/tmp", test_dir))}
np_bed_file <- sprintf("%s/tmp/PAS_ext_ud100_%s.bed", test_dir, "findtail")
np_fasta_file <- sprintf("%s/tmp/PAS_ext_ud100_%s.fasta", test_dir, "findtail")
write.table(np_bed, file = np_bed_file, quote = F, sep = "\t", row.names = F, col.names = F)
system(command = paste0("bedtools getfasta -s -name -fi ",
                        "/media/iceland/share/Index/Genome_index/Human_hg38/GRCh38.primary_assembly.genome.fa",
                        " -bed ", np_bed_file, " -fo ", np_fasta_file))
np_polt(dna.fasta = np_fasta_file)

