
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


source1 <- read.table(sprintf("%s/refSeq/refseq_tts_in_last_exon.bed", argvs[1]), sep = "\t")
use_source1 <- source1[, 1:6] %>% mutate(Name=paste0(V1, ":", V3, ":", V6), score="RefSeq") %>%  select(chr=V1, start=V2, end=V3, Name, score, strand=V6) %>% unique()

source2 <- read.table(sprintf("%s/polyA_DB/UCSC_polyDB_in_last_exon.bed", argvs[1]), sep = "\t")
use_source2 <- source2[, 1:6] %>% mutate(Name=paste0(V1, ":", V3, ":", V6), score="polyA_DB") %>%  select(chr=V1, start=V2, end=V3, Name, score, strand=V6) %>% unique()

source3 <- read.table(sprintf("%s/mws_cs_annotation/mws_cs_info.bed", argvs[1]), sep = "\t")
use_source3 <- source3 %>% mutate(Name=paste0(V1, ":", V3, ":", V6), score="MWS") %>%  select(chr=V1, start=V2, end=V3, Name, score, strand=V6) %>% unique()

source4 <- read.table(sprintf("%s/findpolyAtail/findtail_PASs.bed", argvs[1]), sep = "\t")
use_source4 <- source4[, 1:6] %>% filter(!is.na(V5)) %>% mutate(Name=paste0(V1, ":", V3, ":", V6), score="scRNA_seq") %>%  select(chr=V1, start=V2, end=V3, Name, score, strand=V6) %>% unique()

source5 <- read.table(sprintf("%s/../Annotation/KeepPAS_data_annotation.txt", argvs[1]), sep = "\t", header = T)
use_source5 <- source5 %>% mutate(start=site_position-1, score="Integration", Name=paste0(chr, ":", site_position, ":", strand)) %>%  
    select(chr, start, end=site_position, Name, score, strand) %>% unique()


library(bedtoolsr)
library(pheatmap)

check_data <- use_source5
source1_support <- bt.closest(a = check_data %>% arrange(chr,start), b = use_source1 %>% arrange(chr,start), d = T, s = T)
source1_support <- source1_support$V4[source1_support$V13 <= 40] %>% unique()
source2_support <- bt.closest(a = check_data %>% arrange(chr,start), b = use_source2 %>% arrange(chr,start), d = T, s = T)
source2_support <- source2_support$V4[source2_support$V13 <= 40] %>% unique()
source3_support <- bt.closest(a = check_data %>% arrange(chr,start), b = use_source3 %>% arrange(chr,start), d = T, s = T)
source3_support <- source3_support$V4[source3_support$V13 <= 40] %>% unique()
source4_support <- bt.closest(a = check_data %>% arrange(chr,start), b = use_source4 %>% arrange(chr,start), d = T, s = T)
source4_support <- source4_support$V4[source4_support$V13 <= 40] %>% unique()

tmp_df <- as.data.frame(check_data$Name)
rownames(tmp_df) <- tmp_df$`check_data$Name`
tmp_df$polyAtail <- tmp_df$MWS <- tmp_df$polyDB <- tmp_df$RefSeq <- 0
tmp_df$RefSeq[rownames(tmp_df)%in%source1_support] <- 1
tmp_df$polyDB[rownames(tmp_df)%in%source2_support] <- 1
tmp_df$MWS[rownames(tmp_df)%in%source3_support] <- 1
tmp_df$polyAtail[rownames(tmp_df)%in%source4_support] <- 1
tmp_df <- tmp_df[, -1]
rownames(tmp_df)[apply(tmp_df, 1, sum)==0][1:10]
print(apply(tmp_df, 2, sum)/nrow(tmp_df)*100)
tmp_order_df <- tmp_df %>% arrange(desc(RefSeq), desc(polyDB), desc(MWS), desc(polyAtail)) %>% t()
dim(tmp_order_df)
pheatmap(t(tmp_order_df), cluster_rows = F, cluster_cols = F, breaks = c(0,0.5,1), 
         show_rownames = F, treeheight_col = 15, height = 3, width = 5, 
         border_color = NA, 
         color = colorRampPalette(colors = c("#f9fafb", "#e23b1b"))(2))

polyAtail_PAS <- colnames(tmp_order_df)[apply(tmp_order_df, 2, function(x){x[4]==sum(x)})]
polyAtail_PAS <- lapply(polyAtail_PAS, function(x){
    x_str <- strsplit(x, split = ":")[[1]]
    return(c(x_str[1], as.numeric(x_str[2])-1, as.numeric(x_str[2]), x, ".", x_str[3]))
})
polyAtail_PAS <- do.call('rbind', polyAtail_PAS) %>% as.data.frame()


test_bed <- as.data.frame(polyAtail_PAS)
colnames(test_bed)[1:6] <- c("chr", "start", "end", "name", "score", "strand")

np_bed <- data.frame("chr"=test_bed$chr,
                     "start"=as.numeric(test_bed$end)-101,
                     "end"=as.numeric(test_bed$end)+100,
                     "name"=test_bed$name,
                     "score"=".",
                     "strand"=test_bed$strand)

test_dir <- "example/APA/PAS_source/findpolyAtail/" #argvs[1]
if(!dir.exists(sprintf("%s/tmp", test_dir))){dir.create(sprintf("%s/tmp", test_dir))}
np_bed_file <- sprintf("%s/tmp/PAS_ext_ud100_%s.bed", test_dir, "findtail")
np_fasta_file <- sprintf("%s/tmp/PAS_ext_ud100_%s.fasta", test_dir, "findtail")
write.table(np_bed, file = np_bed_file, quote = F, sep = "\t", row.names = F, col.names = F)
system(command = paste0("bedtools getfasta -s -name -fi ",
                        "/media/iceland/share/Index/Genome_index/Human_hg38/GRCh38.primary_assembly.genome.fa",
                        " -bed ", np_bed_file, " -fo ", np_fasta_file))
np_polt(dna.fasta = np_fasta_file)
