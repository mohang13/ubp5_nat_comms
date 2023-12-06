library(tidyverse)
library(here)
library(DESeq2)
library(tracktables)
library(GenomicRanges)
library(ggpubr)
library(ggtext)


dir(here("input", "all_three_rep", "ATGs"), full.names = T) -> atg_full_list

atg_full_list %>% 
  map(read_tsv, col_names = F) -> atg_list

names(atg_list) <- atg_full_list %>% basename() %>% str_remove_all(".ATGs")

atg_list$`H2Aub-ubp5` %>% pull() -> h2aub_up5_genes

atg_list$`H3K27-ubp5` %>% pull() -> h3k27_up5_genes

#union(intersect(atg_list$`H2Aub-ubp5` %>% pull(), atg_list$`H2Aub-Col` %>% pull()), atg_list$`H2Aub-ubp5` %>% pull()) -> h2aub_up5_genes


here("functions", "theme_publications.R") %>% source()



paste0(c(rep("H2Aub", 6), rep(c("H3K27", "H2Bub"), each = 4)), 
       rep(c("-Col-", "-ubp5-"), by = 7), 
       paste0("rep",c(rep(1:3, each = 2), rep(1:2, each = 2, times = 2)))) %>% 
  enframe(name = NULL, value = "sample") %>% 
  mutate(gtp = sample %>%                         
           str_extract("-.*") %>% 
           str_remove("-[^-]+$") %>% 
           str_remove("-") %>% as_factor()) %>% 
  mutate(chip = sample %>% str_remove("-.*") %>% 
           as_factor()) %>% 
  mutate(replicate = sample %>% 
           str_extract("-[^-]+$") %>% 
           str_remove("-rep") %>% 
           as.numeric() %>%
           as_factor()) %>% 
  mutate(file_no = c(13:18, 13:16, 13:16)) -> chip_design


here("input", "all_three_rep", "read_counts", "new_rep") %>% 
  dir(pattern = "gene_counts", full.names = T) -> file_list 


file_list %>% 
  enframe(name = NULL) %>% 
  mutate(b_name = value %>% basename()) %>% 
  mutate(chip = str_extract(b_name, "-.*") %>% str_remove(".bed") %>% str_remove("-"),
         file_no = b_name %>% str_remove("150") %>% parse_number()) -> file_names


chip_design %>% 
  left_join(file_names) -> file_name_df


file_name_df$value %>% 
  map(read_tsv, col_names = F) -> df_list

names(df_list) <- file_name_df$sample  

df_list %>% 
  imap(~ mutate(.x, name = .y)) %>% 
  map_df(bind_rows) %>% 
  pivot_wider(names_from = name, values_from = X2) %>% 
  column_to_rownames(var = "X1") -> final_count_df


# fot html table
here("input", "all_three_rep", "Araport11_GTF_genes_transposons.May2022.gtf") %>% 
  read_tsv(col_names = F) %>% 
  janitor::clean_names() -> gtf


gtf %>% 
  filter(x3 == "gene") %>% 
  mutate(gene_id = str_extract(x9, "gene_id .*") %>% str_remove("gene_id ") %>% 
           str_extract("((?![0-9]+)[A-Za-z0-9]+)")) %>% 
  select(x1, x4, x5, gene_id) %>% 
  dplyr::rename(seqnames = x1, 
                start = x4, 
                end = x5) -> gtf_new



# H2Aub----

chip_design %>% 
  filter(chip == "H2Aub") -> chip_design_h2aub


#create count matrix

final_count_df %>% 
  select(contains("H2Aub")) %>% 
  as_tibble(rownames = "genes") %>% 
  filter(genes %in% h2aub_up5_genes) %>% 
  column_to_rownames(var = "genes") %>% 
  drop_na() %>% 
  as.matrix() -> cts

############## DDS object ##################

DESeqDataSetFromMatrix(countData = cts,
                       colData = chip_design_h2aub,
                       #single variable design for comparison of ubp5 vs wt
                       design = ~ gtp) -> dds

#relevel so that Col is reference
dds$gtp <- relevel(dds$gtp, ref = "Col")

#perform differential expression analysis

# dds <- DESeq(dds)
# plotDispEsts(dds)

dds <- DESeq(dds, fitType='local')


system("mkdir -p output/all_three_rep/new_deseq2_q0.05/H2Aub")

here("output", "all_three_rep", "new_deseq2_q0.05", "H2Aub") -> save_location

png(here::here(save_location, "dispersion.png"), width = 10, height = 8, units = 'in', res = 300)
plotDispEsts(dds, main= "dispEst: local")
dev.off()

############## Regularized log transformation of Data ##################                        
rld <- rlog(dds)

DESeq2::plotPCA(rld, intgroup = c("gtp", "sample"))
ggsave(here(save_location, "PCA.png"), units="in", width=13, height=8)

res1 <- results(dds, alpha=0.1, cooksCutoff = FALSE, independentFiltering = FALSE,pAdjustMethod = "fdr")



res1 %>% 
  as_tibble(rownames = "gene_id") %>% 
  drop_na() %>%  
  arrange(pvalue) -> ordered_res1

ordered_res1 %>% 
  write_csv(here(save_location, "results.csv"))

png(here::here(save_location, "MA.png"), width = 10, height = 8, units = 'in', res = 300)
DESeq2::plotMA(res1, alpha = 0.01, ylim = c(-4, 4), main="MA plot H2Aub - Col : ubp5", ylab= "log fold change MA plot Col : ubp5") 
dev.off()

# alt hypothesis
res2 <- results(dds, lfcThreshold=0.5, altHypothesis="lessAbs", alpha=0.1, cooksCutoff = FALSE, independentFiltering = FALSE,pAdjustMethod = "fdr")


# par(mfrow=c(2,2),mar=c(2,2,1,1))
# ylim <- c(-2.5,2.5)
# 
# drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
# plotMA(res2, ylim=ylim); drawLines()

res2 %>% 
  as_tibble(rownames = "gene_id") %>% 
  drop_na() %>%  
  arrange(pvalue) -> ordered_res2

ordered_res2 %>% 
  write_csv(here(save_location, "result_althyp.csv"))

# ordered_res2 %>% #pull(gene_id) %>% length()
#   filter(padj < 0.05) %>% 
#   pull(gene_id) -> other_gene_h2aub
# 
# other_gene_h2aub %>% length()

# html table

gtf_new %>% 
  left_join(ordered_res1) %>% 
  drop_na(pvalue, padj) -> joined_gtf

joined_gtf %>% 
  filter(padj < 0.05 & log2FoldChange > 0) -> up_reg

joined_gtf %>% 
  filter(padj < 0.05 & log2FoldChange < 0) -> down_reg

myReport <- makebedtable(up_reg %>% GRanges(), "H2Aub_up.html", basedirectory = save_location)

browseURL(myReport)

myReport <- makebedtable(down_reg %>% GRanges(), "H2Aub_down.html", basedirectory = save_location)

browseURL(myReport)

# H3K27----

chip_design %>% 
  filter(chip == "H3K27") -> chip_design_h3k27

#create count matrix

final_count_df %>% 
  select(contains("H3K27")) %>% 
  as_tibble(rownames = "genes") %>% 
  filter(genes %in% h3k27_up5_genes) %>% 
  column_to_rownames(var = "genes") %>% 
  drop_na() %>% 
  as.matrix() -> cts

############## DDS object ##################

DESeqDataSetFromMatrix(countData = cts,
                       colData = chip_design_h3k27,
                       #single variable design for comparison of ubp5 vs wt
                       design = ~ gtp) -> dds

#relevel so that Col is reference
dds$gtp <- relevel(dds$gtp, ref = "Col")
#perform differential expression analysis

# dds <- DESeq(dds)
# plotDispEsts(dds)

dds <- DESeq(dds, fitType='local')


system("mkdir -p output/all_three_rep/new_deseq2_q0.05/H3K27")

here("output", "all_three_rep", "new_deseq2_q0.05", "H3K27") -> save_location

png(here::here(save_location, "dispersion.png"), width = 5, height = 4, units = 'in', res = 300)
plotDispEsts(dds, main= "dispEst: local")
dev.off()

############## Regularized log transformation of Data ##################                        
rld <- rlog(dds)

DESeq2::plotPCA(rld, intgroup = c("gtp", "sample"))
ggsave(here(save_location, "PCA.png"), units="in", width=13, height=8)

res1 <- results(dds, alpha=0.1, cooksCutoff = FALSE, independentFiltering = FALSE,pAdjustMethod = "fdr")

res1 %>% 
  as_tibble(rownames = "gene_id") %>% 
  drop_na() %>%  
  arrange(pvalue) -> ordered_res1

ordered_res1 %>% 
  write_csv(here(save_location, "results.csv"))

png(here::here(save_location, "MA.png"), width = 5, height = 4, units = 'in', res = 300)
DESeq2::plotMA(res1, alpha = 0.01, ylim = c(-4, 4), main="MA plot H3K27 - Col : ubp5", ylab= "log fold change MA plot Col : ubp5") 
dev.off()
# alt hypothesis
res2 <- results(dds, lfcThreshold=0.5, altHypothesis="lessAbs", alpha=0.1, cooksCutoff = FALSE, independentFiltering = FALSE,pAdjustMethod = "fdr")


# par(mfrow=c(2,2),mar=c(2,2,1,1))
# ylim <- c(-2.5,2.5)
# 
# drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
# plotMA(res2, ylim=ylim); drawLines()

res2 %>% 
  as_tibble(rownames = "gene_id") %>% 
  drop_na() %>%  
  arrange(pvalue) -> ordered_res2

ordered_res2 %>% 
  write_csv(here(save_location, "result_althyp.csv"))


# html table

gtf_new %>% 
  left_join(ordered_res1) %>% 
  drop_na(pvalue, padj) -> joined_gtf

joined_gtf %>% 
  filter(padj < 0.05 & log2FoldChange > 0) -> up_reg

joined_gtf %>% 
  filter(padj < 0.05 & log2FoldChange < 0) -> down_reg

myReport <- makebedtable(up_reg %>% GRanges(), "H3K27_up.html", basedirectory = save_location)

browseURL(myReport)

myReport <- makebedtable(down_reg %>% GRanges(), "H3K27_down.html", basedirectory = save_location)

browseURL(myReport)