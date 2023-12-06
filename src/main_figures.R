library(tidyverse)
library(here)
library(ggtext)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(ggpmisc)

#####
here("functions", "theme_publications.R") %>% source()

read_rds(here("input", "all_atgs.rds")) -> all_atgs
read_rds(here("input", "new_atgs.rds")) -> new_atgs

####

# Fig. 1C----

read_csv(here("input","fret_new.csv"),col_names = F) %>%
  mutate(across(where(is.character), ~na_if(., "FRET E MW"))) %>% 
  mutate(type = ifelse(str_detect(X1,"^[a-zA-Z]"), X1, NA)) %>% 
  fill(type) %>% 
  mutate(X1 = X1 %>% as.numeric(),
         type = type %>% str_replace("speckles", "spec") %>% str_trim() %>% factor()) %>%
  mutate(type = type %>% factor(levels = c("PWO1-GFP+UBP5-mCh spec", "UBP5-GFP+PWO1-mCh spec",
                                           "UBP5-GFP+PWO1-mCh",
                                           "UBP5-GFP", "PWO1-GFP", "PWO1-GFP_mCh spec")) %>% fct_rev()) %>% 
  pivot_longer(cols = -type) %>% 
  select(-name) %>% 
  drop_na() -> data_1c

data_1c %>% 
  write_csv(here("output", "fig1c.csv"))

data_1c %>% 
  count(type) -> count_1c


data_1c %>%
  group_by(type) %>% 
  summarise(mean = value %>% mean,
            sd = value %>% plotrix::std.error()) -> stat_1c

aov(value ~ type, data = data_1c) %>% 
  rstatix::tukey_hsd() -> anova_1c


anova_1c %>% 
  slice(9, 11, 1, 10) -> anova_1c_m


ggplot(data_1c, aes(x = type, y = value)) +
  gg.layers::geom_boxplot2( width = 0.25, width.errorbar = 0.05, 
                            color="grey1", alpha=0.8, fill = "#EF7F31") +
  geom_jitter(data = data_1c, aes(x = type, y = value), width = 0.1, color = "blue", size = 1.5, alpha = 0.4) +
  #stat_pvalue_manual(anova_1c_m, label = "p.adj", y.position = c(32, 30, 28, 26.2), size = 6, bracket.size = 0.7) +
  stat_pvalue_manual(anova_1c_m, label = "p.adj", y.position = c(21, 16, 27, 5), size = 6, bracket.size = 0.7,
                     coord.flip = F, vjust = c(15,8, 4.5, 4.5),
                     hjust = c(-0.35, -2, -2, -0.45, -2, -2, -0.35, -2, -2, -0.5, -0.45, -0.4)
  ) +
  geom_text(data = count_1c, aes(x = type, label = paste0("n = ", n)), y = -4, size = 5, nudge_x = 0.2) +
  scale_y_continuous(breaks = c(-6, -2, 0, 2,  6, 10, 14, 18, 22, 26, 30), limits = c(-6,33)) +
  coord_flip() +
  labs(y = "FRET efficiency (%)",
       x = "",
       title = "") +
  geom_hline(yintercept = 0, alpha = 0.35)+
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(size = 15),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown()) -> p1
ggsave(here("output", 'fig1C.pdf'), units="in", width=14, height=10, dpi = 300, device = cairo_pdf)  
ggsave(here("output", 'fig1C.png'), units="in", width=14, height=10, dpi = 300)
export::graph2ppt(p1, here("output","main_figures.pptx"), width=14, height=10)

# Fig. 2B ----

here("input", "root length+ hypocotyl_Kiruba_Mohan.xlsx") %>% 
  readxl::read_excel(skip = 1) -> df

df %>% 
  #slice(-1) %>% 
  select(1:4) %>% 
  magrittr::set_colnames(c("replicate", "Col-0", "*UBP5-eGFP*", "*ubp5*")) -> root_df


here("input", "root length+ hypocotyl_Kiruba_Mohan.xlsx") %>% 
  readxl::read_excel(skip = 2) %>% 
  #slice(-(1:2)) %>% 
  select(8:11) %>% 
  magrittr::set_colnames(c("replicate", "Col-0", "*UBP5-eGFP*", "*ubp5*")) -> hypocotyl_df


root_df %>% 
  mutate(type = "Root") %>% 
  bind_rows(hypocotyl_df %>%mutate(type = "Hypocotyl")) %>% 
  pivot_longer(-c(type, replicate)) %>% 
  mutate(name = name %>% factor(levels = c("Col-0", "*UBP5-eGFP*", "*ubp5*"))) -> df_final


df_final %>% 
  filter(type == "Root") -> root_long

aov(value ~ name, data = root_long) %>% 
  rstatix::tukey_hsd() -> root_stats

summary(aov(value ~ name, data = root_long))[[1]]$`Pr(>F)`[1] -> root_anova

root_anova[[1]]

df_final %>% 
  filter(type == "Hypocotyl") -> hypocotyl_long

aov(value ~ name, data = hypocotyl_long) %>% 
  rstatix::tukey_hsd() -> hypocotyl_stats

summary(aov(value ~ name, data = hypocotyl_long))[[1]]$`Pr(>F)`[1] -> hypocotyl_anova

root_long %>%  
  ggplot(aes(x = name, y = value)) +
  geom_violin(aes(fill = name)) +
  gg.layers::geom_boxplot2(aes(fill = name), width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  stat_pvalue_manual(root_stats, label = "p.adj", y.position = c(4, 4.5, 5), size = 6, bracket.size = 0.7) +
  scale_fill_manual(values = c("#f2975a", "#94D82D", "#4472C4"))+
  scale_y_continuous(breaks = seq(0,6,1), limits = c(0,5.5)) +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(), 
        strip.text = element_markdown()) +
  labs(y = "Length (cm)",
       x = "",
       fill = "",
       title = "Root length") -> p2
ggsave(here("output", 'fig2B_a.pdf'), units="in", width=7, height=7, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig2B_a.png"),units = "in", width = 7, height = 7, dpi = 300)
export::graph2ppt(p2, here("output","main_figures.pptx"), width=7, height=7, append = T)



hypocotyl_long %>%  
  ggplot(aes(x = name, y = value)) +
  geom_violin(aes(fill = name)) +
  gg.layers::geom_boxplot2(aes(fill = name), width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  stat_pvalue_manual(hypocotyl_stats, label = "p.adj", y.position = c(1.5, 1.75, 2), size = 6, bracket.size = 0.7) +
  scale_fill_manual(values = c("#f2975a", "#94D82D", "#4472C4"))+
  scale_y_continuous(limits = c(0,2.25), breaks = seq(0,2.25, 0.5)) +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(), 
        strip.text = element_markdown()) +
  labs(y = "Length (cm)",
       x = "",
       fill = "",
       title = "Hypocotyl length") -> p3
ggsave(here("output", 'fig2B_b.pdf'), units="in", width=7, height=7, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig2B_b.png"),units = "in", width = 7, height = 7, dpi = 300)
export::graph2ppt(p3, here("output","main_figures.pptx"), width=7, height=7, append = T)  




# Fig. 2D ----

here("input", "enrichment-6.csv") %>% 
  read_csv() %>% 
  janitor::clean_names() %>% 
  mutate(`-log10(FDR)` = -log10(enrichment_fdr),
         pathway = pathway %>% fct_reorder(`-log10(FDR)`)) -> data_2d


data_2d %>% 
  ggplot(aes(y = `-log10(FDR)`, x = pathway, color = fold_enrichment)) +
  geom_segment(aes(xend = pathway, yend = 0)) +
  geom_point(aes(size = n_genes)) +
  coord_flip() +
  scale_color_continuous(low = "blue", high = "red") +
  labs(x = "",
       color = "Fold Enrichment", 
       size = "N. of Genes") +
  theme_minimal(base_family = "Arial", base_size = 12) -> p15

ggsave(here("output", 'fig2D.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig2D.png"),units = "in", width = 13, height = 8, dpi = 300)
export::graph2ppt(p15, here("output","main_figures.pptx"), width=13, height=8, append = T)

# Fig. 2C----
here("input", "seedlings.tsv") %>% 
  read_tsv() %>% 
  janitor::clean_names() -> seedling_rnaseq

read_csv(here("input", "pwo1_targets.csv")) -> pwo1_targets_df

pwo1_targets_df %>% 
  pull(gene) -> pwo1_targets

seedling_rnaseq %>% 
  mutate(pwo1_target = ifelse(gene %in% pwo1_targets, T, F),
         ubp5_target = ifelse(gene %in% all_atgs$UBP5_peaks, T, F)) %>% 
  #filter(padj < 0.05) %>% 
  mutate(sig = case_when(log2fold_change > 1 & padj < 0.05 ~ "Up",
                         log2fold_change < - 1 & padj < 0.05 ~ "Down",
                         TRUE ~ "Non-significant")) -> seedling_df
seedling_df %>% 
  write_csv(here("output", "fig2C_3A.csv"))

seedling_df %>% 
  mutate(type = case_when(log2fold_change > 1 & padj < 0.05 ~ "Up",
                          log2fold_change < - 1 & padj < 0.05 ~ "Down",
                          TRUE ~ NA)) %>% 
  mutate(sig = sig %>% factor(levels = c("Up", "Down", "Non-significant"))) %>% 
  mutate(abs_log2 = abs(log2fold_change)) %>%
  group_by(type) %>%
  arrange(type, desc(abs_log2)) %>%
  mutate(rno = row_number()) %>%
  mutate(dlabel = ifelse(rno %in% 1:10, gene, "")) %>%
  mutate(dlabel = ifelse(sig == "Non-significant", "", dlabel)) %>%
  ggplot(aes(x = log2fold_change, y = -log10(padj), color = sig, label = dlabel)) +
  geom_point() +
  geom_text_repel(max.overlaps = 20, show.legend = F) +
  scale_color_manual(values = c("#f2975a", "#4472C4","gray70")) +
  scale_x_continuous(breaks = seq(seedling_df$log2fold_change %>% min(na.rm = T) %>% floor(), 
                                  seedling_df$log2fold_change %>% max(na.rm = T) %>% ceiling(), 1)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,55,5), limits = c(0, 55)) +
  theme_Publication(base_family = "Arial") +
  labs(y = "-log10(padj)",
       x = "Log2Fold Change",
       color = "",
       title = "Seedlings RNA-seq ubp5") +
  theme(text = element_text(size = 25),
        axis.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown()) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  annotate(geom = 'table', size = 8,
           x=-6,
           y=50,
           label=list(seedling_df %>% mutate(sig = sig %>% factor(levels = c("Up", "Down", "Non-significant"))) %>% 
                        count(sig) %>% 
                        filter(sig != "Non-significant") %>%
                        rename(`log2fc > 1` = sig))) -> p4
ggsave(here("output", 'fig2C.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig2C.png"),units = "in", width = 13, height = 8, dpi = 300)
#export::graph2ppt(p4, here("output","main_figures.pptx"), width=13, height=8, append = T)
ggsave(here("output", 'fig2C.svg'), units="in", width=13, height=8, dpi = 300)

# Fig. 3A----

seedling_df %>% 
  mutate(type = case_when(log2fold_change > 1 & padj < 0.05 ~ "Up",
                          log2fold_change < - 1 & padj < 0.05 ~ "Down",
                          TRUE ~ NA)) %>% 
  mutate(nig = case_when(sig == "Non-significant" ~ "Non-significant",
                         ubp5_target == T ~ "UBP5 target",
                         TRUE ~ sig)) %>% 
  mutate(sig = sig %>% factor(levels = c("Up", "Down", "Non-significant", "UBP5 target"))) %>% 
  mutate(nig = nig %>% factor(levels = c("Up", "Down", "Non-significant", "UBP5 target"))) %>% 
  mutate(abs_log2 = abs(log2fold_change)) %>%
  group_by(type) %>%
  arrange(type, desc(abs_log2)) %>%
  mutate(rno = row_number()) %>%
  mutate(dlabel = ifelse(rno %in% 1:10, gene, "")) %>%
  mutate(dlabel = ifelse(sig == "Non-significant", "", dlabel)) %>% 
  ggplot(aes(x = log2fold_change, y = -log10(padj), color = nig, label = dlabel)) +
  geom_point() +
  geom_text_repel(max.overlaps = 20, show.legend = F) +
  #scale_color_manual(values = c("#73C08C", "#3C62E1", "gray70", "red")) +
  scale_color_manual(values = c("#f2975a", "#4472C4", "gray70", "red")) +
  scale_x_continuous(breaks = seq(seedling_df$log2fold_change %>% min(na.rm = T) %>% floor(), 
                                  seedling_df$log2fold_change %>% max(na.rm = T) %>% ceiling(), 1)) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,55,5), limits = c(0, 55)) +
  theme_Publication(base_family = "Arial") +
  labs(y = "-log10(padj)",
       x = "Log2Fold Change",
       color = "",
       title = "Seedlings RNA-seq ubp5") +
  theme(text = element_text(size = 25),
        axis.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown(),
        title = element_markdown()) +
  geom_vline(xintercept = c(1, -1), linetype = "dashed") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotate(geom = 'table', size = 8,
           x=-8,
           y=53,
           label=list(seedling_df %>% mutate(sig = sig %>% factor(levels = c("Up", "Down", "Non-significant"))) %>% 
                        count(sig, ubp5_target) %>% 
                        filter(sig != "Non-significant") %>%
                        mutate(ubp5_target = ifelse(ubp5_target == T, "Yes", "No")) %>% 
                        rename(`log2fc > 1` = sig,
                               `UBP5 targets` = ubp5_target))) -> p5
ggsave(here("output", 'fig3A.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig3A.png"),units = "in", width = 13, height = 8, dpi = 300)
#export::graph2ppt(p5, here("output","main_figures.pptx"), width=13, height=8, append = T)
ggsave(here("output", 'fig3A.svg'), units="in", width=13, height=8, dpi = 300)

# Fig. 3D----

all_atgs$UBP5_peaks %>% 
  enframe(name = NULL, value = "ubp5_targets") %>%
  mutate(cat = case_when(cat = ubp5_targets %in% new_atgs$`H2Aub-hyper` ~ "Hyper-marked",
                         ubp5_targets %in% new_atgs$`H2Aub-hypo` ~ "Hypo-marked",
                         ubp5_targets %in% new_atgs$`H2Aub-denovo` ~ "De-novo marked",
                         TRUE ~ "Unchanged")) -> data_3d
data_3d %>% 
  write_csv(here("output", "fig3d.csv"))

data_3d %>%   
  count(cat) %>% 
  ggplot(aes(fill = cat, y = n, x = "", label = n)) +
  geom_bar(position="stack", stat="identity", width = 0.6) + 
  geom_text( position = position_stack(vjust = 0.5), color = "black", size = 7) +
  scale_fill_manual(values = c("#f2975a", "#4472C4", "#94D82D", "gray"))+
  theme_Publication(base_family = "Arial") +
  theme(legend.position = "right",
        text = element_text(size = 25)) +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "",
       x = "UBP5 targets (n= 8983)",
       y = "") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 2898, 2898+669, 2898+669+1876, 2898+669+1876+3540)) -> p6
ggsave(here("output", 'fig3D.pdf'), units="in", width=8, height=10, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig3D.png"),units = "in", width = 8, height = 10, dpi = 300)
export::graph2ppt(p6, here("output","main_figures.pptx"), width=7, height=7, append = T)

# Fig. 3G----

here("input", "new_deseq2_q0.05", "H2Aub", "results.csv") %>% 
  read_csv() -> h2aub_deseq

h2aub_deseq %>%
  #filter(gene_id %in% c(all_atgs$UBP5_peaks, all_atgs$`H2Aub-Col`)) %>% 
  mutate(ip = ifelse(gene_id %in% all_atgs$UBP5_peaks, "UBP5 targets", "Non-UBP5 targets")) -> h2aub_pre

# h2aub_pre %>% 
#   count(ip) %>% 
#   mutate(n = paste0("n = ", n))-> h2aub_count

h2aub_pre %>% 
  write_csv(here("output", "fig3g.csv"))

h2aub_pre %>% 
  ggplot(aes(x = ip, y = log2FoldChange)) +
  geom_violin(aes(fill = ip)) +
  gg.layers::geom_boxplot2(aes(fill = ip), width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  #geom_text(data = h2aub_count, aes(x = ip, label = n, y = -3), size = 6) +
  stat_compare_means(label = "p.format", size = 6, label.y = 4, comparisons = list(c("Non-UBP5 targets", "UBP5 targets")
  ),bracket.size = 0.7
  ) +
  scale_y_continuous(limits = c(-4,5), breaks = seq(-4,4,2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#f2975a", "#4472C4"))+
  labs(y = "log<sub>2</sub> fold change(*ubp5*/Col)",
       x = "",
       fill = "",
       title = "H2Aub levels") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown()) -> p7

ggsave(here("output", 'fig3G.pdf'), units="in", width=9, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig3G.png"),units = "in", width = 9, height = 8, dpi = 300)
export::graph2ppt(p7, here("output","main_figures.pptx"), width=7, height=7, append = T)

# Fig. 3C----

intersect(all_atgs$`H2Aub-Col`, all_atgs$`H2Aub-ubp5`) -> intersect_h2aub

here("input", "new_deseq2_q0.05", "H2Aub", "results.csv") %>% 
  read_csv() %>% 
  filter(gene_id %in% intersect_h2aub) %>% 
  mutate(ip = case_when(gene_id %in% new_atgs$`H2Aub-hyper` ~ "Hyper-marked",
                        gene_id %in% new_atgs$`H2Aub-hypo` ~ "Hypo-marked",
                        T ~ "Unchanged") %>% factor(levels = c("Hyper-marked", "Hypo-marked", "Unchanged"))) -> h2aub_pre1

h2aub_pre1 %>% 
  write_csv(here("input", "fig3c.csv"))

h2aub_pre1 %>% 
  count(ip) %>% 
  mutate(n = paste0("n = ", n))-> h2aub_count1

h2aub_pre1 %>% 
  ggplot(aes(x = ip, y = log2FoldChange)) +
  geom_violin(aes(fill = ip)) +
  gg.layers::geom_boxplot2(aes(fill = ip), width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  geom_text(data = h2aub_count1, aes(x = ip, label = n, y = -3), size = 6) +
  stat_compare_means(label = "p.format", size = 6, label.y = c(2.5, 3, 3.5), comparisons = list(c("Hyper-marked", "Hypo-marked"),
                                                                                                c("Hypo-marked", "Unchanged"), 
                                                                                                c("Hyper-marked", "Unchanged")
  ),bracket.size = 0.7
  ) +
  #stat_compare_means(size = 6, label.y = 4.5) +
  scale_y_continuous(limits = c(-4,5), breaks = seq(-4,4,2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#f2975a", "#4472C4", "#94D82D"))+
  labs(y = "log<sub>2</sub> fold change(*ubp5*/Col)",
       x = "",
       fill = "",
       title = "H2Aub") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown()) -> p8
ggsave(here("output", 'fig3C.pdf'), units="in", width=14, height=10, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig3C.png"),units = "in", width = 14, height = 10, dpi = 300)
export::graph2ppt(p8, here("output","main_figures.pptx"), width=14, height=10, append = T)

# Fig. 4D----
here("input", "seedlings.tsv") %>% 
  read_tsv() -> rna_seq_seedlings


rna_seq_seedlings %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::filter(!(log2FoldChange < 1 & log2FoldChange > -1)) %>% 
  dplyr::select(gene, log2FoldChange) -> misregulated

here("input", "new_deseq2_q0.05", "H2Aub", "results.csv") %>% 
  read_csv() %>% 
  dplyr::filter(padj < 0.05) %>% 
  select(gene_id, log2FoldChange) -> h2aub


here("input", "new_deseq2_q0.05", "H3K27", "results.csv") %>% 
  read_csv() %>% 
  dplyr::filter(padj < 0.05) %>% 
  select(gene_id, log2FoldChange) -> h3k27

h2aub %>% 
  dplyr::rename(H2Aub = log2FoldChange) %>% 
  full_join(h3k27 %>% 
              dplyr::rename(H3K27 = log2FoldChange)) %>% 
  full_join(misregulated %>% 
              dplyr::rename(misregulated = log2FoldChange), by = c("gene_id" = "gene"))  -> all_data
all_data %>% 
  select(misregulated, H3K27) %>% 
  drop_na() %>% 
  nrow() -> h3_count

all_data %>% 
  select(misregulated, H2Aub) %>% 
  drop_na() %>% 
  nrow() -> h2a_count

scatter_plot <- function(ggplot_aes, x_name, count){
  ggplot_aes +
    geom_point() +
    #geom_smooth(method = "lm", se=FALSE) +
    scale_x_continuous(limits = c(-2, 3), breaks = seq(-2,3, 1)) +
    scale_y_continuous(limits = c(-10, 10)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    #stat_regline_equation(label.y = c(10, 6.5), aes(label = ..eq.label..), show.legend = F, size = 4.5, color = "#386CB0") +
    #stat_regline_equation(label.y = c(9, 5.5), aes(label = ..rr.label..), show.legend = F, size = 6) +
    scale_colour_Publication() +
    theme_Publication() +
    labs(x = x_name,
         y = "Transcription levels *ubp5*/WT") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_blank(), 
          axis.title.y = element_markdown(),
          plot.title = element_markdown(),
          text = element_text(size = 25)) +
    guides(label = "none") 
}

all_data %>% 
  write_csv(here("output", "fig4d_5i.csv"))

all_data %>%  
  #filter(misregulated > 0 & H2Aub < 0) %>% 
  # #filter(misregulated < 0 & H2Aub > 0) %>% 
  # #filter(misregulated < 0 & H2Aub < 0) %>% 
  # filter(misregu lated > 0 & H2Aub > 0) %>%
  ggplot(aes(x = H2Aub, y = misregulated)) %>% 
  scatter_plot("H2Aub", h2a_count) +
  geom_smooth(method = "loess", se=FALSE) +
  stat_cor(data = all_data %>% filter(misregulated > 0 & H2Aub > 0),
           aes(x = H2Aub, y = misregulated), method = "pearson",  label.x = 1, label.y = 8, size = 5, face = "bold") +
  stat_cor(data = all_data %>% filter(misregulated > 0 & H2Aub < 0),
           aes(x = H2Aub, y = misregulated), method = "pearson",  label.x = -2, label.y = 8, size = 5, face = "bold") +
  stat_cor(data = all_data %>% filter(misregulated < 0 & H2Aub < 0),
           aes(x = H2Aub, y = misregulated), method = "pearson",  label.x = -2, label.y = -7, size = 5, face = "bold") +
  stat_cor(data = all_data %>% filter(misregulated < 0 & H2Aub > 0),
           aes(x = H2Aub, y = misregulated), method = "pearson",  label.x = 1, label.y = -7, size = 5, face = "bold") +
  annotate("text",  x=3, y = 10, label = "I", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=-2, y = 10, label = "II", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=-2, y = -1, label = "III", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=3, y = -1, label = "IV", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=3, y = -9, label = paste0("n = ", h2a_count), vjust=1, hjust=1, fontface = "bold.italic", size = 7) +
  labs(y = "Transcript levels<br>(*ubp5* vs Col-0)",
       x = "H2Aub *ubp5* vs Col-0",
       title = "*ubp5* mis-regulated genes") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        plot.title = element_markdown()) -> p9
ggsave(here("output", 'fig4D.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig4D.png"),units = "in", width = 13, height = 8, dpi = 300)
export::graph2ppt(p9, here("output","main_figures.pptx"), width=13, height=8, append = T)


# Fig. 5I----

all_data %>%  
  ggplot(aes(x = H3K27, y = misregulated)) %>% 
  scatter_plot("H3K27", h3_count) +
  geom_smooth(method = "loess", se=FALSE) +
  stat_cor(data = all_data %>% filter(misregulated > 0 & H3K27 > 0),
           aes(x = H3K27, y = misregulated), method = "pearson",  label.x = 2.5, label.y = 8, size = 5, face = "bold") +
  stat_cor(data = all_data %>% filter(misregulated > 0 & H3K27 < 0),
           aes(x = H3K27, y = misregulated), method = "pearson",  label.x = -2.5, label.y = 8, size = 5, face = "bold") +
  stat_cor(data = all_data %>% filter(misregulated < 0 & H3K27 < 0),
           aes(x = H3K27, y = misregulated), method = "pearson",  label.x = -2.5, label.y = -7, size = 5, face = "bold") +
  stat_cor(data = all_data %>% filter(misregulated < 0 & H3K27 > 0),
           aes(x = H3K27, y = misregulated), method = "pearson",  label.x = 2.5, label.y = -7, size = 5, face = "bold") +
  annotate("text",  x=4.5, y = 10, label = "I", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=-2.5, y = 10, label = "II", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=-2.5, y = -1, label = "III", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  annotate("text",  x=4.5, y = -1, label = "IV", vjust=1, hjust=1, 
           fontface = "bold", size = 7) +
  scale_x_continuous(limits = c(-2.5, 4.5), breaks = seq(-2,5, 1)) +
  annotate("text",  x=4.5, y = -9, label = paste0("n = ", h3_count), vjust=1, hjust=1, fontface = "bold.italic", size = 7)+
  labs(y = "Transcript levels<br>(*ubp5* vs Col-0)",
       x = "H3K27me3 *ubp5* vs Col-0",
       title = "*ubp5* mis-regulated genes") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        plot.title = element_markdown()) -> p10
ggsave(here("output", 'fig5I.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig5I.png"),units = "in", width = 13, height = 8, dpi = 300)
export::graph2ppt(p10, here("output","main_figures.pptx"), width=13, height=8, append = T)

# Fig. 4E----

here("input", "res_althyp.txt") %>% 
  read_tsv() -> alt_hypo_result

alt_hypo_result %>% 
  filter(effect == "NO") %>% 
  pull(1) -> all_other_genes_alt

here("input", "all_genes_average_signal.csv") %>% 
  read_csv() -> all_genes_average_signal

all_genes_average_signal %>% 
  pivot_longer(-gene_id) -> avg_signal_long

avg_signal_long %>% 
  mutate(mis_reg = ifelse(gene_id %in% all_atgs$upregulated, "up", 
                          ifelse(gene_id %in% all_atgs$downregulated, "down", NA))) -> avg_signal_with_mis_reg_genes
avg_signal_with_mis_reg_genes %>% 
  mutate(mis_reg = ifelse(gene_id %in% all_other_genes_alt, "rest", mis_reg)) %>%
  mutate(type = name %>% str_remove(".*-"),
         mark = name %>% str_remove("-.*")) %>% 
  drop_na(mis_reg) -> mis_reg_mod_df_new

mis_reg_mod_df_new %>% 
  mutate(mis_reg = case_when(mis_reg == "down" ~ "Downregulated",
                             mis_reg == "up" ~ "Upregulated", 
                             mis_reg ==  "rest" ~  "Non-misegulated genes")) %>% 
  mutate(type_n = ifelse(type == "ubp5", "*ubp5*", "Col-0")) %>% 
  mutate(type_n = type_n %>% factor(levels = c("Col-0", "*ubp5*")),
         mis_reg = mis_reg %>% factor(levels = c("Downregulated", "Upregulated", "Non-misegulated genes"))) %>%
  filter(str_detect(name, "H2Aub")) -> h2aub_alt_df

h2aub_alt_df %>% 
  select(1:4) %>% 
  write_csv(here("output", "fig4e.csv"))

h2aub_alt_df %>% 
  count(mis_reg, type_n) -> h2aub_alt_count

h2aub_alt_df %>% 
  filter(!value > 20) %>% 
  ggplot(aes(x = type_n, y = value, fill = type)) +
  geom_violin() +
  gg.layers::geom_boxplot2(width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  stat_compare_means(label = "p.format", size = 6, label.y = 19, comparisons = list(c("Col-0", "*ubp5*")),bracket.size = 0.7
  ) +
  # geom_text(data = h2aub_alt_count, aes(x = type_n, label = n, y = 16), size = 6, nudge_x = 0.5) +
  scale_y_continuous(limits = c(0,21)) +
  facet_wrap(~ mis_reg, scales = "free") +
  #scale_fill_manual(values = c("#490E52", "#20688C"))+
  scale_fill_manual(values = c("#f2975a", "#4472C4")) +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown()) +
  labs(y = "Average signal gene body",
       x = "",
       fill = "",
       title = "H2Aub levels") -> p11
ggsave(here("output", 'fig4E.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig4E.png"),units = "in", width = 13, height = 8, dpi = 300)
export::graph2ppt(p11, here("output","main_figures.pptx"), width=13, height=8, append = T)

# Fig. 5H----
mis_reg_mod_df_new %>% 
  mutate(mis_reg = case_when(mis_reg == "down" ~ "Downregulated",
                             mis_reg == "up" ~ "Upregulated", 
                             mis_reg ==  "rest" ~  "Non-misegulated genes")) %>% 
  mutate(type_n = ifelse(type == "ubp5", "*ubp5*", "Col-0")) %>% 
  mutate(type_n = type_n %>% factor(levels = c("Col-0", "*ubp5*")),
         mis_reg = mis_reg %>% factor(levels = c("Downregulated", "Upregulated", "Non-misegulated genes"))) %>%
  filter(str_detect(name, "H3K27")) -> data_5h

data_5h %>% 
  select(1:4) %>% 
  write_csv(here("output", "fig5h.csv"))
data_5h %>% 
  filter(!value > 15) %>% 
  ggplot(aes(x = type_n, y = value, fill = type)) +
  geom_violin() +
  gg.layers::geom_boxplot2(width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  stat_compare_means(label = "p.format", size = 6, label.y = 15.5, comparisons = list(c("Col-0", "*ubp5*")),bracket.size = 0.7,
  ) +
  scale_y_continuous(limits = c(0,17)) +
  facet_wrap(~ mis_reg, scales = "free") +
  #scale_fill_manual(values = c("#490E52", "#20688C"))+
  scale_fill_manual(values = c("#f2975a", "#4472C4")) +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown()) +
  labs(y = "Average signal gene body",
       x = "",
       fill = "",
       title = "H3K27 levels") -> p12
ggsave(here("output", 'fig5H.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig5H.png"),units = "in", width = 13, height = 8, dpi = 300)
export::graph2ppt(p12, here("output","main_figures.pptx"), width=13, height=8, append = T)


# Fig. 5D ----

here("input", "new_deseq2", "H3K27", "results.csv")  %>% 
  read_csv() -> h3k27_deseq

h3k27_deseq %>% 
  filter(gene_id %in% c(all_atgs$UBP5_peaks, all_atgs$`H3K27-Col`)) %>% 
  mutate(ip = ifelse(gene_id %in% all_atgs$UBP5_peaks, "UBP5 targets", "Non-UBP5 targets")) -> h3k27_pre

h3k27_pre %>% 
  write_csv(here("output", "fig5d.csv"))

h3k27_pre %>% 
  count(ip) %>% 
  mutate(n = paste0("n = ", n))-> h3k27_count

h3k27_pre %>% 
  ggplot(aes(x = ip, y = log2FoldChange)) +
  geom_violin(aes(fill = ip)) +
  gg.layers::geom_boxplot2(aes(fill = ip), width = 0.1, width.errorbar = 0.025, 
                           color="grey1", alpha=0.3) +
  geom_text(data = h3k27_count, aes(x = ip, label = n, y = -3), size = 6) +
  stat_compare_means(label = "p.format", size = 6, label.y = 4, comparisons = list(c("Non-UBP5 targets", "UBP5 targets")
  ),bracket.size = 0.7
  ) +
  scale_y_continuous(limits = c(-4,5), breaks = seq(-4,4,2)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#f2975a", "#4472C4"))+
  labs(y = "log<sub>2</sub> fold change(*ubp5*/Col)",
       x = "",
       fill = "",
       title = "H3K27") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown()) -> p13

ggsave(here("output", 'fig5D.pdf'), units="in", width=9, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig5D.png"),units = "in", width = 9, height = 8, dpi = 300)
export::graph2ppt(p13, here("output","main_figures.pptx"), width=9, height=8, append = T)

# Fig. 5C----


read_plot_tsv <- function(file){
  file %>% 
    read_tsv(skip = 1) %>% 
    select(where(~!all(is.na(.x)))) %>% 
    select(-2) %>% 
    pivot_longer(-1) %>% 
    mutate(name = name %>% parse_number()) %>% 
    mutate(rep = bins %>% str_remove(".cov.*") %>% str_extract("rep*.")) %>% 
    mutate(type = bins %>% str_remove("-rep.*") %>% str_extract("-.*") %>% str_remove("-"))
}

tss_plot <- function(df, gene_count, title = "", name = "--",y_axis = ""){
  df %>% 
    ggplot(aes(x = name, y = value, color = bins)) +
    geom_line() +
    geom_line(size = 1.5) + 
    scale_x_continuous(breaks=c(1, 100, 200), labels=c("-1.0kb", "TSS", "1.0kb")) +
    labs(y = y_axis,
         x = "",
         title = title,
         #subtitle = name,
         color = "") +
    scale_colour_Publication() +
    theme_Publication(base_family = "Arial") + 
    theme(text = element_text(size = 25),
          plot.subtitle = element_text(face = "bold")) + 
    annotate("text",  x=Inf, y = Inf, label = paste0("n = ", gene_count), vjust=1, hjust=1, fontface = "bold", size = 6)
}

here("input", "others_TSS_TES_with_ubp5_gfp_ubp5_target_genes.tsv") %>% 
  read_plot_tsv() -> df_n

df_n %>% 
  write_csv(here("output", "fig5c.csv"))

df_n %>% 
  filter(str_detect(bins, "mean") & str_detect(bins, "H3K27")) %>% 
  mutate(bins = bins %>% str_extract("Col|ubp5")) %>% 
  tss_plot(8983, "UBP5-GFP-targets", "H3K27","Average enrichment") +
  scale_x_continuous(breaks=c(1, 100, 200, 300), labels=c("-1.0kb", "TSS", "TES","1.0kb"))-> p14

ggsave(here("output", 'fig5C.pdf'), units="in", width=13, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", "fig5C.png"),units = "in", width = 13, height = 8, dpi = 300)
export::graph2ppt(p14, here("output","main_figures.pptx"), width=13, height=8, append = T)

