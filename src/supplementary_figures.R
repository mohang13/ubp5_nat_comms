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


# supp fig. 1D----
readxl::read_excel(here("input","Raw data_Godwin Manuscript.xlsx"), sheet = 1, skip = 3) -> rdf1d

rdf1d %>% 
  select(1:6) %>% 
  bind_rows(rdf1d %>% 
              select(7:12)) %>% 
  fill(1) %>% 
  select(-Siliques) %>% 
  mutate(type = c(rep("Col-0", 13), rep("ubp5", 13))) %>% 
  pivot_longer(-c("type", `...1`)) %>% 
  janitor::clean_names() %>% 
  mutate(type = type %>% str_replace("ubp5", "*ubp5*") %>% factor(levels = c("*ubp5*", "Col-0")),
         name = name %>% str_replace("total" , "Total ovules") %>% 
           str_replace("unfertilized" , "Unfertilized") %>% 
           factor(levels = c("Normal", "Aborted", "Unfertilized", "Total ovules"))) %>% 
  drop_na() -> data_rdf1d 

data_rdf1d %>% 
  group_by(type, name) %>% 
  summarise(mean = value %>% mean,
            sd = value %>% sd()) -> fdf1d

ggplot(data = fdf1d, aes(x = name, y = mean, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.8), width = 0.25) +
  geom_jitter(data = data_rdf1d, aes(x = name, y = value), position = position_jitterdodge(0.5), color = "blue", show.legend = F) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,65)) +
  scale_fill_manual(values = c("gray60", "#EF7F31")) +
  labs(y = "Number of Ovules",
       x = "",
       fill = "",
       title = "*ubp5* fertilization phenotype") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        axis.text.x = element_markdown(),
        legend.text = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown(),
        plot.title = element_markdown()) 

ggsave(here("output", 'supp_fig1D.pdf'), units="in", width=14, height=10, dpi = 300, device = cairo_pdf)  
ggsave(here("output", 'supp_fig1D.png'), units="in", width=14, height=10, dpi = 300)
ggsave(here("output", 'supp_fig1D.svg'), units="in", width=14, height=10, dpi = 300)

# supp fig. 2A----

readxl::read_excel(here("input","Raw data_Godwin Manuscript.xlsx"), sheet = 2) %>% 
  janitor::clean_names() %>% 
  select(1:3) %>% 
  pivot_longer(-1) %>% #write_csv(here("output", "supp_fig2A.csv"))
  mutate(x1 = x1 %>% str_replace("Cauline leaves", "Cauline<br>leaves") %>% 
           str_replace("Rosette leaves", "Rosette<br>leaves")) %>% 
  mutate(x1 = x1 %>% factor(levels = c("Siliques", "Inflorescence", "Flowers",
                                       "Cauline<br>leaves", "Rosette<br>leaves", "Roots",
                                       "Seedlings")))-> data_2a

data_2a %>%
  group_by(x1) %>% 
  summarise(mean = value %>% mean,
            sd = value %>% sd()) -> stat_2a

aov(value ~ x1, data = data_2a) %>% 
  rstatix::tukey_hsd() %>% 
  filter(group2 == "Seedlings") -> anova_2a


ggplot() +
  geom_col(data = stat_2a, aes(x = x1, y = mean), alpha = 1, width = 0.4, fill = "#EF7F31") +
  geom_errorbar(data = stat_2a, aes(x = x1, ymin = mean - sd, ymax = mean + sd), width = 0.11) +
  geom_jitter(data = data_2a, aes(x = x1, y = value), width = 0.1, color = "blue", size = 1.5) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1, 0.1), limits = c(0,1)) +
  stat_pvalue_manual(anova_2a, label = "p.adj", y.position = c(0.95, 0.8, 0.7, 0.6, 0.5, .4), size = 6, bracket.size = 0.7) +
  labs(y = "Relative *UBP5* expression",
       x = "") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown()) 
ggsave(here("output", 'supp_fig2A.pdf'), units="in", width=15, height=10, dpi = 300, device = cairo_pdf)  
ggsave(here("output", 'supp_fig2A.png'), units="in", width=15, height=10, dpi = 300)
ggsave(here("output", 'supp_fig2A.svg'), units="in", width=14, height=10, dpi = 300)


# supp fig 3A----

readxl::read_excel(here("input","Raw data_Godwin Manuscript.xlsx"), sheet = 4, skip = 2) %>% 
  select(1:4) %>% 
  janitor::clean_names() %>% 
  pivot_longer(-1) %>% 
  mutate(x1 = x1 %>% str_replace("ubp5", "*ubp5*")) %>% 
  mutate(x1 = x1 %>% str_replace("UBP5-eGFP","*UBP5pro::UBP5-eGFP;<br>ubp5*")) %>% 
  mutate(x1 = x1 %>% factor(levels = c("Col-0", "*UBP5pro::UBP5-eGFP;<br>ubp5*", "*ubp5*"))) -> data_3a

data_3a %>% 
  group_by(x1) %>%
  summarise(mean = value %>% mean,
            sd = value %>% sd()) -> stat_3a

aov(value ~ x1, data = data_3a) %>% 
  rstatix::tukey_hsd() -> anova_3a
  
ggplot() +
  geom_col(data = stat_3a, aes(x = x1, y = mean), alpha = 1, width = 0.4, fill = "#EF7F31") +
  geom_errorbar(data = stat_3a, aes(x = x1, ymin = mean - sd, ymax = mean + sd), width = 0.11) +
  geom_jitter(data = data_3a, aes(x = x1, y = value), width = 0.1, color = "blue", size = 1.5) +
  scale_y_continuous(expand = c(0,0), breaks = seq(0,0.24, 0.03), limits = c(0,0.24)) +
  stat_pvalue_manual(anova_3a, label = "p.adj", y.position = c(0.21, 0.22, 0.20), size = 6, bracket.size = 0.7) +
  labs(y = "*UBP5/TIP41*",
       x = "") +
  theme_Publication(base_family = "Arial") +
  theme(text = element_text(size = 25),
        legend.position = "none",
        axis.text.x = element_markdown(),
        strip.text.x = element_markdown(), 
        strip.text = element_markdown(),
        axis.title.y = element_markdown())
ggsave(here("output", 'supp_fig3A.pdf'), units="in", width=9, height=8, dpi = 300, device = cairo_pdf)  
ggsave(here("output", 'supp_fig3A.png'), units="in", width=9, height=8, dpi = 300)
ggsave(here("output", 'supp_fig3A.svg'), units="in", width=9, height=8, dpi = 300)


# Supp. fig 4----

readxl::read_excel(here("input","Raw data_Godwin Manuscript.xlsx"), sheet = 6, skip = 2) %>% 
  select(1:4) %>% 
  pivot_longer(cols = starts_with("Biological"), 
               names_to = "replicate", 
               values_to = "value") %>% 
  drop_na() %>% 
  mutate(condition = rep(c("GC4", "SAMBA", "GAF1", "ACT1", "UPP"), each = 9)) %>%
  select(condition, everything()) %>% 
  mutate(value = value %>% as.numeric()) %>% 
  janitor::clean_names() %>% 
  mutate(x1 = x1 %>% str_replace("Col", "Col-0") %>% str_replace("ubp5", "*ubp5*") %>% 
           str_replace("UBP5-eGFP", "*UBP5-eGFP*") %>% 
           factor(levels = c("Col-0", "*ubp5*","*UBP5-eGFP*")),
         condition = condition %>% factor()) -> data4

data4 %>% 
  group_by(condition, x1) %>%
  summarise(mean = value %>% mean,
            sd = value %>% sd()) -> stat4


fig_4_plot <- function(data_raw, stat_data, value){
  aov(value ~ x1, data = data_raw %>% filter(condition == !!value)) %>%
   rstatix::tukey_hsd() -> anova

  stat_data %>% filter(condition == !!value) %>%
    mutate(max = mean+sd) %>%
    pull(max) %>%
    max() -> max
  
  if(max < 1){
    c(0.1, 0.2, 0.3) -> v
  } else {
    c(0.2, 0.4, 0.6) -> v
  }

  if(max < 1){
    0.4 -> l
  } else {
    0.8 -> l
  }
  
  stat_data %>% dplyr::filter(condition ==!!value) -> stat_filt
  
  data_raw %>% filter(condition == !!value) -> data_filt
  ggplot() +
    geom_col(data = stat_filt, aes(x = x1, y = mean), alpha = 1, width = 0.4, fill = "#EF7F31") +
    geom_errorbar(data = stat_filt %>% filter(condition == !!value), aes(x = x1, ymin = mean - sd, ymax = mean + sd), width = 0.11) +
    geom_jitter(data = data_filt, aes(x = x1, y = value), width = 0.1, color = "blue", size = 1.5) +
    scale_y_continuous(expand = c(0,0), 
                       #breaks = seq(0,0.24, 0.03), 
                       limits = c(0,max + l)) +
    stat_pvalue_manual(anova, label = "p.adj", y.position = c(max + v[1], max +v[3], max + v[2]), size = 6, bracket.size = 0.7) +
    labs(y = "",
         x = "",
         title = paste0("*", value, "*")) +
    theme_Publication(base_family = "Arial") +
    theme(text = element_text(size = 25),
          legend.position = "none",
          axis.text.x = element_markdown(),
          strip.text.x = element_markdown(),
          strip.text = element_markdown(),
          axis.title.y = element_markdown(),
          plot.title = element_markdown(size = 18, face = "bold"))
}

purrr::map(c("GC4", "SAMBA", "GAF1", "ACT1", "UPP"), ~ fig_4_plot(data4, stat4, .x)) -> fig4_plots

(fig4_plots[[1]] + fig4_plots[[2]] + fig4_plots[[5]]) /
  (fig4_plots[[3]]+ fig4_plots[[4]]) +
  plot_annotation(tag_levels = 'A') 
ggsave(here("output", 'supp_fig4.pdf'), units="in", width=15, height=11, dpi = 300, device = cairo_pdf)  
ggsave(here("output", 'supp_fig4.png'), units="in", width=15, height=11, dpi = 300)
ggsave(here("output", 'supp_fig4.svg'), units="in", width=15, height=11, dpi = 300)
