library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(scales)
library(ggpubr)
library(tidyr)
library(readr)
library(rstatix)
library(gghalves)
library(ggbeeswarm)
library(stringr)
source("/fs/ess/PCON0160/ben/bin/mighty.R")

piRNA_seqs = read.delim("/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/piRNAs.txt", colClasses = c("character", "character"), col.names = c("gene", "seq"))

disl2 <- read_delim("/fs/ess/PCON0160/ben/projects/2024_disl2/smRNA_seq_mutants_remove_rRNA/dge/disl2_control_vs_disl2_how27.total_norm.tsv")
p = xy_dge(disl2 %>% filter(feature == "piRNA"), "disl2_control", "disl2_how27", log2_transform = T, MA = F, axmin = 2^-8, axmax = 2^9, fc_max = 6, fc_min = -6, fold_change = 2)
p
ggsave(p, filename = "disl2_piRNA.png", dpi = 300, height = 5, width = 5)

disl2 %>% 
  filter(feature == "piRNA") %>% 
  left_join(piRNA_seqs, by = c("gene")) %>% 
  write.table(., "disl2_vs_control_piRNA.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

pup1 <- read_delim("/fs/ess/PCON0160/ben/projects/2024_disl2/smRNA_seq_mutants_remove_rRNA/dge/old_N2_1_vs_old_pup1.total_norm.tsv")
p = xy_dge(pup1 %>% filter(feature == "piRNA"), "old_N2_1", "old_pup1", log2_transform = T, axmin = 2^-8, axmax = 2^9, fc_max = 6, fc_min = -6, fold_change = 2)
p
ggsave(p, filename = "pup1_piRNA.png", dpi = 300, height = 5, width = 5)

