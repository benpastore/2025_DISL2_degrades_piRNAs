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

colors <- c(
  "A" = "#a3f3ff",    # Light cyan for "A"
  "C" = "#ffd966",    # Light yellow for "C"
  "G" = "#92ef92",    # Light green for "G"
  "T" = "#e6b8e6",    # Light lavender for "T" (U equivalent)
  "TT" = "#d27fe3",   # Light purple for "TT"
  "other" = "#d9d9d9" # Light grey for "Other"
)


conditions = c("disl2_control", "disl2_how27", "tm5337", "old_pup1", "pup1how27", "pup1tm5337")

my_theme = function() {
  
  theme_classic() +
    theme(aspect.ratio = 1,
          axis.ticks = element_line(size = .5, color = "black"),
          axis.ticks.length=unit(-0.10, "cm"),
          text = element_text(size=9, color = "black"),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = -0.5),
          axis.text = element_text(size = 10, color = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0) ,color = "black"),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), color = "black"),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = .5, colour = "black"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.size = unit(1,"line"),
          legend.background = element_blank(),
          legend.key=element_blank()
    )
}

plot_tails = function(tailor_file, conditions, outname) {

  tailor = read_delim(tailor_file)
  
  my_dat = tailor %>% filter(condition %in% conditions)
  my_dat
  my_dat$condition_f = factor(my_dat$condition, levels = conditions)
  my_dat$tail = factor(my_dat$tail_group, levels = c("A", "C", "G", "T", "TT", "other"))
  p1 = ggplot(data = my_dat %>% filter(tail != "other"), aes(x = condition_f, y = mean_tailing_percent, fill = tail)) + 
    geom_bar(stat = 'identity', width = 0.7) + 
    my_theme() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5), aspect.ratio = 1) + 
    scale_fill_manual(values = colors)  
    #scale_y_continuous(limits = c(0,3), breaks = seq(0,3.5,0.5))
  print(p1)
  
  ggsave(p1, filename = paste0(outname,"_tailing.pdf"), dpi = 300, height = 5, width = 5)

}

plot_tails("../piRNA_tailing.tsv" , conditions, "piRNA")
plot_tails("../miRNA_tailing.tsv" , conditions, "miRNA"); ggsave(last_plot(), filename = "miRNA_tailing.pdf", dpi = 300, height = 5, width = 5)
plot_tails("../wago_tailing.tsv" , conditions, "wago"); ggsave(last_plot(), filename = "wago_tailing.pdf", dpi = 300, height = 5, width = 5)
plot_tails("../csr_tailing.tsv" , conditions, "csr"); ggsave(last_plot(), filename = "csr_tailing.pdf", dpi = 300, height = 5, width = 5)
read.delim("../piRNA_tailing.tsv") %>% filter(condition == "disl2_control")

