library(ggpubr)
#library(tidyverse)
library(ggplot2)
#library(lemon)
library(scales)
library(RColorBrewer)
library(ggrepel)
#library(ggpmisc)
library(rstatix)
library(dplyr)
source("/fs/ess/PCON0160/ben/bin/mighty.R")

my_theme = function() {
  
  theme_classic() +
    theme(aspect.ratio = 1,
          axis.ticks = element_line(size = .5, color = "black"),
          axis.ticks.length=unit(-0.10, "cm"),
          text = element_text(size=9, color = "black", family = "Arial"),
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

volcano = function(res, x, y){
  
  res['l2fc'] = res[,paste0(x)]
  res['pv'] = res[,paste0(y)]
  
  res = res %>% 
    mutate(lab = ifelse(l2fc >= 1 & pv < 0.05, "Up",
                        ifelse(l2fc <= -1 & pv < 0.05, "Down", "Other"))) %>% 
    mutate(order = ifelse(lab == "Other", 0, 1)) %>%
    mutate(size = ifelse(lab == "Up" | lab == "Down", .8, 0.6)) %>%
    mutate(alpha = ifelse(lab == "Up" | lab == "Down", 1, 0.3))
  
  
  #cols = c("Up" = "#008000", 
  #         "Down" = "violetred1", 
  #         "Other" = "grey80")
  
  p = ggplot(data = res %>% arrange(order), aes(x = l2fc, y = -log10(pv))) + 
    geom_point(size = 1, alpha = 0.3, color = 'grey70', shape = 19) + 
    scale_size_identity() + 
    scale_alpha_identity() + 
    my_theme() + 
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.size = unit(1,"line"),
          legend.justification = c(0, 1),
          legend.position = c(0, 1),
          legend.background = element_blank(),
          legend.key=element_blank()) + 
    guides(color = guide_legend(override.aes = list(size = 2, shape = 15))) + 
    coord_cartesian() + 
   # scale_color_manual(values = cols) + 
    labs(x = paste0(x), y = paste0(y), color = "") +
    geom_hline(yintercept = 1.33333,lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_vline(xintercept = -1,lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_vline(xintercept = 1,lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)
  return(p)
}  


my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

data = read.delim("prg1_turboID.txt", sep = "\t")

data = data %>% 
  mutate(prg1_lfc = log2((prg1_1 + prg1_2)/(control_1 + control_2+0.1))) %>% 
  rowwise() %>% 
  mutate(prg1_pv = my.t.test.p.value(c(control_1, control_2), c(prg1_1, prg1_2), var.equal = T, alternative = 'less'))

data %>% filter(prg1_lfc >= 2 & prg1_pv < 0.05)

p = volcano(data, "prg1_lfc", "prg1_pv") +
  geom_point(data = data %>% filter(Alternate.ID == "cid-1"), aes(x = prg1_lfc, y = -log10(prg1_pv)), color = "magenta", size = 2) + 
  geom_point(data = data %>% filter(Alternate.ID == "prg-1"), aes(x = prg1_lfc, y = -log10(prg1_pv)), color = "blue", size = 2) + 
  geom_point(data = data %>% filter(Alternate.ID == "disl-2"), aes(x = prg1_lfc, y = -log10(prg1_pv)), color = "blue", size = 2) + 
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4,by=1)) + 
  scale_x_continuous(limits = c(-3,10), breaks = seq(-4,10,by=2))

p

ggsave(p, filename = "prg1TurboID.pdf", dpi = 300, height = 5, width = 5, device = cairo_pdf)



