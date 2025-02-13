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

dat = read.delim("../time_course_result/20241219/master_tables/20241219.aligned.v0.m1000.count_total_norm.tsv")
conditions = read.delim("../time_course_result/20241219/samples/replicates.csv", sep = ",", colClasses = c("character", "character"))

piRNA = dat %>% 
  filter(feature == "piRNA" | feature == "all_pre_miRNA") %>% 
  select(-biotype, -class) %>% 
  pivot_longer(!c(gene, feature), names_to = 'sample', values_to = 'rpm') %>% 
  mutate(sample = as.character(sample)) %>% 
  left_join(conditions, by = c("sample" = "simple_name")) %>% 
  group_by(condition, sample, feature) %>% 
  summarise(rpm = sum(rpm)) %>% 
  group_by(condition, feature) %>% 
  summarise(M = mean(rpm), S = sd(rpm)) %>% 
  ungroup() %>% 
  mutate(condition = as.character(condition)) %>% 
  mutate(genotype = ifelse(grepl("disl2", condition), "disl2", "wt")) %>% 
  mutate(time_post_depletion = as.character(gsub("disl2_", "", condition))) %>% 
  filter(time_post_depletion %in% c("TP2", "TP4", "TP6", "TP8", "TP14"))

piRNA$time_post_depletion_f = factor(piRNA$time_post_depletion, levels = c("TP2", "TP4", "TP6", "TP8", "TP14"))
piRNA$genotype_f = factor(piRNA$genotype, levels = c("wt", "disl2"))
piRNA = piRNA %>% drop_na()

p = ggplot(data = piRNA %>% filter(feature == "piRNA") %>% filter(genotype_f == "wt"), aes(x = time_post_depletion_f, y = M)) + 
  geom_bar(stat = 'identity', aes(fill = genotype_f), position = position_dodge(0.9)) + 
  my_theme() + 
  geom_errorbar(aes(x = time_post_depletion, fill = genotype_f, ymin = M - S, ymax = M + S), position = position_dodge(0.9), width = 0.3) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5), aspect.ratio = 1) + 
  scale_y_continuous(limits = c(0,60000), breaks = seq(0,60000,by = 20000)) +
  labs(y = "piRNA RPM")
p  
ggsave(p, filename = 'piRNA_depletion.pdf', dpi = 300, height = 5, width = 5)

p = ggplot(data = piRNA %>% filter(feature == "all_pre_miRNA") %>% filter(genotype_f == "wt"), aes(x = time_post_depletion_f, y = M)) + 
  geom_bar(stat = 'identity', aes(fill = genotype_f), position = position_dodge(0.9)) + 
  my_theme() + 
  geom_errorbar(aes(x = time_post_depletion, fill = genotype_f, ymin = M - S, ymax = M + S), position = position_dodge(0.9), width = 0.3) +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5), aspect.ratio = 1) + 
  labs(y = "RPM")
p  
ggsave(p, filename = 'miRNA_depletion.pdf', dpi = 300, height = 5, width = 5)


# replicates
sort(colnames(dat))
dat = dat %>% filter(feature == "piRNA")

xy_scat(dat %>% filter(feature == "piRNA"), "TOFU5_Degron_rep1_TP3_S9_R1", "TOFU5_Degron_rep2_TP3_S10_R1", log2_transform = T, fold_change = F, axmin = 2^-8, axmax = 2^10, guidelines = F)
ggsave(last_plot(), filename = "replicates_0hr.png", dpi = 300, height = 5, width = 5)

x = cor.test(dat$TOFU5_Degron_rep1_TP3_S9_R1, dat$TOFU5_Degron_rep2_TP3_S10_R1)
cor(dat$TOFU5_Degron_rep1_TP3_S9_R1, dat$TOFU5_Degron_rep2_TP3_S10_R1)

xy_scat(dat %>% filter(feature == "piRNA"), "TOFU5_Degron_rep1_TP4_S13_R1", "TOFU5_Degron_rep2_TP4_S14_R1", log2_transform = T, fold_change = F, axmin = 2^-8, axmax = 2^10, guidelines = F)
ggsave(last_plot(), filename = "replicates_2hr.png", dpi = 300, height = 5, width = 5)

xy_scat(dat %>% filter(feature == "piRNA"), "TOFU5_Degron_rep1_TP5_S17_R1", "TOFU5_Degron_rep2_TP5_S18_R1", log2_transform = T, fold_change = F, axmin = 2^-8, axmax = 2^10, guidelines = F)
ggsave(last_plot(), filename = "replicates_4hr.png", dpi = 300, height = 5, width = 5)

xy_scat(dat %>% filter(feature == "piRNA"), "TOFU5_Degron_rep1_TP6_S21_R1", "TOFU5_Degron_rep2_TP6_S22_R1", log2_transform = T, fold_change = F, axmin = 2^-8, axmax = 2^10, guidelines = F)
ggsave(last_plot(), filename = "replicates_6hr.png", dpi = 300, height = 5, width = 5)

xy_scat(dat %>% filter(feature == "piRNA"), "TOFU5_Degron_rep1_TP8_S29_R1", "TOFU5_Degron_rep2_TP8_S30_R1", log2_transform = T, fold_change = F, axmin = 2^-8, axmax = 2^10, guidelines = F)
ggsave(last_plot(), filename = "replicates_12hr.png", dpi = 300, height = 5, width = 5)


# plot RNA decay
piRNA_seqs = read_delim(file = "/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/piRNAs.txt", col_names = c("gene", "seq"), delim = "\t")

classes = read_delim("/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/piRNA_classes.tsv", col_names = c("gene", "biotype", "class"), delim = "\t") %>% distinct()

classes %>% 
  left_join(piRNA_seqs, by = "gene") %>% 
  write.table("piRNA_classes_seqs.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

plot_rna_decay <- function(data) {
  # Subset the data for the selected sequence
  #data_seq <- data[data$seq == sequence_id, ]
  for (i in 1:nrow(data)) {
    
    data_seq = data[i, ]
    # Extract parameters for the selected sequence
    A <- data_seq$A[1]
    k <- data_seq$k[1]
    seq <- data_seq$gene[1]
    uncert_A <- data_seq$uncert_A[1]
    uncert_k <- data_seq$uncert_k[1]

    # Create time points for the fitted curve
    time_fit <- seq(0,16,0.1)
    
    # Calculate fitted curve and uncertainty bounds
    fitted_curve <- A * exp(-k * time_fit) 
    fitted_upper <- (A + uncert_A) * exp(-(k - uncert_k) * time_fit) 
    fitted_lower <- (A - uncert_A) * exp(-(k + uncert_k) * time_fit) 
    
    # Create a data frame for the fit data and uncertainty bounds
    fit_data <- data.frame(
      time = time_fit,
      percent_remain = fitted_curve,
      fitted_upper = fitted_upper,
      fitted_lower = fitted_lower,
      type = "Fitted Curve"
    )
    
    # Add a column to the observed data to indicate it's empirical data
    print(data_seq)
    data_seq$type <- "Empirical Data"
    #data_seq$time = c(0,1,2,4,6,8,12,14)
    
    data_long = data_seq %>% 
      dplyr::select(-gene, -A, -half_life, -k, -Rsq, -uncert_A, -uncert_k, -starting_val, -type, -half_life_quantile, -log2_half_life, -last4, -seq, -class) %>% 
      pivot_longer(everything(), names_to = "time", values_to = "percent_remain") %>% 
      mutate(type = "Empirical Data") %>% 
      mutate(time = gsub("X", "", time))
    data_long$time = as.numeric(data_long$time)
    
    # Combine the empirical data and fit data into a single data frame for plotting
    combined_data <- rbind(
      data_long[, c("time", "percent_remain", "type")],
      fit_data[, c("time", "percent_remain", "type")]
    )
    
    combined_data['seq'] = seq
    fit_data['seq'] = seq
    
    if (i == 1) {
      combined_data_res = combined_data
      fit_data_res = fit_data
    } else {
      combined_data_res = rbind(combined_data_res, combined_data)
      fit_data_res = rbind(fit_data_res, fit_data)
    }
    
  }  
  
  # Plot the data points, fitted curve, and uncertainty
  print(combined_data_res)
  ggplot(combined_data_res %>% filter(type == "Empirical Data"), aes(x = time, y = percent_remain)) +
    geom_point(aes(color = seq), size = 2) + 
    geom_line(data = combined_data_res %>% filter(type == "Fitted Curve"), aes(x = time, y = percent_remain, color = seq)) +
    #geom_line(data = combined_data_res %>% filter(type == "Fitted Curve"), aes(color = seq, x = time, y = percent_remain)) +
    #geom_ribbon(data = fit_data_res, aes(fill = seq, x = time, ymin = fitted_lower, ymax = fitted_upper), alpha = 0.2) + 
    my_theme() + 
    theme(aspect.ratio = 1) + 
    scale_x_continuous(limits = c(0, 14), breaks = c(0,2,4,6,8,12)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.25)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8))
}


plot_rna_decay_genotypes <- function(data) {
  # Subset the data for the selected sequence
  #data_seq <- data[data$seq == sequence_id, ]
  for (i in 1:nrow(data)) {
    
    data_seq = data[i, ]
    # Extract parameters for the selected sequence
    A <- data_seq$A[1]
    k <- data_seq$k[1]
    seq <- data_seq$genotype[1]
    uncert_A <- data_seq$uncert_A[1]
    uncert_k <- data_seq$uncert_k[1]
    
    # Create time points for the fitted curve
    time_fit <- seq(0,16,0.1)
    
    # Calculate fitted curve and uncertainty bounds
    fitted_curve <- A * exp(-k * time_fit) 
    fitted_upper <- (A + uncert_A) * exp(-(k - uncert_k) * time_fit) 
    fitted_lower <- (A - uncert_A) * exp(-(k + uncert_k) * time_fit) 
    
    # Create a data frame for the fit data and uncertainty bounds
    fit_data <- data.frame(
      time = time_fit,
      percent_remain = fitted_curve,
      fitted_upper = fitted_upper,
      fitted_lower = fitted_lower,
      type = "Fitted Curve"
    )
    
    # Add a column to the observed data to indicate it's empirical data
    print(data_seq)
    data_seq$type <- "Empirical Data"
    #data_seq$time = c(0,1,2,4,6,8,12,14)
    
    data_long = data_seq %>% 
      dplyr::select(-gene, -genotype, -A, -half_life, -k, -Rsq, -uncert_A, -uncert_k, -starting_val, -type, -half_life_quantile, -log2_half_life, -last4, -seq, -class, -last3) %>% 
      pivot_longer(everything(), names_to = "time", values_to = "percent_remain") %>% 
      mutate(type = "Empirical Data") %>% 
      mutate(time = gsub("X", "", time))
    data_long$time = as.numeric(data_long$time)
    
    # Combine the empirical data and fit data into a single data frame for plotting
    combined_data <- rbind(
      data_long[, c("time", "percent_remain", "type")],
      fit_data[, c("time", "percent_remain", "type")]
    )
    
    combined_data['seq'] = seq
    fit_data['seq'] = seq
    
    if (i == 1) {
      combined_data_res = combined_data
      fit_data_res = fit_data
    } else {
      combined_data_res = rbind(combined_data_res, combined_data)
      fit_data_res = rbind(fit_data_res, fit_data)
    }
    
  }  
  
  # Plot the data points, fitted curve, and uncertainty
  print(combined_data_res)
  ggplot(combined_data_res %>% filter(type == "Empirical Data"), aes(x = time, y = percent_remain)) +
    geom_point(aes(color = seq), size = 2) + 
    geom_line(data = combined_data_res %>% filter(type == "Fitted Curve"), aes(x = time, y = percent_remain, color = seq)) +
    #geom_line(data = combined_data_res %>% filter(type == "Fitted Curve"), aes(color = seq, x = time, y = percent_remain)) +
    #geom_ribbon(data = fit_data_res, aes(fill = seq, x = time, ymin = fitted_lower, ymax = fitted_upper), alpha = 0.2) + 
    my_theme() + 
    theme(aspect.ratio = 1) + 
    scale_x_continuous(limits = c(0, 14), breaks = c(0,2,4,6,8,12)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.25)) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 8))
}

decay = read.delim("./wildtype_time_course.tsv", sep = "\t")

rsq = ggplot(data = decay, aes(x = Rsq)) + 
  geom_density(color = 'blue', fill = 'lightblue') + 
  theme_bw() + 
  ggtitle("R-squared") + 
  theme(aspect.ratio = 0.5)

halflife = ggplot(data = decay %>% filter(Rsq > 0.9), aes(x = half_life)) + 
  geom_density(color = 'blue', fill = 'lightblue') + 
  theme_bw() + 
  ggtitle("half-life") + 
  theme(aspect.ratio = 0.5)

decayrate = ggplot(data = decay %>% filter(Rsq > 0.9), aes(x = k)) + 
  geom_density(color = 'blue', fill = 'lightblue') + 
  theme_bw() + 
  ggtitle("Decay rate") + 
  theme(aspect.ratio = 0.5)

ggarrange(rsq, halflife, decayrate, ncol = 1)

# check raw staistics 
raw_decay = decay %>% 
  mutate(log2_half_life = log2(half_life)) %>% 
  mutate(half_life_quantile = as.factor(ntile(log2_half_life, 4)))

p = ggplot(data = raw_decay %>% filter(Rsq > 0.1), aes(x = half_life)) + 
  geom_histogram(aes(fill = half_life_quantile), bins = 75) + 
  scale_x_continuous(trans = 'log2') + #, limits = c(0.5,20), breaks = c(0.5, 1,2,4,8,16)) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = c("1" = "red", "2" = "mediumspringgreen", "3" = "mediumorchid", "4" = 'blue'))
p


decay_filt = decay %>% 
  filter(Rsq >= 0.80) %>%
  filter(starting_val >= 3) %>% 
  mutate(log2_half_life = log2(half_life)) %>% 
  mutate(half_life_quantile = as.factor(ntile(log2_half_life, 4))) %>% 
  left_join(piRNA_seqs, by = c("gene")) %>% 
  mutate(last4 = str_sub(seq,-5,-1)) %>% 
  mutate(last3 = str_sub(seq,-3,-1))

write.table(decay_filt, "wt_decay_processed.tsv", quote = F, col.names = T, sep = "\t", row.names = F)

decay_filt %>% top_frac(0.01, wt = k)

p = ggplot(data = decay_filt, aes(x = half_life)) + 
  geom_histogram(aes(fill = half_life_quantile), bins = 75) + 
  scale_x_continuous(trans = 'log2', limits = c(0.5,20), breaks = c(0.5, 1,2,4,8,16)) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = c("1" = "red", "2" = "mediumspringgreen", "3" = "mediumorchid", "4" = 'blue'))
p
ggsave(p, filename = "half-life-hist.pdf", dpi = 300, height = 7, width = 7, device = cairo_pdf)

nrow(decay_filt)
mean(decay_filt$half_life)
median(decay_filt$half_life)
range(decay_filt$half_life)

plot_quant = function() {

  s1 = decay_filt %>% 
    filter(half_life_quantile == 1) %>% 
    sample_n(1)
  
  s2 = decay_filt %>% 
    filter(half_life_quantile == 2) %>% 
    filter(Rsq>0.9) %>%
    sample_n(1)
  
  s3 = decay_filt %>% 
    filter(half_life_quantile == 3) %>% 
    filter(Rsq>0.9) %>% 
    sample_n(1)
  
  s4 = decay_filt %>% 
    filter(half_life_quantile == 4) %>% 
    filter(grepl("G$",last3) | grepl("C$", last3)) %>% 
    top_frac(-0.01, wt = k) %>% 
    sample_n(1)
  
  s = rbind(s1, s2, s3, s4)
  s = s %>% mutate(gene = paste0(gene," ",last3)) %>% select(-last3)
  plot_rna_decay(s)
  
}

p = plot_quant(); p
ggsave(p, filename = "regression_examples.pdf", dpi = 300, height = 7, width = 7, device = cairo_pdf)


########################
# disl-2 data
#######################
disl2_decay = read.delim("./disl2_time_course.tsv", sep = "\t")
disl2_decay_filt = disl2_decay %>% 
  filter(Rsq >= 0.80) %>% 
  filter(starting_val >= 3) %>% 
  mutate(log2_half_life = log2(half_life)) %>% 
  mutate(half_life_quantile = as.factor(ntile(log2_half_life, 4))) %>% 
  left_join(piRNA_seqs, by = c("gene")) %>% 
  mutate(last4 = str_sub(seq,-5,-1)) %>% 
  mutate(last3 = str_sub(seq,-3,-1)) 

nrow(disl2_decay_filt)
mean(disl2_decay_filt$half_life)
median(disl2_decay_filt$half_life)
range(disl2_decay_filt$half_life)

ggplot(data = disl2_decay_filt, aes(x = half_life)) + 
  geom_histogram(aes(fill = half_life_quantile), bins = 150) + 
  #scale_x_continuous(trans = 'log2') + #, limits = c(0.5,20), breaks = c(0.5, 1,2,4,8,16)) + 
  my_theme() + 
  theme(aspect.ratio = 1) + 
  scale_fill_manual(values = c("1" = "red", "2" = "mediumspringgreen", "3" = "mediumorchid", "4" = 'blue'))


disl2_hist = disl2_decay_filt %>% 
  select(half_life, gene) %>% 
  rename(dist = half_life) %>% 
  mutate(group = "disl2")

wt_hist = decay_filt %>% 
  select(half_life, gene) %>% 
  rename(dist = half_life) %>% 
  mutate(group = "wt")

hist_dat = rbind(disl2_hist, wt_hist)
hist_dat = hist_dat %>% left_join(classes, by = c("gene"))

p = ggplot(data = hist_dat, aes(x = dist)) + 
  geom_density(aes(fill = group), alpha = 0.5) + 
  scale_x_continuous(trans = 'log2') +
  my_theme() + 
  theme(aspect.ratio = 1)
p

compare_hl = wt_hist %>% 
  rename(wt_hl = dist) %>% 
  left_join(disl2_hist, by = c("gene")) %>% 
  rename(disl2_hl = dist) %>% 
  select(-group.x, -group.y) %>% 
  mutate(half_life_fc = log2(disl2_hl/wt_hl)) %>% 
  left_join(piRNA_seqs, by = c("gene")) %>% 
  left_join(classes, by = c("gene"))

# correlation between piRNA level and half-life
prgip_all = read_delim("../prg1_IP_result/20250113/master_tables/20250113.aligned.v0.m1000.count_total_norm.tsv")

prg1_IP = prgip_all %>% 
  filter(biotype == "piRNA") %>% 
  select(gene, N2prgip_17_uni) %>% 
  rename(prg1_IP = N2prgip_17_uni) %>% 
  group_by(gene) %>% 
  top_n(1, wt = prg1_IP) %>% 
  ungroup()

decay_filt_prgIP = decay_filt %>% 
  left_join(prg1_IP, by = c("gene"))

nrow(decay_filt)
nrow(decay_filt_prgIP)

decay_filt_prgIP %>% group_by(half_life_quantile) %>% count()
decay_filt %>% filter(grepl("21ur-12848", gene))


p = plot_boxplot(dat = decay_filt_prgIP, counts_col = "prg1_IP", samples_col = "half_life_quantile", ylog2 = T, dots = T, ymin = 2^-5, ymax = 2^15)
p
ggsave(p, filename = "piRNA_level_vs_half_life.png", dpi = 300, height = 5, width = 5)
ggplot(data = decay_filt_prgIP, aes(x = log2(half_life), y = log2(prg1_IP))) + 
  geom_point(size = 0.4) + 
  stat_cor()



