
############# PACKAGE TEST ############# 
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE,repos='http://cran.us.r-project.org')
    if(!require(x,character.only = TRUE)) stop ("Failed to install the package. Please check the internet access or update your R if it is too old.")
  }
}


########### colors ################
orange = "#e69f00"
skyblue = "#56b4e9"
bluegreen = "#009273"
yellow = "#f0e442"
blue = "#0072b2"
red = "#d55e00"
purple = "#cc79a7"

############# THEME ############# 
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

my_theme_free_aspect = function() {
    
    theme_classic() +
    theme(
          axis.ticks = element_line(size = .5, color = "black"),
          axis.ticks.length=unit(-0.13, "cm"),
          text = element_text(size=10, color = "black"),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = -0.5),
          axis.text = element_text(size = 10, color = "black"),
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0) ,color = "black"),
          axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), color = "black"),
          panel.grid = element_blank(),
          axis.line = element_line(size = .5, colour = "black")
    ) 
}

############# GLOBAL PARAMETERS ############# 
BIGPOINT=2
MEDPOINT=1
SMALLPOINT=0.6
ALPHAPOINT=0.7
LINETYPE="dashed"
LINECOLOR="grey60"
LINEALPHA = 1
LINEWIDTH = 0.4
FC=2
PVALUE=0.05
RPM = 0


############# Helper Functions ###############
z_score = function(N, x, s){
  
  z = (N - x) / s
  return(z)
  
}



############# HELPER FUNCTIONS ############# 
pick_pt_size = function(res){
    
    if (nrow(res)<1000){
    point_size = BIGPOINT
    } else if (nrow(res) < 3000){
    point_size = MEDPOINT
    } else {
    point_size = SMALLPOINT
    }
    
    return(point_size)
}

#############   Volcano    ############# 
volcano = function(res, x, y){
  
  res['l2fc'] = res[,paste0(x)]
  res['pv'] = res[,paste0(y)]
  
  res = res %>% 
    mutate(lab = ifelse(l2fc >= 1 & pv < 0.05, "Up",
                        ifelse(l2fc <= -1 & pv < 0.05, "Down", "Other"))) %>% 
    mutate(order = ifelse(lab == "Other", 0, 1)) %>%
    mutate(size = ifelse(lab == "Up" | lab == "Down", .8, 0.6)) %>%
    mutate(alpha = ifelse(lab == "Up" | lab == "Down", 1, 0.3))
    
  
  cols = c("Up" = "springgreen3", 
           "Down" = "violetred1", 
           "Other" = "grey80")
  
  p = ggplot(data = res %>% arrange(order), aes(x = l2fc, y = -log10(pv), color = lab)) + 
    geom_point(aes(size = size), alpha = 1) + 
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
    scale_color_manual(values = cols) + 
    labs(x = paste0(x), y = paste0(y), color = "") +
    geom_hline(yintercept = 1.33333,lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_vline(xintercept = -1,lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_vline(xintercept = 1,lty = LINETYPE, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)
  return(p)
}  


############# X vs Y plots ############# 
xy_plot = function(res, x, y, cols, labs, rho, m, M, log2_transform = F, log10_transform = F, cor_method = "Pearson's"){
  
    print(labs)
  
    print(cols)
  
    if (is.null(labs)) { 
      res$lab = "all"
      cols = c("all" = "grey50")
      labs = c(paste0("all ", nrow(res)))
    }
  
    p = ggplot(data = res %>% arrange(order), aes(x = X, y = Y)) + 
    #geom_point(aes(alpha = alpha, size = point_size, color = lab), shape = 16) + 
    geom_point(shape = 19, aes(alpha = alpha, color = lab, size = point_size)) + 
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
      guides(
        color = guide_legend(
          override.aes = list(
            size = 3,  # Adjust size for legend swatches
            shape = 15  # Box shape
          )
        )
      ) + 
    coord_cartesian() + 
    scale_color_manual(values = cols, labels = labs, drop = F) + 
    #scale_fill_manual(values = cols, labels = labs) + 
    labs(x = paste0(x), y = paste0(y), color = "") 
    
    if (log2_transform == TRUE) {
      p = p + 
        scale_y_continuous(
        limits = c(m, M),
        trans = "log2",
        #labels = trans_format("log2", math_format(2^.x)),
        labels = trans_format("log2", math_format(.x)),
        breaks = trans_breaks("log2", function(x) 2^x)
      ) + 
      scale_x_continuous(
        limits = c(m, M),
        trans = "log2",
        #labels = trans_format("log2", math_format(2^.x)),
        labels = trans_format("log2", math_format(.x)),
        breaks = trans_breaks("log2", function(x) 2^x)
      )
    }
    
    if (log10_transform == TRUE) {
      p = p + 
        scale_y_continuous(
          limits = c(m, M),
          trans = "log10",
          labels = trans_format("log10", math_format(.x)),
          breaks = trans_breaks("log10", function(x) 10^x)
        ) + 
        scale_x_continuous(
          limits = c(m, M),
          trans = "log10",
          labels = trans_format("log10", math_format(.x)),
          breaks = trans_breaks("log10", function(x) 10^x)
        ) + 
        annotation_logticks(sides = 'bl', base = 10)
    }
    
    
    
    if (rho != "None"){
        p = p + labs(subtitle = paste0(cor_method, " rho: ", round(rho, 5)))
    } 
    
    return(p)
}

xy_dge = function(res, x, y, axmin = 2^-10, axmax = 2^15, deseq = FALSE, fold_change = FALSE, pvalue = FALSE, log2_transform = FALSE, log10_transform = FALSE, min_abund = FALSE, MA = FALSE, fc_min = -10, fc_max = 10, cor_method = 'pearson') {
    
    if (log2_transform == T & log10_transform == T) {
      print("Cannot do both log2 and log10 transformation defaulting to log2....")
      log10_transform = F
    }
    point_size = pick_pt_size(res)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    
    cor = cor.test(res$X, res$Y, method = cor_method)
    rho = cor(res$X, res$Y, method = cor_method)
    
    if (deseq) {
      res = res %>% select(X, Y, log2FoldChange, padj, baseMean)
      res = res %>% dplyr::rename(pvalue = padj, basemean = baseMean) 
    } else {
      if (MA == TRUE) {
        res = res %>% select(X, Y, log2FoldChange, pvalue, basemean)
      } else {
        res = res %>% select(X, Y, log2FoldChange, pvalue, basemean)
      }
    }
    
    res = res  %>% filter(X > 0 | Y > 0)

    if (fold_change == FALSE){
      fc = FC
    } else {
      fc = fold_change
    }
    
    if (pvalue == FALSE){
      pv = PVALUE
    } else {
      pv = pvalue
    }
  
    if (min_abund != FALSE) {
      res = res %>% 
        mutate(lab = ifelse(log2FoldChange >= log2(fc) & pvalue < pv & Y >= min_abund, "Up",
                            ifelse(log2FoldChange <= log2(1/fc) & pvalue < pv & X >= min_abund, "Down", "None"))) 
    } else {
      res = res %>% 
        mutate(lab = ifelse(log2FoldChange >= log2(fc) & pvalue < pv & X >= min_abund, "Up",
                            ifelse(log2FoldChange <= log2(1/fc) & pvalue < pv & Y >= min_abund , "Down", "None"))) 
    }
    
    res = res %>% 
      mutate(lab = ifelse(is.na(lab), "None", lab)) %>% 
      mutate(color = ifelse(lab == "Up", "blue", 
                            ifelse(lab == "Down", "magenta", "grey70"))) %>% 
      mutate(order = ifelse(lab == "None", 0, 1)) %>%
      mutate(point_size = ifelse(lab == "Up" | lab == "Down", point_size, point_size/2)) %>%
      mutate(alpha = ifelse(lab == "Up" | lab == "Down", 0.5, 0.4))
    
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    up = res %>% filter(lab == "Up") 
    down = res %>% filter(lab == "Down") 
    unchanged = res %>% filter(lab == "None")
    
   # cols = c("Up" = "blue", 
  #         "Down" = "magenta", 
  #         "None" = "grey70")
    
  #  labs = c(
  #    "blue" =  paste0("Fold Change >= ",fc," & p < ",pv," (n =",nrow(up),")"),
  #    "magenta" = paste0("Fold Change <= ", round(1/fc, 2)," & p < ",pv," (n =",nrow(down),")"),
  #    "grey70" = "Other")
    
    
    cols = c("Up" = "blue",
             "Down" = "magenta",
             "None" = "grey70")
    
    res$lab <- factor(res$lab,levels = c("Up", "Down", "None"))
    
    labs = c(paste0("Fold Change >= ",fc," & p < ",pv," (n =",nrow(up),")"),
             paste0("Fold Change <= ", round(1/fc, 2)," & p < ",pv," (n =",nrow(down),")"),
             "Other")

    
    dummy_data <- data.frame(
      X = c(NA, NA, NA),
      Y = c(NA, NA, NA),
      log2FoldChange = c(NA, NA, NA),
      pvalue = c(NA, NA, NA),
      basemean = c(NA, NA, NA),
      lab = c("Up", "Down", "None"),
      color = c("blue", "magenta", "grey70"),
      order = c(1, 1, 0),  # Maintain order for aesthetics
      point_size = c(1, 1, 1),  # Default size for legend points
      alpha = c(1, 1, 1)  # Default alpha for legend points
    )
    
    print(res %>% head(3))
    print(dummy_data)
    
    # Combine dummy data with the real dataset
    res <- rbind(res, dummy_data)
    
    if (MA == TRUE) {
      p = ggplot(data = res %>% arrange(order), aes(x = basemean, y = log2FoldChange, color = lab, size = point_size, alpha = alpha)) + 
        geom_point() + 
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
        scale_color_manual(values = cols, labels = labs) + 
        labs(x = paste0("Basemean"), y = paste0("log2(FoldChange)"), color = "", title = paste0(x, " vs ", y)) + 
        scale_x_continuous(
          limits = c(axmin, axmax),
          trans = "log2",
          labels = trans_format("log2", math_format(2^.x)),
          breaks = trans_breaks("log2", function(x) 2^x)
        ) + 
        scale_y_continuous(
          limits = c(fc_min, fc_max)
        ) + 
        labs(subtitle = paste0(cor_method,"'s rho: ", round(rho, 2)))
      
      p = p +
        geom_hline(yintercept = 0, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
        geom_hline(yintercept = -1, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
        geom_hline(yintercept = 1, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
    } else {
      #p = ggplot(data = res %>% arrange(order), aes(x = X, y = Y, color = color)) + 
      #  geom_point(aes(size = point_size, alpha = alpha)) + 
      #  scale_size_identity() + 
      #  scale_alpha_identity() + 
      #  scale_color_identity(guide = guide_legend(), label = labs) + 
      #  my_theme() + 
      #  theme(legend.title = element_text(size = 8),
      #        legend.text = element_text(size = 8),
      #        legend.key.size = unit(1,"line"),
      #        legend.justification = c(0, 1),
      #        legend.position = c(0, 1),
      #        legend.background = element_blank(),
      #        legend.key = element_blank()) + 
      #  coord_cartesian() + 
      #  labs(x = paste0(x), y = paste0(y), color = "") + 
      #  geom_abline(slope = 1, intercept = 0, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
      p = xy_plot(res, x, y, cols, labs, rho, axmin, axmax, log2_transform, log10_transform) + 
        geom_abline(slope = 1, intercept = 0, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)
      if (log2_transform == TRUE) {
        p = p + 
          geom_abline(slope = 1, intercept = log2(fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
          geom_abline(slope = 1, intercept = log2(1/fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
      } else if (log10_transform == TRUE){
        p = p + 
          geom_abline(slope = 1, intercept = log10(fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
          geom_abline(slope = 1, intercept = log10(1/fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
      } else {
        p = p + 
          geom_abline(slope = 1, intercept = fc, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
          geom_abline(slope = 1, intercept = 1/fc, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)  
      }
    }
    
    if (min_abund != F) {
      p = p + 
        #geom_segment(aes(x = 0, xend = min_abund, y = min_abund, yend = min_abund), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, size = 0.1) + 
        #geom_segment(aes(y = 0, yend = min_abund, x = min_abund, xend = min_abund), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, size = 0.1)
        geom_hline(yintercept = min_abund, lty = "dashed", color = "grey80", alpha = LINEALPHA, lwd = LINEWIDTH) +
        geom_vline(xintercept = min_abund, lty = "dashed", color = "grey80", alpha = LINEALPHA, lwd = LINEWIDTH) 
      
      
    }

    return(p)
}

MA_dge = function(res, x, y, fc_min = -10, fc_max = 10, axmin = 2^-10, axmax = 2^15, deseq = FALSE, fold_change = FALSE, pvalue = FALSE, min_abund = FALSE) {
  
  point_size = pick_pt_size(res)
  
  res['X'] = res[ , paste0(x)]
  res['Y'] = res[ , paste0(y)]
  
  cor = cor.test(res$X, res$Y, method = "pearson")
  rho = cor(res$X, res$Y, method = "pearson")
  
  if (deseq) {
    res = res %>% select(X, Y, basemean, log2FoldChange, padj)
    res = res %>% dplyr::rename(pvalue = padj)
  } else {
    res = res %>% select(X, Y, basemean, log2FoldChange, pvalue)
  }
  
  if (fold_change == FALSE){
    fc = FC
  } else {
    fc = fold_change
  }
  
  if (pvalue == FALSE){
    pv = PVALUE
  } else {
    pv = pvalue
  }
  
  if (min_abund != FALSE) {
    res = res %>% 
      mutate(lab = ifelse(log2FoldChange >= fc & pvalue < pv & X >= min_abund, "Up",
                          ifelse(log2FoldChange <= -1*fc & pvalue < pv & Y >= min_abund , "Down", "None"))) 
  } else {
    res = res %>% 
      mutate(lab = ifelse(log2FoldChange >= fc & pvalue < pv & Y >= RPM, "Up",
                          ifelse(log2FoldChange <= -1*fc & pvalue < pv & X >= RPM , "Down", "None"))) 
  }
  
  res = res %>% 
    mutate(lab = ifelse(is.na(lab), "None", lab)) %>% 
    mutate(order = ifelse(lab == "None", 0, 1)) %>%
    mutate(point_size = ifelse(lab == "Up" | lab == "Down", 2, 0.70)) %>%
    mutate(alpha = ifelse(lab == "Up" | lab == "Down", 0.5, 0.5))
  
  res[which(res$X == 0), 'X'] = xmin
  res[which(res$Y == 0), 'Y'] = xmin
  
  up = res %>% filter(lab == "Up") 
  down = res %>% filter(lab == "Down") 
  unchanged = res %>% filter(lab == "None")
  
  cols = c("Up" = "#008000", 
           "Down" = "red2", 
           "None" = "grey80")
  
  labs = c(
    paste0("log2FC >= ",1*fc," & p < ",pv," (n =",nrow(up),")"),
    paste0("log2FC <= ",-1*fc," & p < ",pv," (n =",nrow(down),")"),
    "Other")
  
  p = ggplot(data = res %>% arrange(order), aes(x = basemean, y = log2FoldChange, color = lab, size = point_size, alpha = alpha)) + 
      geom_point() + 
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
      scale_color_manual(values = cols, labels = labs) + 
      labs(x = paste0("Basemean"), y = paste0("log2(FoldChange)"), color = "", title = paste0(x, " vs ", y)) + 
      scale_x_continuous(
        limits = c(xmin, xmax),
        trans = "log2",
        labels = trans_format("log2", math_format(2^.x)),
        breaks = trans_breaks("log2", function(x) 2^x)
      ) + 
      scale_y_continuous(
        limits = c(ymin, ymax)
      ) 
  
  p = p +
    geom_hline(yintercept = 0, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_hline(yintercept = -1, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
    geom_hline(yintercept = 1, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
  
  return(p)
}


xy_scat = function(res, x, y, axmin = 2^-10, axmax = 2^15, log2_transform = FALSE, log10_transform = FALSE, guidelines = TRUE, fold_change = FALSE, cor_method = 'pearson'){
  
    if (log2_transform == T & log10_transform == T) {
      print("Cannot do both log2 and log10 transformation defaulting to log2....")
      log10_transform = F
    }
    
    point_size = pick_pt_size(res)
    
    if (fold_change == FALSE){
      fc = FC
    } else {
      fc = fold_change
    }
    
    print(fc)
    
    res['X'] = res[ , paste0(x)]
    res['Y'] = res[ , paste0(y)]
    res = res %>% select(X, Y) %>% replace_na(list(X = 0, Y = 0))
    
    cor = cor.test(res$X, res$Y, method = cor_method)
    rho = cor(res$X, res$Y, method = cor_method)
    #rho = 1
    res['xy'] = res$Y / res$X
    res[which(res$X == 0), 'X'] = axmin
    res[which(res$Y == 0), 'Y'] = axmin
    
    res = res %>% 
      mutate(lab = ifelse(xy >= fc, "up", 
                          ifelse(xy <= (1/fc), "down", "other")))
    
    print(res %>% head())
    
    
    #res['lab'] = "All"
    
    #cols = c("up" = "#008000", 
    #         "down" = "red2", 
    #         "other" = "grey80")
    cols = c("up" = "black", "down" = "black", "other" = "black")
    
    res['order'] = 1

    #cols = c("All" = "black")
    #labs = c(paste0("All N = ",nrow(res)))
    
    #labs = c(
    #  paste0("log2FC >= ",1*fc, " ", nrow(res %>% filter(lab == 'up'))),
    #  paste0("log2FC <= ",-1*fc," ", nrow(res %>% filter(lab == 'down'))),
    #  "Other")
    
    res['point_size'] = point_size
    res['alpha'] = 0.5
    
    p = xy_plot(res, x, y, cols, labs = NULL, rho, axmin, axmax, log2_transform, log10_transform, cor_method = cor_method)
    
    if (guidelines == T){
      p = p + 
        geom_abline(slope = 1, intercept = 0, color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
      
      if (fold_change != FALSE) {
        if (log2_transform == TRUE) {
          p = p + 
            geom_abline(slope = 1, intercept = log2(fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
            geom_abline(slope = 1, intercept = log2(1/fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
        } else if (log10_transform == TRUE){
          p = p + 
            geom_abline(slope = 1, intercept = log10(fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
            geom_abline(slope = 1, intercept = log10(1/fc), lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) 
        } else {
          p = p + 
            geom_abline(slope = 1, intercept = fc, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH) + 
            geom_abline(slope = 1, intercept = -1*fc, lty = "dashed", color = LINECOLOR, alpha = LINEALPHA, lwd = LINEWIDTH)  
        }
      }
    }
    
    return(p)
    
}
