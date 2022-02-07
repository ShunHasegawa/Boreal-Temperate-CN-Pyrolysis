# create axis labels from RDA obj for PCA plot ----------------------------

get_PCA_axislab <- function(pcaobj){
  d        <- summary(pcaobj)$cont$importance  # get eigen values
  pro_expl <- round(d[2, ] * 100, 1)
  axes     <- paste0(names(pro_expl), " (", pro_expl, "%)")
  return(axes)
}




# run RDA and generate figure ---------------------------------------------

rda_bysite <- function(x, ylimad = 1.08){
  d_sp  <- decostand(select(x, aromatic:s_lignin), method = "hellinger")
  d_rda <- rda(d_sp ~ Treatment, x)
  
  # site scores
  site_d <- data.frame(scores(d_rda, choices = 1, display = "sites", scaling = 3)) %>% 
    bind_cols(x) 
  site_p <- ggplot(site_d, aes(x = TrtID, y = RDA1))+
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_boxplot(outlier.colour = "white", size = .3)+
    geom_jitter(width = .1, alpha = .7)+
    labs(x = NULL, y = get_PCA_axislab(d_rda))+
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "none",
          strip.text.x = element_text(size = 4))
  
  
  # Sp score
  sp_d   <- data.frame(scores(d_rda, choices = 1, display = "species", scaling = 3)) %>% 
    mutate(pyr_comp = row.names(.),
           pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
           pyr_comp_ed = mapvalues(pyr_comp, 
                                   c("aromatic", "carbohydrate", "g_lignin", "Long_chain_aliphatic", "N_comp", "Others", "Phenol", "s_lignin"),
                                   c("Aromatic", "Carbohydrate", "G lignin", "Long-chain aliphatic", "N comp.", "Others", "Phenol", "S lignin"))) %>% 
    arrange(RDA1)
  
  sp_p <- ggplot(sp_d, aes(x = -1, y = RDA1, label = pyr_comp_ed)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = -1.1) +
    geom_point(aes(x = -1.1), size = 1) +
    geom_text(hjust  = 0,  size = 2) +
    labs(x = NULL, y = NULL) +
    lims(x = c(-1.15, 1.11)) +
    science_theme +
    theme(panel.border = element_blank(), 
          axis.text.x  = element_blank(), 
          axis.text.y  = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank())
  
  ylim <- c(min(c(site_d$RDA1, sp_d$RDA1)), max(c(site_d$RDA1, sp_d$RDA1))) * ylimad
  
  # merge figs
  rda_p <- ggarrange(site_p + lims(y = ylim), 
                     sp_p   + lims(y = ylim),
                     ncol = 2, widths = c(2, .9))
  
  return(list(model = d_rda, fig = rda_p))
  
}


rda_bybiome <- function(x, ylimad = 1.08){
  d_sp  <- decostand(select(x, aromatic:s_lignin), method = "hellinger")
  d_rda <- rda(d_sp ~ Treatment * Site + Condition(Site), x)
  
  # site scores
  site_d <- data.frame(scores(d_rda, choices = 1, display = "sites", scaling = 3)) %>% 
    bind_cols(x) 
  site_p <- ggplot(site_d, aes(x = TrtID, y = RDA1))+
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_boxplot(outlier.colour = "white", size = .3)+
    geom_jitter(width = .1, alpha = .7)+
    labs(x = NULL, y = get_PCA_axislab(d_rda))+
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 5),
          legend.position = "none",
          strip.text.x = element_text(size = 4))
  
  
  # Sp score
  sp_d   <- data.frame(scores(d_rda, choices = 1, display = "species", scaling = 3)) %>% 
    mutate(pyr_comp = row.names(.),
           pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
           pyr_comp_ed = mapvalues(pyr_comp, 
                                   c("aromatic", "carbohydrate", "g_lignin", "Long_chain_aliphatic", "N_comp", "Others", "Phenol", "s_lignin"),
                                   c("Aromatic", "Carbohydrate", "G lignin", "Long-chain aliphatic", "N comp.", "Others", "Phenol", "S lignin"))) %>% 
    arrange(RDA1)
  
  sp_p <- ggplot(sp_d, aes(x = -1, y = RDA1, label = pyr_comp_ed)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = -1.1) +
    geom_point(aes(x = -1.1), size = 1) +
    geom_text(hjust  = 0,  size = 2) +
    labs(x = NULL, y = NULL) +
    lims(x = c(-1.15, 1.11)) +
    science_theme +
    theme(panel.border = element_blank(), 
          axis.text.x  = element_blank(), 
          axis.text.y  = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank())
  
  ylim <- c(min(c(site_d$RDA1, sp_d$RDA1)), max(c(site_d$RDA1, sp_d$RDA1))) * ylimad
  
  # merge figs
  rda_p <- ggarrange(site_p + lims(y = ylim), 
                     sp_p   + lims(y = ylim),
                     ncol = 2, widths = c(2, .9))
  
  return(list(model = d_rda, fig = rda_p))
  
}
