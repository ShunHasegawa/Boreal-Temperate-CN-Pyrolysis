


spect_litter_prop_alcntr <- ddply(spect_litter_prop, .(Site), function(x){
  # Fertilised
  ddf <- filter(x, treatment == "Fertilised") %>% 
    mutate(Site_Trt = paste(Site, Trt, sep = "_"))
  
  # Fertilised treatment for each site
  st_val <- unique(ddf$Site_Trt)
  
  # Replicate a data matric of control for each treatment
  ddc <- filter(x, treatment == "Control")
  dd_new <- ldply(st_val, function(y) mutate(ddc, Site_Trt = y)) %>% 
    bind_rows(ddf)
  
  return(dd_new)
})
names(spect_litter_prop_alcntr)
litter_pyr_sp <- decostand(select(spect_litter_prop_alcntr, aromatic:s_lignin), method = "hellinger")
litter_pyr_rda <- rda(litter_pyr_sp ~ Trt * Site_Trt + Condition(Site_Trt), spect_litter_prop_alcntr)

anova(litter_pyr_rda, strata = spect_litter_prop_alcntr$Site_Trt)
summary(litter_pyr_rda)
plot(litter_pyr_rda)
ordisurf(litter_pyr_rda, log(spect_litter_prop_alcntr$wC))
litter_rda_pyr_df <- data.frame(scores(litter_pyr_rda)$sites) %>% 
  bind_cols(spect_litter_prop_alcntr) 


ggplot(litter_rda_pyr_df, aes(x = Trt, y = RDA1))+
  # geom_point(aes(col = Site))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site), outlier.colour = "white")+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site_Trt, scales = "free_x", space = "free_x")+
  # lims(y = c(-.9, 1.5)) +
  labs(x = NULL, y = get_PCA_axislab(litter_pyr_rda))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

cntr_mean <- litter_rda_pyr_df %>% 
  filter(treatment == "Control") %>% 
  group_by(Site_Trt) %>% 
  summarise(RDA1_cntr = mean(RDA1)) %>% 
  ungroup()

litter_rda_std <- litter_rda_pyr_df %>% 
  filter(treatment == "Fertilised") %>% 
  left_join(cntr_mean) %>% 
  mutate(std_RDA1 = RDA1 - RDA1_cntr,
         Total_N = ifelse(Site == "Rosinedal", 1000, Total_N),
         Nadd_yr = ifelse(Site == "Rosinedal", 71.4, Nadd_yr))
summary(litter_rda_std)

cntr_CN <- litter_rda_pyr_df %>% 
  filter(treatment == "Control") %>% 
  group_by(Site) %>% 
  summarise(CN = mean(CN)) %>% 
  mutate(std_RDA1 = 0)


ggplot(trt_litter, aes(x = Trt, y = wN))+
  geom_boxplot(aes(col = Site), outlier.colour = "white")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")

with(litter_rda_std, cor.test(RDA1, wN))

ggplot(litter_rda_std, aes(x = log(Total_N), y = std_RDA1))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site, group = I(Total_N * Nadd_yr)), outlier.colour = "white")


  

ggplot(litter_rda_std, aes(x = Trt, y = std_RDA1))+
  # geom_point(aes(col = Site))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site), outlier.colour = "white")+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  # lims(y = c(-.9, 1.5)) +
  labs(x = NULL, y = get_PCA_axislab(litter_pyr_rda))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
