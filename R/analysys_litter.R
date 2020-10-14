# Litter analysis ---------------------------------------------------------

# Svartberget, Aheden, Flakaliden —litter—
comp_litter <- ldply(c('Aheden'      = "Data/Compound_Litter_Aheden.csv",
                       'Svartberget' = "Data/Compound_Litter_Svartberget.csv",
                       'Flakaliden'  = "Data/Compound_Litter_Flakaliden.csv"),
                     read.csv, .id = "Site") %>%
  select(Site, Window, Group)

spect_litter <- ldply(c('Aheden'      = "Data/Spectra_Litter_Aheden.csv",
                        'Svartberget' = "Data/Spectra_Litter_Svartberget.csv",
                        'Flakaliden'  = "Data/Spectra_Litter_Flakaliden.csv"),
                      function(x){
                        d <- read.csv(x) %>% 
                          gather(key = "variable", "value", starts_with("X")) 
                        return(d)
                      },
                      .id = "Site") %>% 
  select(-RI, -RT_s) %>% 
  left_join(comp_litter) %>% 
  mutate(Layer = "Litter",
         fileid = as.character(as.numeric(gsub("X", "", variable)))) %>% 
  group_by(Site, variable, Group, Layer, fileid) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()


# proportion for biproducts (CO2, sulfur comp) and unk
spect_litter_allprop <- spect_litter %>% 
  group_by(Site, variable) %>% 
  mutate(prop = value/sum(value)) %>% 
  group_by(Site, Group) %>% 
  summarise(prop = mean(prop))
filter(spect_litter_allprop, Group %in% c("co2", "sulfur_comp", "unk"))

# Calculate prop without biproducts and unk
spect_litter_prop <- spect_litter %>% 
  filter(!(Group %in% c("CO2", "sulfur_comp", "unk"))) %>% 
  mutate(grp = ifelse(Group %in% c("n-diketone", "n-ketone", "n-alkane", "n-alkene", "n-alkanal", "n-FA"), "Long_chain_aliphatic",
                      ifelse(Group %in% c("chlorophyll", "steroid", "vitamin"), "Others", 
                             as.character(Group))),
         grp = gsub("-", "_", grp)) %>% 
  group_by(fileid, Layer, Site, variable) %>% 
  mutate(prop = value/sum(value)) %>% 
  group_by(fileid, Layer, Site, grp) %>% 
  summarise(prop = sum(prop)) %>% 
  ungroup() %>% 
  spread(grp, prop) %>% 
  bind_rows(Rosinedal_litter_raw) %>% 
  left_join(trt_dd) %>% 
  mutate(Trt  = factor(Trt, 
                       levels = c("Control", "Fertilised", "N1", 
                                  "N2", paste0("N", c(3, 6, 12, 50), "kg"))),
         lcratio = (g_lignin + s_lignin + Phenol)/carbohydrate) 



# RDA ---------------------------------------------------------------------
litter_pyr_sp <- decostand(select(spect_litter_prop, aromatic:s_lignin), method = "hellinger")
litter_pyr_rda <- rda(litter_pyr_sp ~ treatment * Site + Condition(Site), spect_litter_prop)

summary(litter_pyr_rda)

anova(litter_pyr_rda, strata = spect_litter_prop$Site)
anova(litter_pyr_rda, strata = spect_litter_prop$Site, by = "margin")
summary(litter_pyr_rda)
plot(litter_pyr_rda)
ordispider(litter_pyr_rda, interaction(spect_litter_prop$Site, spect_litter_prop$Trt), label = TRUE, 
           cex = .4, col = rep(palette()[1:5], 2))

litter_rda_pyr_df <- data.frame(scores(litter_pyr_rda, choices = 1, display = "sites", scaling = 3)) %>% 
  bind_cols(spect_litter_prop) 
litter_rda_pyr_site_p <- ggplot(litter_rda_pyr_df, aes(x = Trt, y = -RDA1))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site), outlier.colour = "white", size = .3)+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  lims(y = c(-.4, .7)) +
  labs(x = NULL, y = paste("L horizon", get_PCA_axislab(litter_pyr_rda)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Sp score
pyr_litter_rda_sp   <- data.frame(scores(litter_pyr_rda, choices = 1, display = "species", scaling = 3)) %>% 
  mutate(pyr_comp = row.names(.),
         pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
         pyr_comp_ed = mapvalues(pyr_comp, 
                                 c("aromatic", "carbohydrate", "g_lignin",
                                   "Long_chain_aliphatic", "N_comp", "Others",
                                   "Phenol", "s_lignin"),
                                 c("Aromatic", "Carbohydrate", "G lignin",
                                   "Long-chain aliphatic", "N comp.", "Others",
                                   "Phenol", "S lignin"))) %>% 
  arrange(RDA1)


litter_rda_pyrsp_p <- ggplot(pyr_litter_rda_sp, aes(x = -1, y = -RDA1, label = pyr_comp_ed)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -1.1) +
  geom_point(aes(x = -1.1), size = 1) +
  geom_text(hjust  = 0,  size = 4) +
  labs(x = NULL, y = NULL) +
  lims(x = c(-1.15, 1.11), y = c(-.4, .7)) +
  science_theme +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
litter_rda_pyrsp_p

litter_rda_pyr_p <- ggarrange(litter_rda_pyr_site_p, litter_rda_pyrsp_p, ncol = 2, widths = c(2, .9))
litter_rda_pyr_p

save(litter_rda_pyr_df, pyr_litter_rda_sp, litter_pyr_sp, litter_pyr_rda,
     file = "Output/Data/Pyr_RDA.RDA")


ggplot(litter_rda_pyr_df, aes(x = d13C, y = RDA1, col = Site))+
  geom_point()+
  geom_smooth(method = lm)


# . Standardise RDA ---------------------------------------------------------
litter_rda_pyr_df_cntr <- litter_rda_pyr_df %>%
  filter(treatment == "Control") %>% 
  group_by(Site) %>% 
  summarise(RDA1_cntr = mean(RDA1)) %>% 
  ungroup()

litter_rda_pyr_df_std <- litter_rda_pyr_df %>% 
  filter(treatment == "Fertilised") %>% 
  left_join(litter_rda_pyr_df_cntr) %>% 
  mutate(RDA1_std = RDA1 - RDA1_cntr) %>% 
  group_by(Site, Trt, Total_N, Nadd_yr, Vegetation, Duration) %>% 
  summarise_at(.vars = vars(RDA1_std, CN), .funs = funs(M = mean, SE = se))

ggplot(litter_rda_pyr_df_std, aes(x = Total_N, y = RDA1_std_M))+
  geom_point(aes(col = Site))+
  geom_errorbar(aes(ymin = RDA1_std_M - RDA1_std_SE, ymax = RDA1_std_M + RDA1_std_SE, col = Site))+
  geom_hline(yintercept = 0, linetype = "dotted")




# NMDS --------------------------------------------------------------------

litter_pyr_nmds <- metaMDS(litter_pyr_sp, distance = "bray")

pdf(file = "Output/Figs/d13C_analysis/Pyrolysis_NMDS_Lhoriz.pdf", width = 6, height = 6)
plot(litter_pyr_nmds, main = "Pyrolysis L horizon")
abline(v = 0, h = 0, col = "gray40", lty = 2)
ordispider(litter_pyr_nmds, paste(spect_litter_prop$Site, spect_litter_prop$Trt, sep = ":"),
           label = TRUE, 
           cex = .4, col = c(rep(palette()[1], 5), rep(palette()[2], 2),
                             rep(palette()[3], 3), rep(palette()[4], 2)))
plot(envfit(litter_pyr_nmds ~ CN + d13C + wN + wC, spect_litter_prop), col = "purple")
text(litter_pyr_nmds, display = "species", col = "brown")
dev.off()



# d13C --------------------------------------------------------------------

litter_pyr_rda_13c <- rda(litter_pyr_sp ~ d13C + wN + Condition(Site), spect_litter_prop)
anova(litter_pyr_rda_13c, by = "margin")

pdf(file = "Output/Figs/d13C_analysis/Pyrolysis_RDA_Lhoriz_d13C_wN.pdf", width = 6, height = 6)
plot(litter_pyr_rda_13c, scaling = 3, main = "Pyrolysis (L horizon) RDA vs. d13C, wN")
text(litter_pyr_rda_13c, display = "species", col = "brown", scaling = 3)
dev.off()

# Picea abies
litter_PA <- filter(spect_litter_prop, Vegetation == "Picea abies")
litter_PA_sp <-decostand(select(litter_PA, aromatic:s_lignin), method = "hellinger")
litter_pyr_rda_13c_pa <- rda(litter_PA_sp ~ d13C + wN + Condition(Site), litter_PA)
anova(litter_pyr_rda_13c_pa)
anova(litter_pyr_rda_13c_pa, by = "axis")

pdf(file = "Output/Figs/d13C_analysis/Pyrolysis_RDA_Lhoriz_Pabies_d13C_wN.pdf", width = 6, height = 6)
plot(litter_pyr_rda_13c_pa, scaling = 3, main = "Pyrolysis (L horizon, Picea abies) RDA vs. d13C, wN")
dev.off()

# Pynus sylvestris
litter_PS <- filter(spect_litter_prop, Vegetation == "Pynus sylvestris")
litter_PA_sp <-decostand(select(litter_PS, aromatic:s_lignin), method = "hellinger")
litter_pyr_rda_13c_ps <- rda(litter_PA_sp ~ d13C + wN + Condition(Site), litter_PS)
anova(litter_pyr_rda_13c_ps)
anova(litter_pyr_rda_13c_ps, by = "axis")

pdf(file = "Output/Figs/d13C_analysis/Pyrolysis_RDA_Lhoriz_Psylvestris_d13C_wN.pdf", width = 6, height = 6)
plot(litter_pyr_rda_13c_ps, scaling = 3, main = "Pyrolysis (L horizon, Pynus sylvestris) RDA vs. d13C, wN")
dev.off()
