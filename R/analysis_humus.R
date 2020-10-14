# Humus analysis ---------------------------------------------------------

# NOTE!! X273_1 and X268_1 should be checked if they are the same sample (but
# more quantity) as X273 and X268, respectively. There is a small chance that I
# have messed up the order and labelled wrongly.....

# Svartberget, Aheden, Flakaliden —humus—
comp_humus <- ldply(c('Aheden'      = "Data/Compound_humus_Aheden.csv",
                       'Svartberget' = "Data/Compound_humus_Svartberget.csv",
                       'Flakaliden'  = "Data/Compound_humus_Flakaliden.csv"),
                     read.csv, .id = "Site") %>%
  select(Site, Window, Group)

spect_humus <- ldply(c('Aheden'      = "Data/Spectra_humus_Aheden.csv",
                       'Svartberget' = "Data/Spectra_humus_Svartberget.csv",
                       'Flakaliden'  = "Data/Spectra_humus_Flakaliden.csv"),
                     function(x){
                       d <- read.csv(x) %>% 
                         gather(key = "variable", "value", starts_with("X")) 
                       return(d)
                     },
                     .id = "Site") %>% 
  select(-RI, -RT_s) %>% 
  left_join(comp_humus) %>% 
  # X273 showed week signal, so realysed and IDed as X273_1
  filter(variable != "X273") %>% 
  mutate(variable = mapvalues(variable, "X273_1", "X273")) %>% 
  mutate(Layer = "Humus",
         fileid = as.character(as.numeric(gsub("X", "", variable)))) %>% 
  group_by(Site, variable, Group, Layer, fileid) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()

# proportion for biproducts (CO2, sulfur comp) and unk
spect_humus_allprop <- spect_humus %>% 
  group_by(Site, variable) %>% 
  mutate(prop = value/sum(value)) %>% 
  group_by(Site, Group) %>% 
  summarise(prop = mean(prop))
filter(spect_humus_allprop, Group  %in% c("co2", "sulfur_comp", "unk"))

# Calculate prop without biproducts and unk
spect_humus_prop <- spect_humus %>% 
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
  bind_rows(Rosinedal_humus_raw) %>% 
  left_join(trt_dd) %>% 
  mutate(Trt  = factor(Trt, 
                       levels = c("Control", "Fertilised", "N1", 
                                  "N2", paste0("N", c(3, 6, 12, 50), "kg"))),
         lcratio = (g_lignin + s_lignin + Phenol)/carbohydrate) 



# RDA ---------------------------------------------------------------------
humus_pyr_sp <- decostand(select(spect_humus_prop, aromatic:s_lignin), method = "hellinger")
humus_pyr_rda <- rda(humus_pyr_sp ~ treatment * Site + Condition(Site), spect_humus_prop)
# humus_pyr_rda <- rda(humus_pyr_sp ~ Trt * Site + Condition(Site), spect_humus_prop)
summary(humus_pyr_rda)

anova(humus_pyr_rda, strata = spect_humus_prop$Site)
anova(humus_pyr_rda, strata = spect_humus_prop$Site, by = "margin")
anova(humus_pyr_rda, strata = spect_humus_prop$Site, by = "axis")
summary(humus_pyr_rda)
plot(humus_pyr_rda)
ordispider(humus_pyr_rda, interaction(spect_humus_prop$Site, spect_humus_prop$Trt), label = TRUE, 
           cex = .4, col = rep(palette()[1:5], 2))

humus_rda_pyr_df <- data.frame(scores(humus_pyr_rda, choices = 1, display = "sites", scaling = 3)) %>% 
  bind_cols(spect_humus_prop) 

humus_rda_pyr_site_p <- ggplot(humus_rda_pyr_df, aes(x = Trt, y = -RDA1))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site), outlier.colour = "white", size = .3)+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  lims(y = c(-.4, .51)) +
  labs(x = NULL, y = paste("F/H horizon", get_PCA_axislab(humus_pyr_rda)))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Sp score
pyr_humus_rda_sp   <- data.frame(scores(humus_pyr_rda, choices = 1, display = "species", scaling = 3)) %>% 
  mutate(pyr_comp = row.names(.),
         pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
         pyr_comp_ed = mapvalues(pyr_comp, 
                                 c("aromatic", "carbohydrate", "g_lignin", "Long_chain_aliphatic", "N_comp", "Others", "Phenol", "s_lignin"),
                                 c("Aromatic", "Carbohydrate", "G lignin", "Long-chain aliphatic", "N comp.", "Others", "Phenol", "S lignin"))) %>% 
  arrange(RDA1)

humus_rda_pyrsp_p <- ggplot(pyr_humus_rda_sp, aes(x = -1, y = -RDA1, label = pyr_comp_ed)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -1.1) +
  geom_point(aes(x = -1.1), size = 1) +
  geom_text(hjust  = 0,  size = 4) +
  labs(x = NULL, y = NULL) +
  lims(x = c(-1.15, 1.11), y = c(-.4, .51)) +
  science_theme +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
humus_rda_pyrsp_p

humus_rda_pyr_p <- ggarrange(humus_rda_pyr_site_p, humus_rda_pyrsp_p, ncol = 2, widths = c(2, .9))
humus_rda_pyr_p

save(humus_rda_pyr_df, pyr_humus_rda_sp, humus_pyr_sp, humus_pyr_rda,
     file = "Output/Data/Pyr_RDA.RDA")




# NMDS --------------------------------------------------------------------

humus_pyr_nmds <- metaMDS(humus_pyr_sp, distance = "bray")
plot(humus_pyr_nmds, main = "Pyrolysis F/H horizon")
abline(v = 0, h = 0, col = "gray40", lty = 2)
ordispider(humus_pyr_nmds, paste(spect_humus_prop$Site, spect_humus_prop$Trt, sep = ":"),
           label = TRUE, 
           cex = .4, col = c(rep(palette()[1], 5), rep(palette()[2], 2),
                             rep(palette()[3], 3), rep(palette()[4], 2)))
plot(envfit(humus_pyr_nmds ~ CN + d13C + wN + wC, spect_humus_prop), col = "purple")
text(humus_pyr_nmds, display = "species", col = "brown")




# d13C --------------------------------------------------------------------

humus_pyr_rda_13c <- rda(humus_pyr_sp ~ d13C + wN + Condition(Site), spect_humus_prop)
anova(humus_pyr_rda_13c, by = "margin")
plot(humus_pyr_rda_13c, scaling = 3)
text(humus_pyr_rda_13c, display = "species", col = "brown", scaling = 3)


# Picea abies
humus_PA <- filter(spect_humus_prop, Vegetation == "Picea abies")
humus_PA_sp <-decostand(select(humus_PA, aromatic:s_lignin), method = "hellinger")
humus_pyr_rda_13c_pa <- rda(humus_PA_sp ~ d13C + Condition(Site), humus_PA)
anova(humus_pyr_rda_13c_pa)
anova(humus_pyr_rda_13c_pa, by = "axis")
plot(humus_pyr_rda_13c_pa)


# Pynus sylvestris
humus_PS <- filter(spect_humus_prop, Vegetation == "Pynus sylvestris")
humus_PA_sp <-decostand(select(humus_PS, aromatic:s_lignin), method = "hellinger")
humus_pyr_rda_13c_ps <- rda(humus_PA_sp ~ d13C + Condition(Site), humus_PS)
anova(humus_pyr_rda_13c_ps)
anova(humus_pyr_rda_13c_ps, by = "axis")
plot(humus_pyr_rda_13c_ps, scaling = 3)
