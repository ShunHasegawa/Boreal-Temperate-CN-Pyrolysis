
# Boreal forest -----------------------------------------------------------

# NOTE!! X273_1 and X268_1 should be checked if they are the same sample (but
# more quantity) as X273 and X268, respectively. There is a small chance that I
# have messed up the order and labeled wrongly.....

# Svartberget, Aheden, Flakaliden —humus—
br_comp_humus <- ldply(c('Aheden'      = "Data/Compound_humus_Aheden.csv",
                         'Svartberget' = "Data/Compound_humus_Svartberget.csv",
                         'Flakaliden'  = "Data/Compound_humus_Flakaliden.csv"),
                     read.csv, .id = "Site") %>%
  select(Site, Window, Group)

br_spect_humus <- ldply(c('Aheden'      = "Data/Spectra_humus_Aheden.csv",
                          'Svartberget' = "Data/Spectra_humus_Svartberget.csv",
                          'Flakaliden'  = "Data/Spectra_humus_Flakaliden.csv"),
                     function(x){
                       d <- read.csv(x) %>% 
                         gather(key = "variable", "value", starts_with("X")) 
                       return(d)
                     },
                     .id = "Site") %>% 
  select(-RI, -RT_s) %>% 
  left_join(br_comp_humus) %>% 
  # X273 showed weak signal, so re-analysed and IDed as X273_1
  filter(variable != "X273") %>% 
  mutate(variable = mapvalues(variable, "X273_1", "X273")) %>% 
  mutate(fileid = as.character(as.numeric(gsub("X", "", variable)))) %>% 
  group_by(Site, variable, Group, fileid) %>% 
  summarise(value = sum(value), .groups = "drop")


# US samples --------------------------------------------------------------
tm_comp <- ldply(c('BR' = "Data/US_Soil/Compound/Compound_US_BR.csv",
                   'CA' = "Data/US_Soil/Compound/Compound_US_CA.csv",
                   'FE' = "Data/US_Soil/Compound/Compound_US_FE.csv",
                   'HF' = "Data/US_Soil/Compound/Compound_US_HF.csv",
                   'KA' = "Data/US_Soil/Compound/Compound_US_KA.csv",
                   'ME' = "Data/US_Soil/Compound/Compound_US_ME.csv",
                   'NH' = "Data/US_Soil/Compound/Compound_US_NH.csv"),
                 read.csv, .id = "Site") %>%
  select(Site, Window, Group)

tm_spect <- ldply(c('BR' = "Data/US_Soil/Spectra/Spectra_BR.csv",
                    'CA' = "Data/US_Soil/Spectra/Spectra_CA.csv",
                    'FE' = "Data/US_Soil/Spectra/Spectra_FE.csv",
                    'HF' = "Data/US_Soil/Spectra/Spectra_HF.csv",
                    'KA' = "Data/US_Soil/Spectra/Spectra_KA.csv",
                    'ME' = "Data/US_Soil/Spectra/Spectra_ME.csv",
                    'NH' = "Data/US_Soil/Spectra/Spectra_NH.csv"),
                  function(x){
                    d <- read.csv(x) %>% 
                      gather(key = "variable", "value", starts_with("BR")) 
                    return(d)
                  },
                  .id = "Site") %>% 
  select(-RI, -RT_s) %>% 
  left_join(tm_comp) %>% 
  mutate(fileid = as.character(as.numeric(gsub("BR|_.", "", variable)))) %>% 
  group_by(Site, variable, Group, fileid) %>% 
  summarise(value = sum(value), .groups = "drop")


# merge US (temperate) and Boreal samples
spect_tb <- rbind(br_spect_humus, tm_spect)
dim(spect_tb)

# proportion for biproducts (CO2, sulfur comp) and unk
spect_allprop <- spect_tb %>% 
  group_by(Site, variable) %>% 
  mutate(prop = value/sum(value)) %>% 
  group_by(Site, Group) %>% 
  summarise(prop = mean(prop))
filter(spect_allprop, Group  %in% c("co2", "sulfur_comp", "unk"))


# Calculate prop without biproducts and unk and merge with Rosinedal
spect_tb_prop <- spect_tb %>% 
  filter(!(Group %in% c("CO2", "sulfur_comp", "unk"))) %>% 
  mutate(grp = ifelse(Group %in% c("n-diketone", "n-ketone", "n-alkane", "n-alkene", "n-alkanal", "n-FA"), "Long_chain_aliphatic",
                      ifelse(Group %in% c("chlorophyll", "steroid", "vitamin"), "Others", 
                             as.character(Group))),
         grp = gsub("-", "_", grp)) %>% 
  group_by(fileid, Site, variable) %>% 
  mutate(prop = value/sum(value)) %>% 
  group_by(fileid, Site, grp) %>% 
  summarise(prop = sum(prop), .groups = "drop") %>% 
  spread(grp, prop, fill = 0) %>% 
  mutate(g_lignin = g_lignin + gs_lignin) %>% 
  select(-gs_lignin) %>% 
  bind_rows(select(Rosinedal_humus_raw, -Horizon))
rowSums(select(spect_tb_prop, aromatic:s_lignin))


# compute lignin carbohydrate ratios and merge with envrironmental variables
dim(spect_tb_prop)
dim(trt_bt_d)
spect_tb_d <- spect_tb_prop %>%  
  mutate(lcratio = (g_lignin + s_lignin + Phenol)/carbohydrate,
         gsratio = g_lignin / s_lignin) %>% 
  left_join(trt_bt_d) %>%
  mutate(Site = factor(Site, levels = site_order)) %>% 
  filter(!(Site %in% c("BR", "KA")))

# Check outliers
plot(spect_tb_d$lcratio)
which.max(spect_tb_d$lcratio)
spect_tb_d[23, ]
spect_tb_d <- filter(spect_tb_d, fileid != 211)



# RDA ---------------------------------------------------------------------
tb_pyr_sp <- decostand(select(spect_tb_d, aromatic:s_lignin), method = "hellinger")
tb_pyr_rda <- rda(tb_pyr_sp ~ Treatment * Site + Condition(Site), spect_tb_d)
summary(tb_pyr_rda)

anova(tb_pyr_rda, strata = spect_tb_d$Site)
anova(tb_pyr_rda, strata = spect_tb_d$Site, by = "margin")
anova(tb_pyr_rda, strata = spect_tb_d$Site, by = "axis")
summary(tb_pyr_rda)
plot(tb_pyr_rda)
ordispider(tb_pyr_rda, interaction(spect_tb_d$Site, spect_tb_d$TrtID), label = TRUE, cex = .4)


# site scores
tb_rda_pyr_df <- data.frame(scores(tb_pyr_rda, choices = 1, display = "sites", scaling = 3)) %>% 
  bind_cols(spect_tb_d) 

tb_rda_pyr_site_p <- ggplot(tb_rda_pyr_df, aes(x = TrtID, y = -RDA1))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site), outlier.colour = "white", size = .3)+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  labs(x = NULL, y = get_PCA_axislab(tb_pyr_rda))+
  lims(y = c(-.46, .44))+
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 5),
        legend.position = "none",
        strip.text.x = element_text(size = 4))

# Sp score
pyr_tb_rda_sp   <- data.frame(scores(tb_pyr_rda, choices = 1, display = "species", scaling = 3)) %>% 
  mutate(pyr_comp = row.names(.),
         pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)]),
         pyr_comp_ed = mapvalues(pyr_comp, 
                                 c("aromatic", "carbohydrate", "g_lignin", "Long_chain_aliphatic", "N_comp", "Others", "Phenol", "s_lignin"),
                                 c("Aromatic", "Carbohydrate", "G lignin", "Long-chain aliphatic", "N comp.", "Others", "Phenol", "S lignin"))) %>% 
  arrange(RDA1)

tb_rda_pyrsp_p <- ggplot(pyr_tb_rda_sp, aes(x = -1, y = -RDA1, label = pyr_comp_ed)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -1.1) +
  geom_point(aes(x = -1.1), size = 1) +
  geom_text(hjust  = 0,  size = 2) +
  labs(x = NULL, y = NULL) +
  lims(x = c(-1.15, 1.11), y = c(-.46, .45)) +
  science_theme +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
tb_rda_pyrsp_p

tb_rda_pyr_p <- ggarrange(tb_rda_pyr_site_p, tb_rda_pyrsp_p, ncol = 2, widths = c(2, .9))
tb_rda_pyr_p
ggsavePP(filename = "Output/Figs/Pyr_RDA_BorealTemperate", tb_rda_pyr_p, 6.5, 3)


# . By Site ---------------------------------------------------------------

rda_bysite_res <- dlply(spect_tb_d, .(Site), rda_bysite)
llply(rda_bysite_res, function(x) anova(x$model))

# save RDA figs
l_ply(names(rda_bysite_res), function(x){
  p <- rda_bysite_res[[x]]$fig
  ggsavePP(filename = paste0("Output/Figs/Pyr_RDA_", x), p, 3, 3)
})



# NMDS --------------------------------------------------------------------

tb_pyr_nmds <- metaMDS(tb_pyr_sp, distance = "euc", try = 50)
stressplot(tb_pyr_nmds)
pdf(file = "Output/Figs/Pyr_NMDS_BorealTemperate.pdf", width = 6, height = 6)
plot(tb_pyr_nmds, main = "Pyrolysis F/H horizon")
abline(v = 0, h = 0, col = "gray40", lty = 2)
ordispider(tb_pyr_nmds, spect_tb_d$TrtID, label = TRUE, cex = .4)
plot(envfit(tb_pyr_nmds ~ CNratio + d13C + wN + wC, spect_tb_d), col = "purple")
text(tb_pyr_nmds, display = "species", col = "brown")
dev.off()



# Figs -----------------------------------------------------

# Lignin:carbohydrate
bt_lc_p <- ggplot(spect_tb_d, aes(x = TrtID, y = lcratio))+
  geom_boxplot(aes(col = Site), outlier.colour = "white", size = .3)+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 5),
        legend.position = "none",
        strip.text.x = element_text(size = 4))+
  labs(x = NULL, y = "Lignin:Carbohydrate ratio")
ggsavePP("Output/Figs/bt_LC_ratio", bt_lc_p, 6.5, 3)

# N compound
bt_ncmp_p <- ggplot(spect_tb_d, aes(x = TrtID, y = N_comp * 100))+
  geom_boxplot(aes(col = Site), outlier.colour = "white", size = .3)+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 5),
        legend.position = "none",
        strip.text.x = element_text(size = 4)) +
  labs(x = NULL, y = "Pyr-N compound (%)")
ggsavePP("Output/Figs/bt_Pyr_Ncomp", bt_ncmp_p, 6.5, 3)




# Analysis -----------------------------------------------------------
# Linear model
source("R/analysis_BorealTemperate_LM.R") 

# meta analysis
source("R/analysis_BorealTemperate_meta.R") 




# List of all compounds and group -------------------------------------------------

Ros_comp_all <- ldply(list(Litter = litter_spec_all, Humus = humus_spec_all), .id = "Horizon") %>% 
  select(comp, grp) %>% 
  distinct() %>% 
  rename(Comp = comp, Group = grp)

comp_list <- ldply(c('BR' = "Data/US_Soil/Compound/Compound_US_BR.csv",
                     'CA' = "Data/US_Soil/Compound/Compound_US_CA.csv",
                     'FE' = "Data/US_Soil/Compound/Compound_US_FE.csv",
                     'HF' = "Data/US_Soil/Compound/Compound_US_HF.csv",
                     'KA' = "Data/US_Soil/Compound/Compound_US_KA.csv",
                     'ME' = "Data/US_Soil/Compound/Compound_US_ME.csv",
                     'NH' = "Data/US_Soil/Compound/Compound_US_NH.csv",
                     'Aheden'      = "Data/Compound_humus_Aheden.csv",
                     'Svartberget' = "Data/Compound_humus_Svartberget.csv",
                     'Flakaliden'  = "Data/Compound_humus_Flakaliden.csv"),
                   read.csv, .id = "Site") %>% 
  select(Comp, Group) %>% 
  bind_rows(Ros_comp_all) %>% 
  distinct() %>% 
  arrange(Group, Comp)

write.csv(comp_list, "Output/Tables/All_Pyr_Compound_list.csv", row.names = FALSE)

