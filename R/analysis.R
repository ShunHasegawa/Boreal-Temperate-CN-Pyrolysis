source("R/packages.R")
source("R/functions.R")
source("R/generic_functions.R")

sitecols <- brewer.pal(4, "Dark2")
load("Data/IRMS.RData")
irms_full <- irms_full %>% 
  mutate_at(.vars = c("Layer", "Sample_ID"), .funs = funs(as.character))
fileid_d <- read.csv("Data/Pyrolysis_sample_list.csv", sep = ";") %>% 
  mutate_all(.funs = funs(as.character)) %>% 
  select(Layer, fileid, Sample_ID) 
trt_dd <- full_join(irms_full, fileid_d) %>% 
  mutate(treatment = ifelse(Trt == "Control", "Control", "Fertilised"))
trt_litter <- trt_dd %>% 
  filter(Layer == "Litter")
trt_humus <- trt_dd %>% 
  filter(Layer == "Humus")

# Rosinedal 2018
load("Data/Pyrolysis_Rosinedal_Oct2018_litter_spectrum_prop.RData")
names(litter_spec_all)
names(litter_spec_all)
Rosinedal_env <- litter_spec_all %>%
  select(variable, treatment, location, wN, d15N, d13C, CNratio, leaf_d15N) %>% 
  group_by(variable, treatment) %>% 
  mutate_at(.vars = vars("wN", "d15N", "d13C", "CNratio", "leaf_d15N"), .funs = funs(mean))

Rosinedal_litter_raw <- litter_spec_all %>% 
  filter(location2 %in% c("fertilised:inside", "control")) %>% 
  mutate(grp = ifelse(grp %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene"), "Long_chain_aliphatic",
                      ifelse(grp %in% c("chlorophyll", "steroid", "vitamin"), "Others", grp)),
         treatment = ifelse(location2 == "fertilised:inside", "Fertilised", "Control")) %>% 
  rename(Layer = layer) %>% 
  select(grp, comp, variable, value, Layer, treatment, variable, CNratio, -location2) %>% 
  group_by(grp, treatment, CNratio, Layer) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(grp, value) %>% 
  droplevels() %>% 
  rename(CN = CNratio) %>% 
  mutate(Site = "Rosinedal", 
         Trt = treatment)


# Svartberget, Aheden, Flakaliden —litter—
comp_litter <- ldply(c('Aheden'      = "Data/Compound_Litter_Aheden.csv",
                       'Svartberget' = "Data/Compound_Litter_Svartberget.csv",
                       'Flakaliden'  = "Data/Comound_Litter_Flakaliden.csv"),
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
  left_join(trt_dd)

# proportion for biproducts (CO2, sulfur comp) and unk
spect_litter_allprop <- spect_litter %>% 
  group_by(Site, variable) %>% 
  mutate(prop = value/sum(value)) %>% 
  group_by(Site, Group) %>% 
  summarise(prop = mean(prop))
filter(spect_litter_allprop, Group  %in% c("co2", "sulfur_comp", "unk"))

# Calculate prop withough biproducts and unk
spect_litter_prop <- spect_litter %>% 
  filter(!(Group %in% c("co2", "sulfur_comp", "unk"))) %>% 
  group_by(Site, variable, Layer, Trt) %>% 
  mutate(prop = value/sum(value),
         grp = ifelse(Group %in% c("n-diketone", "n-ketone", "n-alkane", "n-alkene", "n-alkanal", "n-FA"), "Long_chain_aliphatic",
                      ifelse(Group %in% c("chlorophyll", "steroid", "vitamin"), "Others", 
                             as.character(Group))),
         grp = gsub("-", "_", grp)) %>% 
  group_by(Site, variable, grp, CN, treatment, Layer, Trt) %>% 
  summarise(prop = sum(prop)) %>% 
  ungroup() %>% 
  select(Layer, Site, variable, CN, grp, treatment, Trt, prop) %>% 
  spread(grp, prop) %>% 
  bind_rows(Rosinedal_litter_raw) %>% 
  mutate(Trt  = factor(Trt, 
                       levels = c("Control", "Fertilised", "IN_100", "N1", 
                                  "N2", paste0("N", c(3, 6, 12, 50), "kg"))))



# RDA ---------------------------------------------------------------------
litter_sp <- decostand(select(spect_litter_prop, aromatic:s_lignin), method = "hellinger")
litter_rda <- rda(litter_sp ~ treatment * Site + Condition(Site), spect_litter_prop)

anova(litter_rda, strata = spect_litter_prop$Site)
summary(litter_rda)
plot(litter_rda)
ordispider(litter_rda, interaction(spect_litter_prop$Site, spect_litter_prop$treatment), label = TRUE, 
           cex = .4, col = rep(palette()[1:5], 2))

litter_rda_df <- data.frame(scores(litter_rda)$sites) %>% 
  bind_cols(spect_litter_prop) 
litter_rda_pyr_site_p <- ggplot(litter_rda_df, aes(x = Trt, y = -RDA1))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_boxplot(aes(col = Site), outlier.colour = "white")+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = sitecols)+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  lims(y = c(-.9, 1.5)) +
  labs(x = NULL, y = get_PCA_axislab(litter_rda))


# Sp score
pyr_litter_rda_sp   <- data.frame(scores(litter_rda)$species) %>% 
  mutate(pyr_comp = row.names(.),
         pyr_comp = factor(pyr_comp, levels = pyr_comp[order(RDA1)])) %>% 
  arrange(RDA1)


pyr_comp_ed <- c("G lignin", "N comp.", "Phenol", "Aromatic", "Others", "S lignin", "Long-chain aliphatic", "Carbohydrate")
litter_rda_pyrsp_p <- ggplot(pyr_litter_rda_sp, aes(x = -1, y = -RDA1*5, label = pyr_comp_ed)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -1.1) +
  geom_point(aes(x = -1.1), size = 1) +
  geom_text(hjust  = 0,  size = 4) +
  labs(x = "", y = "") +
  lims(x = c(-1.15, 1.11), y = c(-.9, 1.5)) +
  science_theme +
  theme(panel.border = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank())
litter_rda_pyrsp_p

litter_rda_pyr_p <- ggarrange(litter_rda_pyr_site_p, litter_rda_pyrsp_p, ncol = 2, widths = c(2, .9))






