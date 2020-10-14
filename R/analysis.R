source("R/packages.R")
source("R/functions.R")
source("R/generic_functions.R")
sitecols <- brewer.pal(4, "Dark2")


# Load data ---------------------------------------------------------------

# Rosinedal 2018
load("Data/Pyrolysis_Rosinedal_Oct2018_litter_spectrum_prop.RData")
load("Data/Pyrolysis_Rosinedal_Oct2018_humus_spectrum_prop.RData")
Ros_spec_all <- ldply(list(Litter = litter_spec_all, Humus = humus_spec_all), .id = "Layer")
names(Ros_spec_all)

names(Ros_spec_all)
Rosinedal_env <- Ros_spec_all %>%
  rename(CN = CNratio,
         DW = soilweight) %>% 
  filter(location2 %in% c("fertilised:inside", "control")) %>%
  mutate(Vegetation = "Pynus sylvestris",
         Site = "Rosinedal",
         treatment = mapvalues(treatment, c("control", "fertilised"), c("Control", "Fertilised")),
         Trt = treatment,
         TrtID = paste(Site, Trt, sep = "_"),
         Vegetation = "Pynus sylvestris",
         Nadd_yr = ifelse(treatment == "Control", 0, 73.08),
         Start_yr = 2006,
         End_yr = 2018,
         Duration = 12,
         Total_N = ifelse(treatment == "Control", 0, 950)) %>% 
  select(wN, d15N, FN, wC, d13C, FC, Layer, Site, Nadd_yr, Start_yr, End_yr,Vegetation, Trt, DW, Duration, Total_N, CN, fileid, treatment, TrtID) %>% 
  distinct() # remove duplicated rows

# check there is no duplicate in fileid
all(!duplicated(Rosinedal_env$fileid))

# IRMS and other sites
load("Data/IRMS.RData")
irms_full <- irms_full %>% 
  mutate_if(is.factor, as.character)
fileid_d <- read.csv("Data/Pyrolysis_sample_list.csv", sep = ";") %>% 
  mutate_all(.funs = funs(as.character)) %>% 
  select(Layer, fileid, Sample_ID) 

# Merge
trt_dd <- full_join(irms_full, fileid_d) %>% 
  mutate(treatment = ifelse(Trt == "Control", "Control", "Fertilised")) %>% 
  select(-Sample_ID, -Plate.ID, -Plot) %>% 
  bind_rows(Rosinedal_env)




# Multivariate analysis ---------------------------------------------------

# Rosinedal spectum data
Rosinedal_spec_raw <- Ros_spec_all %>% 
  filter(location2 %in% c("fertilised:inside", "control")) %>%
  mutate(grp = ifelse(grp %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "Long_chain_aliphatic",
                      ifelse(grp %in% c("chlorophyll", "steroid", "vitamin", "hopanoid"), "Others", grp)),
         Site = "Rosinedal",
         fileid = as.character(fileid)) %>% 
  select(fileid, Layer, Site, grp, value) %>% 
  group_by(Layer, Site, grp, fileid) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  spread(grp, value)
rowSums(select(Rosinedal_spec_raw, aromatic:s_lignin))
Rosinedal_litter_raw <- filter(Rosinedal_spec_raw, Layer == "Litter")
Rosinedal_humus_raw <- filter(Rosinedal_spec_raw, Layer == "Humus")


source("R/analysys_litter.R")
source("R/analysys_humus.R")


# Merge RDA figures
rda_pyr_p <- ggarrange(litter_rda_pyr_site_p + 
                         theme(axis.text.x = element_blank(),
                               plot.margin = margin(t = .2, r = 0, b = .5, l = .5, unit = "line")),
                       litter_rda_pyrsp_p +
                         theme(plot.margin = margin(t = .2, r = 0, b = .5, l = .2, unit = "line")),
                       humus_rda_pyr_site_p + 
                         theme(strip.background = element_blank(),
                               strip.text = element_blank(),
                               plot.margin = margin(t = 0, r = 0, b = .5, l = .5, unit = "line")),  
                       humus_rda_pyrsp_p +
                         theme(plot.margin = margin(t = 0, r = 0, b = .5, l = .2, unit = "line")), 
                       ncol = 2, widths = c(2, .9), nrow = 2, heights = c(1, 1))
rda_pyr_p
ggsavePP(filename = "Output/Figs/RDA_Pyrolysis_HumusLitter", rda_pyr_p, 
       width = 8, height= 5)



# Lignin:Carbohydrate ratios ----------------------------------------------

lcr_litter <- spect_litter_prop %>% 
  mutate(lcr = g_lignin/carbohydrate)
ggplot(lcr_litter, aes(x = Trt, y = lcr))+
  geom_boxplot(aes(col = Site))+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")

ggplot(lcr_litter, aes(x = Trt, y = N_comp))+
  geom_boxplot(aes(col = Site))+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")

lcr_litter_cntr <- lcr_litter %>% 
  filter(treatment == "Control") %>% 
  group_by(Site) %>% 
  summarise(lcr_cntr = mean(lcr),
            CN_cntr = mean(CN)) %>% 
  ungroup()

lcr_litter_fert <- lcr_litter %>% 
  filter(treatment == "Fertilised") %>% 
  left_join(lcr_litter_cntr) %>% 
  mutate(logRR = log(lcr/lcr_cntr))

ggplot(lcr_litter_fert, aes(x = log(Total_N), y = logRR))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_boxplot(aes(col = Site, group = Total_N))

ggplot(lcr_litter_fert, aes(x = Duration, y = logRR))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_boxplot(aes(col = Site, group = Duration))

lcr_humus <- spect_humus_prop %>% 
  mutate(lcr = g_lignin/carbohydrate)
ggplot(lcr_humus, aes(x = Trt, y = lcr))+
  geom_boxplot(aes(col = Site))+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")

ggplot(lcr_humus, aes(x = Trt, y = N_comp))+
  geom_boxplot(aes(col = Site))+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")


lcr_humus_cntr <- lcr_humus %>% 
  filter(treatment == "Control") %>% 
  group_by(Site) %>% 
  summarise(lcr_cntr = mean(lcr),
            CN_cntr = mean(CN)) %>% 
  ungroup()

lcr_humus_fert <- lcr_humus %>% 
  filter(treatment == "Fertilised") %>% 
  left_join(lcr_humus_cntr) %>% 
  mutate(logRR = log(lcr/lcr_cntr))

ggplot(lcr_humus_fert, aes(x = log(Total_N), y = logRR))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_boxplot(aes(col = Site, group = Total_N))


lcr_full <- bind_rows(lcr_litter_fert, lcr_humus_fert) %>% 
  mutate(layer2 = mapvalues(Layer, c("Humus", "Litter"), c("F/H horizon", "L horizon")),
         Layer2 = factor(layer2, levels = c("L horizon", "F/H horizon")))
lcr_fig <- ggplot(lcr_full, aes(x = Total_N, y = logRR))+
  geom_hline(yintercept = 0, linetype = "dotted")+
  geom_boxplot(aes(col = Site, group = Total_N), size = .3)+
  facet_grid(Layer2 ~ .)+
  scale_color_manual(values = sitecols)+
  labs(x = "Total added N (kg)", y = "log RR of ligin:carbohydrate")+
  theme(legend.position = "none")
ggsavePP("Output/Figs/lignin_carbohydrate_ratio", lcr_fig, width = 6, height = 6)


cn_all <- ggplot(lcr_full, aes(x = Total_N, y = CN))+
  geom_hline(aes(yintercept = CN_cntr, col = Site), alpha = .7)+
  geom_boxplot(aes(col = Site, group = Total_N))+
  facet_grid(Layer2 ~ .)+
  scale_color_manual(values = sitecols)+
  labs(x = "Total added N (kg)", y = "C:N ratio")
ggsavePP("Output/Figs/CNratio_totalN", cn_all, width = 6, height = 4)