source("R/packages.R")
source("R/functions.R")
source("R/generic_functions.R")
sitecols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#999914", "#a65628", "#f781bf", "#999999")


# Load data ---------------------------------------------------------------


# . Rosinedal 2018 --------------------------------------------------------

load("Data/Pyrolysis_Rosinedal_Oct2018_litter_spectrum_prop.RData")
load("Data/Pyrolysis_Rosinedal_Oct2018_humus_spectrum_prop.RData")
Ros_spec_all <- ldply(list(Litter = litter_spec_all, Humus = humus_spec_all), .id = "Horizon") %>% 
  # filter(location2 %in% c("fertilised:inside", "control")) %>% 
  mutate(fileid     = as.character(fileid),
         Treatment  = mapvalues(treatment, c("control", "fertilised"), c("Control", "Fertilised")),
         Vegetation = "Pynus sylvestris",
         Site       = "Rosinedal",
         TrtID      = paste(Site, Treatment, sep = "_"),
         Group      = ifelse(grp %in% c("carboxylic_acid", "ketone", "n_alkane", "n_alkene", "other_aliphatics"), "Long_chain_aliphatic",
                             ifelse(grp %in% c("chlorophyll", "steroid", "vitamin", "hopanoid"), "Others", grp))) %>% 
  select(Horizon, Site, fileid, Treatment, TrtID, Group, comp, value, wN, wC, d15N, d13C)

# Spectrum data
Rosinedal_spec_raw <- Ros_spec_all %>% 
  select(Site, Horizon, fileid, Group, value) %>% 
  group_by(Site, Horizon, fileid, Group) %>% 
  summarise(value = sum(value), .groups = "drop") %>% 
  spread(Group, value)

Rosinedal_litter_raw <- filter(Rosinedal_spec_raw, Horizon == "Litter")
Rosinedal_humus_raw  <- filter(Rosinedal_spec_raw, Horizon == "Humus")



# .  Environmental variables for the other sites --------------------------

# site information
site_dd <- read.csv("Data/Site_info.csv") %>% 
  mutate(TrtID = paste(Site, Trt, sep = "_")) %>% 
  select(-Trt)

# Boreal forests

# Rosinedal
Rosinedal_env <- Ros_spec_all %>%
  select(Horizon, Site, fileid, Treatment, TrtID, wN, wC, d15N, d13C) %>% 
  distinct() # remove duplicated rows

# check there is no duplicate in fileid
all(!duplicated(Rosinedal_env$fileid))
xtabs(~ Treatment + Horizon, Rosinedal_env)


# Pyr sample id for the other boreal sites
fileid_d <- read.csv("Data/Pyrolysis_sample_list.csv", sep = ";") %>% 
  rename(Horizon = Layer) %>% 
  mutate_all(.funs = list(as.character)) %>% 
  select(Horizon, fileid, Sample_ID) 

# env variables
load("Data/IRMS.RData")
trt_br_d <- irms_full %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(Treatment = ifelse(Trt == "Control", "Control", "Fertilised")) %>% 
  rename(Horizon = Layer) %>% 
  left_join(fileid_d) %>% 
  filter(Site != "Rosinedal_OF") %>%
  select(fileid, Site, Horizon, Treatment, TrtID, wN, wC, d15N, d13C) %>% 
  bind_rows(Rosinedal_env)
  
# Temperate forests (US samples)
us_sample_list <- read.csv("Data/US_Pyr_list.csv") %>% 
  mutate_all(as.character)

trt_tm_d <- read.csv("Data/US_sample_list.csv") %>% 
  filter(Horizon == "O") %>% 
  rename(Site = Site_ID) %>%
  mutate(Replicate  = as.character(Replicate),
         Treatment  = ifelse(Treatment == "A", "Control", "Fertilised"),
         TrtID      = paste(Site, Treatment, sep = "_")) %>% 
  left_join(us_sample_list) %>% 
  select(Site, Horizon, fileid, Treatment, TrtID, wC, wN, d15N, d13C)

# merge boreal and tempearte forest site information
site_order <- c("Aheden", "Flakaliden", "Rosinedal", "Svartberget", "CA", "FE", "HF", "ME", "NH", "BR", "KA")
trt_order  <- c("Aheden_Control", "Aheden_N3kg", "Aheden_N6kg", "Aheden_N12kg", "Aheden_N50kg", 
                "Flakaliden_Control", "Flakaliden_Fertilised",
                "Rosinedal_Control" , "Rosinedal_Fertilised", 
                "Svartberget_Control", "Svartberget_N1", "Svartberget_N2", 
                "CA_Control", "CA_Fertilised",
                "FE_Control", "FE_Fertilised",
                "HF_Control", "HF_Fertilised",
                "ME_Control", "ME_Fertilised",
                "NH_Control", "NH_Fertilised",
                "BR_Control", "KA_Control")
trt_bt_d <- trt_br_d %>% 
  filter(Horizon == "Humus") %>%   # Only the sample from the humus (F/H) horizon will be used for the anlysis with the US samples
  bind_rows(trt_tm_d) %>% 
  left_join(site_dd) %>% 
  mutate(CNratio = wC/wN)%>% 
  select(fileid, Biome, Site, Vegetation, Treatment, TrtID, wC, wN, d15N, d13C, CNratio, everything(), -Horizon) %>% 
  mutate(Site = factor(Site,levels = site_order),
         TrtID = factor(TrtID, levels = trt_order))




# Multivariate analysis ---------------------------------------------------


source("R/analysys_litter.R")
source("R/analysis_humus.R")


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
names(lcr_litter)


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


lc_all <- rbind(lcr_humus, lcr_litter) %>% 
  mutate(Horizon = ifelse(Layer == "Humus", "F/H horizon", "L horizon"),
         Horizon = factor(Horizon, levels = c("L horizon", "F/H horizon"))) %>% 
  group_by(Horizon, Site, Trt) %>% 
  summarise_at(.vars = vars(lcr), .funs = list(M = mean, SE = se, N = get_n))

lcratio_p <- ggplot(lc_all, aes(x = Trt, y = M))+
  geom_bar(aes(fill = Site), stat = "identity")+
  geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width = .2, size = .5)+
  facet_grid(Horizon ~ Site, scales = "free", space = "free_x")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  labs(x = NULL, y = "Lignin:Carbohydrate ratio")+
  scale_fill_manual(values = sitecols)
lcratio_p
ggsavePP(filename = "Output/Figs/LC_ratio_bySite", width = 6, height = 6,
         plot = lcratio_p)
