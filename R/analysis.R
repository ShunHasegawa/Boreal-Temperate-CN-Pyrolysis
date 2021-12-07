source("R/packages.R")
source("R/functions.R")
source("R/generic_functions.R")
sitecols <- c("#e41a1c", "#377eb8", "#4daf4a", "#cab2d6", "#984ea3", "#ff7f00", "#999914", "#a65628", "#f781bf", "#999999")


# Load data ---------------------------------------------------------------


# . Rosinedal 2018 --------------------------------------------------------

load("Data/Pyrolysis_Rosinedal_Oct2018_litter_spectrum_prop.RData")
load("Data/Pyrolysis_Rosinedal_Oct2018_humus_spectrum_prop.RData")
Ros_spec_all <- ldply(list(Litter = litter_spec_all, Humus = humus_spec_all), .id = "Horizon") %>% 
  filter(location2 %in% c("fertilised:inside", "control")) %>%
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
fileid_d <- read.csv("Data/Pyrolysis_sample_list.csv") %>% 
  rename(Horizon = Layer) %>% 
  mutate_all(.funs = list(as.character)) %>% 
  select(Horizon, fileid, Sample_ID) 


# . env variables ---------------------------------------------------------

# Rosinedal OF
load("Data/IRMS_RosinedalOF_2019_2020.RData") 
ro_irms20_d <- IRMS_raw_d %>% 
  filter(Year == 2020) %>% 
  mutate(Horizon = mapvalues(Horizon, c("FH", "L"), c("Humus", "Litter")))
RosinedalOF_env2020 <- read.csv("Data/Pyrolysis_sample_list_RosinedalOF2020.csv") %>% 
  rename(plot = Plot,
         fileid = ID) %>% 
  mutate(fileid = as.character(fileid),
         treatment = substr(plot, 1, 2)) %>% 
  left_join(ro_irms20_d) %>% 
  filter(treatment %in% c("FC", "IN")) %>% 
  mutate(Site      = "RosinedalOF",
         Treatment = mapvalues(treatment, c("FC", "IN"), c("Control", "Fertilised")),
         TrtID     = paste(Site, Treatment, sep = "_")) %>% 
  select(fileid, Site, Horizon, Treatment, TrtID, wN, wC, d15N, d13C)
some(trt_br_d)

# The other sites from the boreal forests (Åheden, Svartberget, Flakaliden)
load("Data/IRMS.RData") 
trt_br_d <- irms_full %>% 
  mutate_if(is.factor, as.character) %>% 
  mutate(Treatment = ifelse(Trt == "Control", "Control", "Fertilised")) %>% 
  rename(Horizon = Layer) %>% 
  left_join(fileid_d) %>% 
  filter(Site != "Rosinedal_OF") %>% # This is from 2019 and not 2020 so remove
  select(fileid, Site, Horizon, Treatment, TrtID, wN, wC, d15N, d13C) %>% 
  bind_rows(Rosinedal_env) %>% 
  bind_rows(RosinedalOF_env2020)
  
  
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
site_order <- c("Aheden", "Flakaliden", "Rosinedal", "RosinedalOF", "Svartberget", "CA", "FE", "HF", "ME", "NH", "BR", "KA")
trt_order  <- c("Aheden_Control", "Aheden_N3kg", "Aheden_N6kg", "Aheden_N12kg", "Aheden_N50kg", 
                "Flakaliden_Control", "Flakaliden_Fertilised",
                "Rosinedal_Control" , "Rosinedal_Fertilised", 
                "RosinedalOF_Control", "RosinedalOF_Fertilised",
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


# source("R/analysys_litter.R")
# this was previously analysed but now data objects were modified and so the
# script needs updating to be run.

# load, process and alysed the samples from the humus together with the US samples
source("R/analysis_BorealTemperate.R")



# Figs for Environmental variables ----------------------------------------

# outlier was identified above
filter(trt_bt_d, fileid == 211)
trt_bt_d <- filter(trt_bt_d, fileid != 211)

# CN catio
bt_cn_P <- ggplot(trt_bt_d, aes(x = TrtID, y = CNratio))+
  geom_boxplot(aes(col = Site), outlier.colour = "white", size = .3)+
  geom_jitter(aes(col = Site), width = .1, alpha = .7)+
  scale_color_manual(values = c(sitecols, "black", "black"))+
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 5),
        legend.position = "none",
        strip.text.x = element_text(size = 4))+
  labs(y = "C:N ratio", x = NULL)
ggsavePP("Output/Figs/Environment/CNratio", bt_cn_P, 6.5, 3)

# N addition
bt_nad_p <- ggplot(filter(site_dd, Treatment == "Fertilised"), aes(x = N_year, y = N_rate))+
  geom_text(aes(label = Site, col = Site, size = N_added))+
  scale_color_manual(values = c(sitecols, "black", "black"), guide = FALSE)+
  labs(x = "N added period (year)", y = "N added rate (kg ha-1 year-1)")+
  lims(x = c(6, 37))
ggsavePP("Output/Figs/Environment/Nadd", bt_nad_p, 5, 4)
