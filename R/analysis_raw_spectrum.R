# Analysing raw spectrum data

litterall_raw <- read.table("Data/Spectra_litter_all.txt", header = TRUE) %>% 
  select(-Window, -RI, -Annotated) %>% 
  mutate(total = rowSums(.[, -1])) %>% 
  mutate_at(.vars = vars(-RT_s), .funs = funs(./total)) %>%
  select(-total) %>% 
  gather(key = "fileid", value = "area", starts_with("X")) %>% 
  mutate(fileid = gsub(pattern = "X", "", fileid),
         layer = "Litter") %>% 
  filter(fileid != "109") %>% 
  mutate(fileid = as.character(as.numeric(gsub(pattern = "_1", "", fileid)))) %>% 
  left_join(trt_litter) %>% 
  spread(RT_s, area)

litterall_sp <- decostand(select(litterall_raw, starts_with("Win")), method = "hellinger")
pcoa_litter_comp <- rda(litterall_sp ~ CN  + Condition(Site), litterall_raw)
anova(pcoa_litter_comp, by = "margin")
summary(pcoa_litter_comp)
plot(pcoa_litter_comp)

pcoa_litter_comp <- rda(litterall_sp)
litter_site <- data.frame(scores(pcoa_litter_comp)$sites) %>% 
  bind_cols(litterall_raw)

ggplot(litter_site, aes(x = CN, y = RDA1, col = Site))+
  geom_point(size = 4)+
  ggtitle("Litter")

ggplot(litter_site, aes(x = PC1, y = PC2, col = Site))+
  geom_point(size = 4)+
  ggtitle("Litter")





pcoa_litter_comp <- rda(litterall_sp)
fit <- envfit(pcoa_litter_comp, litterall_raw$CN, perm = 999)
plot(pcoa_litter_comp)
plot(fit)


?envfit

pcoa_litter_comp <- cmdscale(d = vegdist(litterall_sp, method = "bray"),  eig = TRUE, k = 3)
pcoa_litter_comp_res <- data.frame(pcoa_litter_comp$points) %>% 
  rename(Axis1 = X1, 
         Axis2 = X2,
         Axis3 = X3) %>%
  bind_cols(litterall_raw)


p_litter_comp_pcoa <- ggplot(pcoa_litter_comp_res, aes(x = Axis1, y = Axis2, col = Site, size = CN)) +
  geom_point(alpha = .7) +
  ggtitle("Litter compounds")
p_litter_comp_pcoa




# humu --------------------------------------------------------------------

humusall_raw <- read.table("Data/Spectra_humus_all.txt", header = TRUE) %>% 
  select(-Window, -RI, -Annotated) %>% 
  mutate(total = rowSums(.[, -1])) %>% 
  mutate_at(.vars = vars(-RT_s), .funs = funs(./total)) %>%
  select(-total) %>% 
  gather(key = "fileid", value = "area", starts_with("X")) %>% 
  mutate(fileid = gsub(pattern = "X", "", fileid),
         layer = "humus") %>% 
  filter(!(fileid %in% c("249", "250", "252", "253", "255", "262", "268", "273"))) %>% 
  mutate(fileid = as.character(as.numeric(gsub(pattern = "_1", "", fileid)))) %>% 
  left_join(trt_humus) %>% 
  spread(RT_s, area)

humusall_sp <- decostand(select(humusall_raw, starts_with("Win")), method = "hellinger")
pcoa_humus_comp <- rda(humusall_sp ~ CN  + Condition(Site), humusall_raw)
anova(pcoa_humus_comp, by = "margin")
summary(pcoa_humus_comp)
plot(pcoa_humus_comp)
humus_site <- data.frame(scores(pcoa_humus_comp)$sites) %>% 
  bind_cols(humusall_raw)

ggplot(humus_site, aes(x = CN, y = RDA1, col = Site))+
  geom_point(size = 4)+
  ggtitle("humus")

pcoa_humus_comp <- rda(humusall_sp)
humus_site <- data.frame(scores(pcoa_humus_comp)$sites) %>% 
  bind_cols(humusall_raw)
which.min(humus_site$PC1[humus_site$Site == "Aheden"])
humus_site[humus_site$Site == "Aheden", ][15, ]
humus_site %>% 
  filter(fileid == "273")

ggplot(humus_site, aes(x = PC1, y = PC2, col = Site))+
  geom_point(size = 4)+
  ggtitle("humus")





pcoa_humus_comp <- rda(humusall_sp)
fit <- envfit(pcoa_humus_comp, humusall_raw$CN, perm = 999)
plot(pcoa_humus_comp)
plot(fit)



pcoa_humus_comp <- cmdscale(d = vegdist(humusall_sp, method = "bray"),  eig = TRUE, k = 3)
pcoa_humus_comp_res <- data.frame(pcoa_humus_comp$points) %>% 
  rename(Axis1 = X1, 
         Axis2 = X2,
         Axis3 = X3) %>%
  bind_cols(humusall_raw)


p_humus_comp_pcoa <- ggplot(pcoa_humus_comp_res, aes(x = Axis1, y = Axis2, col = Site, size = CN)) +
  # geom_point(alpha = .7) +
  geom_text(aes(label = fileid))+
  ggtitle("humus compounds")
p_humus_comp_pcoa

