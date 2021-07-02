
# Lignin:Carbohydrate response ratio

# Prepare df ---------------------------------------------------------------------
names(spect_tb_d)
lc_trt_d <- spect_tb_d %>% 
  filter(Treatment == "Fertilised") %>% 
  select(Biome, Site, Vegetation, Treatment, lcratio, TrtID, wC, wN, d15N, d13C, CNratio, Site_name,
         N_deposition, N_rate, N_year, N_added, Latitude, Longitude, MAT, MAP,
         Sampling_year, N_year_start, N_year_end, Fertiliser) %>% 
  rename(Fertilised = lcratio)

lc_cnt_d <- spect_tb_d %>% 
  filter(Treatment == "Control") %>% 
  select(Site, lcratio, CNratio) %>% 
  group_by(Site) %>% 
  summarise_at(.vars = vars(lcratio, CNratio), .funs = list(mean)) %>% 
  rename(Control = lcratio,
         CN_cnt  = CNratio)
  
lc_rr_d <- lc_trt_d %>% 
  left_join(lc_cnt_d) %>% 
  mutate(lrr   = log(Fertilised/Control), # log response ratio
         CN_rr = CNratio/CN_cnt) %>% 
  select(Site, TrtID, Fertilised, Control, lrr, CN_rr, everything())
some(lc_rr_d)


# analysis ----------------------------------------------------------------

# inspect the data
pairs(lrr ~ N_rate + N_year + N_added + CNratio + CN_cnt + MAT + MAP, data = lc_rr_d, panel = panel.smooth)
boxplot(lrr ~ Vegetation, lc_rr_d)
boxplot(lrr ~ Biome, lc_rr_d)
t1 <- tree::tree(lrr ~ N_rate + N_year + CN_cnt + CN_rr, lc_rr_d)
t2 <- tree::tree(lrr ~ N_rate + N_year + CNratio, lc_rr_d)
t3 <- tree::tree(lrr ~ N_added + CN_cnt + CN_rr, lc_rr_d)
t4 <- tree::tree(lrr ~ N_added + CNratio, lc_rr_d)
plot(t1)
text(t1)
plot(t2)
text(t2)

tp1 <- rpart(lrr ~ N_rate + N_year + CN_cnt + CN_rr, lc_rr_d)
tp2 <- rpart(lrr ~ N_rate + N_year + CNratio, lc_rr_d)
tp3 <- rpart(lrr ~ N_added + CN_cnt + CN_rr, lc_rr_d)
tp4 <- rpart(lrr ~ N_added + CNratio, lc_rr_d)
plot(tp1)
text(tp1)

par(mfrow = c(2, 2))
g1 <- gam(lrr ~ s(N_rate, k = 4) + s(N_year, k = 4) + s(CN_cnt, k = 4) + s(CN_rr, k = 4), data = lc_rr_d)
plot(g1)
g2 <- gam(lrr ~ s(N_added, k = 4) + s(CNratio, k = 4) + s(CN_rr, k = 4), data = lc_rr_d)
plot(g2)

lmm1 <- lmer(lrr ~ N_rate + N_year + CN_cnt + CN_rr + (1|Site), lc_rr_d)
Anova(lmm1, test.statistic = "F")
qqresidPlot(lmm1)

lmm2 <- lmer(lrr ~ N_added + CNratio + (1|Site), lc_rr_d)
Anova(lmm2, test.statistic = "F")
qqresidPlot(lmm2)

# ancova
anc1 <- lmer(log(Fertilised) ~ log(Control) + N_added + CNratio + (1|Site), lc_rr_d)
summary(anc1)
Anova(anc1, test.statistic = "F")
qqresidPlot(anc1)



# Sumamrise by treatment replicates to remove pseudo replicates ------------------------------------------------------------------
names(lc_rr_d)
lc_rr_smmry <- lc_rr_d %>% 
  group_by(Site, TrtID, Biome, Vegetation, Treatment, Site_name, N_deposition, N_rate, N_year, N_added, Latitude, Longitude, MAT, MAP, Sampling_year, N_year_start, N_year_end, Fertiliser, Control, CN_cnt) %>% 
  summarise_at(.vars = vars(Fertilised, CNratio), .funs = list(mean)) %>% 
  ungroup() %>% 
  mutate(lrr = log(Fertilised/Control),
         CN_rr = CNratio/CN_cnt)

# inspect
pairs(lrr ~ N_rate + N_year + N_added + CNratio + CN_cnt + MAT + MAP, data = lc_rr_smmry, panel = panel.smooth)
g1 <- gam(lrr ~ s(N_added, k = 4) + s(CNratio, k = 4) + Vegetation, data = lc_rr_smmry)
par(mfrow = c(2, 2))
plot(g1)
g2 <- gam(lrr ~ s(N_year, k = 4) + s(N_rate, k = 4) + s(CNratio, k = 4) + Vegetation, data = lc_rr_smmry)
plot(g2)

lm3 <- lm(lrr ~ N_added + CNratio + I(N_added^2), lc_rr_smmry)
lm4 <- update(lm3, ~ .-I(N_added^2))
lm5 <- lm(lrr ~ Biome, lc_rr_smmry)
anova(lm3, lm4)
AICc(lm3, lm4, lm5)
summary(lm3)
anova(lm3)
plot(lm3)
visreg(lm3, xvar = "N_added")
