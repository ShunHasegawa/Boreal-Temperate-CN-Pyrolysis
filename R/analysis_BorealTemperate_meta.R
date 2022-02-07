
# Prepare df --------------------------------------------------------------

# mean, sd and N
lcmean_d1 <- spect_tb_d %>% 
  select(Site, lcratio, Treatment, TrtID) %>% 
  group_by(Site, Treatment, TrtID) %>% 
  summarise_at(.vars = vars(lcratio), .funs = list(M = mean, N = get_n, SD = sd)) %>% 
  ungroup()

# control
lcmean_c <- lcmean_d1 %>% 
  filter(Treatment == "Control") %>% 
  select(Site, M, N, SD) %>% 
  rename(c_M = M, c_N = N, c_SD = SD)

# fertilised
lcmean_f <- lcmean_d1 %>% 
  filter(Treatment != "Control") %>% 
  rename(f_M = M, f_N = N, f_SD = SD) %>% 
  left_join(lcmean_c)  %>% 
  arrange(TrtID)


# treatment df
trt_bt_smmry_d <- trt_bt_d %>%  
  group_by(Site, Treatment, TrtID) %>%
  summarise_at(.vars = vars(CNratio), .funs = list(mean)) %>% 
  ungroup()
trt_bt_smmry_c <- trt_bt_smmry_d %>% 
  filter(Treatment == "Control") %>% 
  rename(CN_cnt = CNratio) %>% 
  select(Site, CN_cnt)
trt_bt_smmry_f <- trt_bt_smmry_d %>% 
  filter(Treatment == "Fertilised") %>%
  left_join(trt_bt_smmry_c) %>% 
  left_join(site_dd) %>% 
  mutate(CN_lrr = log(CNratio/CN_cnt))
  
  
# Log Ratio of Means (ROM)
rom <- escalc(measure = "ROM",
              n1i     = f_N, 
              n2i     = c_N,
              m1i     = f_M,
              m2i     = c_M,
              sd1i    = f_SD,
              sd2i    = c_SD,
              slab    = TrtID,
              data    = lcmean_f) %>% 
  left_join(trt_bt_smmry_f) %>% 
  mutate(lci = yi - 1.96 * sqrt(vi), uci = yi + 1.96 * sqrt(vi)) %>% 
  # turn numeric into factor
  mutate(f_nad = quantileCut(N_added, 3),
         f_nrt = quantileCut(N_rate , 3),
         f_nyr = quantileCut(N_year , 4),
         f_rcn = quantileCut(CN_lrr , 3),
         f_cnr = quantileCut(CNratio, 3),
         sdn = 1:n()) %>% 
  # size will be use for plotting; the larger, the smaller variance
  mutate(biom_col = ifelse(Biome == "Boreal", "gray20", "red"),
         wi = 1/sqrt(vi),
         size = .5 + 3 * (wi - min(wi))/(max(wi) - min(wi)))

res <- rma(yi, vi, data = rom)
model_performance(res)
forest(res, main = "Lignin:Carbohrate", cex = .7)




# Meta-analysis -----------------------------------------------------------

# no moderator
res0 <- rma.mv(yi, vi, data = rom, random = ~1|Site)
forest(res0)
summary(res0)
model_performance(res0)

# with moderator and random factors
mr0  <- rma.mv(yi, vi, data = rom, random = ~1|Site)

# . N_added, CN_cnt, CN_lrr
mr1  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_cnt + CN_lrr, random = ~1|Site)
mr2  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_cnt,          random = ~1|Site)
mr3  <- rma.mv(yi, vi, data = rom, mod = ~ N_added +          CN_lrr, random = ~1|Site)
mr4  <- rma.mv(yi, vi, data = rom, mod = ~           CN_cnt + CN_lrr, random = ~1|Site)
AICc(mr0, mr1, mr2, mr3, mr4)
mr5  <- rma.mv(yi, vi, data = rom, mod = ~ N_added,                   random = ~1|Site)
mr6  <- rma.mv(yi, vi, data = rom, mod = ~           CN_cnt,          random = ~1|Site)
AICc(mr0, mr2, mr5, mr6)
summary(mr5)

# . N_year, N_rate, CN_cnt, CN_lrr
mr7   <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_cnt + CN_lrr, random = ~1|Site)
mr8   <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_cnt,          random = ~1|Site)
mr9   <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate +          CN_lrr, random = ~1|Site)
mr10  <- rma.mv(yi, vi, data = rom, mod = ~ N_year +          CN_cnt + CN_lrr, random = ~1|Site)
mr11  <- rma.mv(yi, vi, data = rom, mod = ~          N_rate + CN_cnt + CN_lrr, random = ~1|Site)
AICc(mr0, mr7, mr8, mr9, mr10, mr11)
mr12  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,                   random = ~1|Site)
mr13  <- rma.mv(yi, vi, data = rom, mod = ~ N_year +          CN_cnt,          random = ~1|Site)
mr13  <- rma.mv(yi, vi, data = rom, mod = ~ N_year +                   CN_lrr, random = ~1|Site)
mr14  <- rma.mv(yi, vi, data = rom, mod = ~          N_rate + CN_cnt,          random = ~1|Site)
mr15  <- rma.mv(yi, vi, data = rom, mod = ~          N_rate +          CN_lrr, random = ~1|Site)
mr16  <- rma.mv(yi, vi, data = rom, mod = ~                   CN_cnt + CN_lrr, random = ~1|Site)
AICc(mr0, mr12, mr13, mr14, mr15, mr16)
mr17  <- rma.mv(yi, vi, data = rom, mod = ~ N_year, random = ~1|Site)
mr18  <- rma.mv(yi, vi, data = rom, mod = ~ N_rate, random = ~1|Site)
mr19  <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt, random = ~1|Site)
mr20  <- rma.mv(yi, vi, data = rom, mod = ~ CN_lrr, random = ~1|Site)
AICc(mr0, mr17, mr18, mr19, mr20)
summary(mr0)

# . N_added, CNratio, CN_lrr
mr21  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CNratio + CN_lrr, random = ~1|Site)
mr22  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CNratio,          random = ~1|Site)
mr23  <- rma.mv(yi, vi, data = rom, mod = ~ N_added +           CN_lrr, random = ~1|Site)
mr24  <- rma.mv(yi, vi, data = rom, mod = ~           CNratio + CN_lrr, random = ~1|Site)
AICc(mr0, mr21, mr22, mr23, mr24)
mr25  <- rma.mv(yi, vi, data = rom, mod = ~ CNratio, random = ~1|Site)
AICc(mr0, mr22, mr25, mr5)
summary(mr5)
summary(mr22)

# . N_year, N_rate, CNratio, CN_lrr
mr26 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CNratio + CN_lrr, random = ~1|Site)
mr27 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CNratio,          random = ~1|Site)
mr28 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate +           CN_lrr, random = ~1|Site)
mr29 <- rma.mv(yi, vi, data = rom, mod = ~ N_year +          CNratio + CN_lrr, random = ~1|Site)
mr30 <- rma.mv(yi, vi, data = rom, mod = ~          N_rate + CNratio + CN_lrr, random = ~1|Site)
AICc(mr0, mr26, mr27, mr28, mr29, mr30)

mr31  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,                    random = ~1|Site)
mr32  <- rma.mv(yi, vi, data = rom, mod = ~ N_year +          CNratio,          random = ~1|Site)
mr32  <- rma.mv(yi, vi, data = rom, mod = ~ N_year +                    CN_lrr, random = ~1|Site)
mr33  <- rma.mv(yi, vi, data = rom, mod = ~          N_rate + CNratio,          random = ~1|Site)
mr34  <- rma.mv(yi, vi, data = rom, mod = ~          N_rate +           CN_lrr, random = ~1|Site)
mr35  <- rma.mv(yi, vi, data = rom, mod = ~                   CNratio + CN_lrr, random = ~1|Site)
mr36  <- rma.mv(yi, vi, data = rom, mod = ~ CNratio,  random = ~1|Site)
AICc(mr0, mr31, mr32, mr33, mr34, mr35, mr36)

# . Best model
AICc(mr0, mr22, mr5)
summary(mr0)
summary(mr22)
summary(mr5)




# without random factor ---------------------------------------------------
mm0 <- rma.mv(yi, vi, data = rom)

# . CN_cnt, N_year, N_rate, CN_lrr 
mm1 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_year + N_rate + CN_lrr)
mm2 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_year + N_rate)
mm3 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_year          + CN_lrr)
mm4 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt          + N_rate + CN_lrr)
mm5 <- rma.mv(yi, vi, data = rom, mod = ~          N_year + N_rate + CN_lrr)
AICc(mm0, mm1, mm2, mm3, mm4, mm5)
summary(mm2)

mm6 <- rma.mv(yi, vi, data = rom, mod = ~          N_year + N_rate)
mm7 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt          + N_rate)
mm8 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_year)
AICc(mm0, mm2, mm6, mm7, mm8)
summary(mm8)

mm9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year)
mm10 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt)
AICc(mm0, mm8, mm9, mm10)
summary(mm10)

# . CN_cnt, N_added, CN_lrr 
mm11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_added + CN_lrr)
mm12 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_added)
mm13 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + CN_lrr)
mm14 <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_lrr)
AICc(mm0, mm11, mm12, mm13, mm14)
summary(mm11)

mm15 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_added)
mm16 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt +           CN_lrr)
mm17 <- rma.mv(yi, vi, data = rom, mod = ~          N_added + CN_lrr)
AICc(mm0, mm11, mm15, mm16, mm17)
summary(mm11)
mm11_lm <- lm(yi ~ CN_cnt + N_added + CN_lrr, rom)
car::vif(mm11_lm)

#. N_year, CNratio, N_rate
mm18 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio + N_year + N_rate)
AICc(mm0, mm18)

mm19 <- rma.mv(yi, vi, data = rom, mod = ~           N_year + N_rate)
mm20 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio          + N_rate)
mm21 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio + N_year)
AICc(mm18, mm19, mm20, mm21)

mm22 <- rma.mv(yi, vi, data = rom, mod = ~ N_year)
mm23 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio)
AICc(mm21, mm22, mm23)


#. N_added, CNratio
mm24 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio + N_added)
mm25 <- rma.mv(yi, vi, data = rom, mod = ~ N_added)
AICc(mm23, mm24, mm25)


# best model
AICc(mm11, mm10, mm23, mm24)
summary(mm24)


# figure with predicted values
range(rom$CNratio)
range(rom$N_added)
mm24_pred_cn <- predict(mm24, cbind(seq(18, 49, length.out = 100), mean(rom$N_added)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi = pred,
         CNratio = X.CNratio,
         N_added = X.N_added)
mm24_pred_an <- predict(mm24, cbind(mean(rom$CNratio), seq(48, 2000, length.out = 100)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi = pred,
         CNratio = X.CNratio,
         N_added = X.N_added)

mm24_CN_p <- ggplot(rom, aes(x = CNratio, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm24_pred_cn, col = "blue")+
  geom_line(data = mm24_pred_cn, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm24_pred_cn, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = "none")+
  labs(x = "C:N ratio", y = "Log RR of lignin:carbohydrate")+
  theme(legend.position   = c(.2, .85),
        legend.title      = element_blank(),
        legend.background = element_blank())

mm24_an_p <- ggplot(rom, aes(x = N_added, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm24_pred_an, col = "blue")+
  geom_line(data = mm24_pred_an, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm24_pred_an, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = "none")+
  labs(x = "Total added N (kg ha-1)", y = "Log RR of lignin:carbohydrate")+
  theme(legend.position   = c(.2, .85),
        legend.title      = element_blank(),
        legend.background = element_blank())

mm24_p <- cbind(ggplotGrob(mm24_CN_p), 
                ggplotGrob(mm24_an_p))  
# grid.newpage()
# grid.draw(mm24_p)
ggsavePP("Output/Figs/metaanalysis_model5", mm24_p, 6.5, 3.5)




# without random factor + Biome ---------------------------------------------------

ggplot(rom, aes(x = log(N_added), y = yi, col = Biome))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(rom, aes(x = N_added, y = yi, col = Biome))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(rom, aes(x = N_added, y = yi, col = Vegetation))+
  facet_grid(. ~ Biome)+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)

ggplot(rom, aes(x = N_rate, y = yi, col = Vegetation))+
  facet_grid(. ~ Biome)+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)
ggplot(rom, aes(x = N_year, y = yi, col = Biome))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(rom, aes(x = log(N_year), y = yi, col = Biome))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(rom, aes(x = CN_lrr, y = yi, col = Biome))+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(rom, aes(x = MAT, y = yi, col = Biome))+
  geom_point()+
  geom_smooth(method = "lm")

names(rom)
boxplot(CN_cnt ~ Biome, rom)
boxplot(N_deposition ~ Biome, rom)
boxplot(N_rate ~ Biome, rom)
boxplot(N_year ~ Biome, rom)
boxplot(N_added ~ Biome, rom)
boxplot(CN_lrr ~ Biome, rom)


mb0 <- rma.mv(yi, vi, data = rom)

# . CN_cnt, N_year, N_rate, CN_lrr 
mb1 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year + N_rate + CN_lrr + CN_cnt)
AICc(mb0, mb1)
mb2 <- rma.mv(yi, vi, data = rom, mod = ~ Biome          + N_rate + CN_lrr + CN_cnt)
mb3 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year          + CN_lrr + CN_cnt)
mb4 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year + N_rate          + CN_cnt)
mb5 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year + N_rate + CN_lrr)
AICc(mb0, mb2, mb3, mb4, mb5)

mb6 <- rma.mv(yi, vi, data = rom, mod = ~ Biome          + N_rate + CN_lrr)
mb7 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year          + CN_lrr)
mb8 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year + N_rate)
AICc(mb0, mb6, mb7, mb8)

mb9  <- rma.mv(yi, vi, data = rom, mod = ~ Biome          + N_rate)
mb10 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year)
AICc(mb0, mb9, mb10)

mb11 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_added + I(N_added^2))
mb12 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_added)
mb13 <- rma.mv(yi, vi, data = rom, mod = ~         N_added + I(N_added^2))
AICc(mb11, mb12, mb13)


mb11 <- rma.mv(yi, vi, data = rom, mod = ~ Biome)

# without Flakaliden
rom_noFl <- filter(rom, Site != "Flakaliden")
mf0 <- rma.mv(yi, vi, data = rom_noFl, mod = ~ Biome + N_year + N_rate)
mf1 <- rma.mv(yi, vi, data = rom_noFl, mod = ~ Biome + N_year)
mf2 <- rma.mv(yi, vi, data = rom_noFl, mod = ~ Biome          +  N_rate)
mf3 <- rma.mv(yi, vi, data = rom_noFl, mod = ~         N_year + N_rate)
AICc(mf0, mf1, mf2, mf3)
mf4 <- rma.mv(yi, vi, data = rom_noFl, mod = ~ Biome)
mf5 <- rma.mv(yi, vi, data = rom_noFl, mod = ~ N_year)
AICc(mf0, mf1, mf4, mf5)
summary(mf4)

ma1 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_rate)
summary(ma1)

# # figure with predicted values
# range(rom$N_year)
# unique(rom$Biome)
# coef(summary(mb2))
# mb2_pred_bim <- predict(mb2, cbind(c(0, 1), mean(rom$N_year)), addx = TRUE) %>% 
#   data.frame(.) %>%
#   rename(yi = pred,
#          N_year = X.N_year) %>% 
#   mutate(Biome = ifelse(X.BiomeTemperate == 0, "Boreal", "Temperate"))
# mb2_pred_nyr <- predict(mb2, cbind(c(0, 1), rep(seq(7, 33, length.out = 100), each = 2)), addx = TRUE) %>% 
#   data.frame(.) %>%
#   rename(yi = pred,
#          N_year = X.N_year) %>% 
#   mutate(Biome = ifelse(X.BiomeTemperate == 0, "Boreal", "Temperate"))
# 
# mb2_Nyr_p <- ggplot(rom, aes(x = N_year, y = yi))+
#   geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
#   geom_line(aes(col = Biome), data = mb2_pred_nyr)+
#   geom_line(data = mb2_pred_nyr, aes(y = ci.lb, col = Biome), linetype = "dashed")+
#   geom_line(data = mb2_pred_nyr, aes(y = ci.ub, col = Biome), linetype = "dashed")+
#   geom_point(aes(fill = Biome, size = wi), alpha = .5, shape = 21)+
#   scale_fill_manual(values = c("black", "red"))+
#   scale_color_manual(values = c("black", "red"))+
#   scale_size_continuous(range = c(1, 10), guide = "none")+
#   ylim(-.3, 1.1)+
#   theme(legend.position = "none")+
#   labs(x = "N-added period (year)", y = "Log RR of lignin:carbohydrate")+
#   theme(legend.position   = c(.2, .85),
#         legend.title      = element_blank(),
#         legend.background = element_blank())
# 
# ggsavePP("Output/Figs/metaanalysis_biome_Nyr", mb2_Nyr_p, 3.5, 3.5)


# By Biome ----------------------------------------------------------------
res_bm <- rma.mv(yi, vi, data = rom, mods = ~ Biome - 1, random = ~ 1|Site)
summary(res_bm)

bm_N <- rom %>% 
  group_by(Biome) %>% 
  summarise(N = n())
res_bm_cfd <- coef(summary(res_bm)) %>% 
  data.frame(.) %>%
  bind_cols(bm_N)
with(res_bm_cfd, forest(x = estimate, sei = se, slab= Biome,
                        ilab=paste0("(", N ,")"), ilab.xpos=-.4, 
                        main = "By biome"))




# By Vegetation -----------------------------------------------------------

res_vg <- rma.mv(yi, vi, data = rom, mods = ~ Vegetation - 1, random = ~ 1|Site)
summary(res_vg)

vg_N <- rom %>% 
  group_by(Vegetation) %>% 
  summarise(N = n())
res_vg_cfd <- coef(summary(res_vg)) %>% 
  data.frame(.) %>%
  bind_cols(vg_N)
with(res_vg_cfd, forest(x = estimate, sei = se, slab= Vegetation,
                        ilab=paste0("(", N ,")"), ilab.xpos=-.4))




# By Fertiliser -----------------------------------------------------------


res_ft <- rma.mv(yi, vi, data = rom, mods = ~ Fertiliser - 1, random = ~ 1|Site)
summary(res_ft)

ft_N <- rom %>% 
  group_by(Fertiliser) %>% 
  summarise(N = n())
res_ft_cfd <- coef(summary(res_ft)) %>% 
  data.frame(.) %>%
  bind_cols(ft_N)
with(res_ft_cfd, forest(x = estimate, sei = se, slab= Fertiliser,
                        ilab=paste0("(", N ,")"), ilab.xpos=-.4))




# By year -----------------------------------------------------------------

plot(yi ~ f_nyr, rom)

res_nyr <- rma.mv(yi, vi, data = rom, mods = ~ f_nyr - 1, random = ~ 1|Site)
nyr_N <- rom %>% 
  group_by(f_nyr) %>% 
  summarise(N = n())
res_nyr_cfd <- coef(summary(res_nyr)) %>% 
  data.frame(.) %>%
  bind_cols(nyr_N)
with(res_nyr_cfd, forest(x = estimate, sei = se, slab= f_nyr,
                        ilab=paste0("(", N ,")"), ilab.xpos=-.6))




# By N rate ---------------------------------------------------------------

plot(yi ~ f_nrt, rom)
res_nrt <- rma.mv(yi, vi, data = rom, mods = ~ f_nrt - 1, random = ~ 1|Site)
nrt_N <- rom %>% 
  group_by(f_nrt) %>% 
  summarise(N = n())
res_nrt_cfd <- coef(summary(res_nrt)) %>% 
  data.frame(.) %>%
  bind_cols(nrt_N)
with(res_nrt_cfd, forest(x = estimate, sei = se, slab= f_nrt,
                         ilab=paste0("(", N ,")"), ilab.xpos=-.6))




# By added N --------------------------------------------------------------

plot(yi ~ f_nad, rom)
rom %>% 
  select(Site, TrtID, f_nad) %>% 
  distinct() %>% 
  arrange(f_nad)
res_nad <- rma.mv(yi, vi, data = rom, mods = ~ f_nad - 1, random = ~ 1|Site)
nad_N <- rom %>% 
  group_by(f_nad) %>% 
  summarise(N = n())
res_nad_cfd <- coef(summary(res_nad)) %>% 
  data.frame(.) %>%
  bind_cols(nad_N)
with(res_nad_cfd, forest(x = estimate, sei = se, slab= f_nad,
                         ilab=paste0("(", N ,")"), ilab.xpos=-.6))




# By CN ------------------------------------------------

plot(yi ~ f_cnr, rom)
res_cnr <- rma.mv(yi, vi, data = rom, mods = ~ f_cnr -1 , random = ~ 1|Site)
summary(res_cnr)
cnr_N <- rom %>% 
  group_by(f_cnr) %>% 
  summarise(N = n())
res_cnr_cfd <- coef(summary(res_cnr)) %>% 
  data.frame(.) %>%
  bind_cols(cnr_N)
with(res_cnr_cfd, forest(x = estimate, sei = se, slab= f_cnr,
                         ilab=paste0("(", N ,")"), ilab.xpos=-.3))


# By reseponse ratio of CN ------------------------------------------------

plot(yi ~ f_rcn, rom)
res_rcn <- rma.mv(yi, vi, data = rom, mods = ~ f_rcn -1 , random = ~ 1|Site)
summary(res_rcn)
rcn_N <- rom %>% 
  group_by(f_rcn) %>% 
  summarise(N = n())
res_rcn_cfd <- coef(summary(res_rcn)) %>% 
  data.frame(.) %>%
  bind_cols(rcn_N)
with(res_rcn_cfd, forest(x = estimate, sei = se, slab= f_rcn,
                         ilab=paste0("(", N ,")"), ilab.xpos=-.3))




# Forest plots ------------------------------------------------------------


pdf(file = "Output/Figs/Forestplot_all.pdf", width = 6, height = 3.5)
forest(res, main = "Lignin:Carbohrate", cex = .7, top = 1)
with(res_bm_cfd, forest(x = estimate, sei = se, slab= Biome     , ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By biome"))
with(res_vg_cfd, forest(x = estimate, sei = se, slab= Vegetation, ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By vegetation"))
with(res_ft_cfd, forest(x = estimate, sei = se, slab= Fertiliser, ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By fertiliser"))
with(res_nyr_cfd, forest(x = estimate, sei = se, slab= f_nyr    , ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By No. of years"))
with(res_nrt_cfd, forest(x = estimate, sei = se, slab= f_nrt    , ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By N add rate"))
with(res_nad_cfd, forest(x = estimate, sei = se, slab= f_nad    , ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By added N"))
with(res_cnr_cfd, forest(x = estimate, sei = se, slab= f_cnr    , ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By CN ratio"))
with(res_rcn_cfd, forest(x = estimate, sei = se, slab= f_rcn    , ilab=paste0("(", N ,")"), ilab.xpos=-.6, top = 1, cex = .7, main = "By log response ratio of CN"))
dev.off()

pdf(file = "Output/Figs/Forestplot_byBiome.pdf", width = 5, height = 2)
par(mar = c(4, 4, 2, 3))
with(res_bm_cfd, forest(x = estimate, sei = se, slab= Biome     , ilab=paste0("(", N ,")"), ilab.xpos=-.3, top = 1, cex = .7, main = NULL))
dev.off()

pdf(file = "Output/Figs/Forestplot_bySite.pdf", width = 5, height = 4)
forest(res, main = "Lignin:Carbohrate", cex = .7, top = 1)
dev.off()

# Subset analysis ---------------------------------------------------------

bmr1  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_lrr,        random = ~1|Site, subset = Biome == "Boreal")
bmr2  <- rma.mv(yi, vi, data = rom, mod = ~ N_added,                 random = ~1|Site, subset = Biome == "Boreal")
bmr3  <- rma.mv(yi, vi, data = rom, mod = ~ CN_lrr,                  random = ~1|Site, subset = Biome == "Boreal")
bmr4  <- rma.mv(yi, vi, data = rom,                                  random = ~1|Site, subset = Biome == "Boreal")
bmr5  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_lrr,random = ~1|Site, subset = Biome == "Boreal")
bmr6  <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + CN_lrr,         random = ~1|Site, subset = Biome == "Boreal")
bmr7  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + CN_lrr,         random = ~1|Site, subset = Biome == "Boreal")
bmr8  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,         random = ~1|Site, subset = Biome == "Boreal")
bmr9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year,                  random = ~1|Site, subset = Biome == "Boreal")
bmr10 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate,                  random = ~1|Site, subset = Biome == "Boreal")
bmr11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt,                random = ~1|Site, subset = Biome == "Boreal")
AICc(bmr1, bmr2, bmr3, bmr4, bmr5, bmr6, bmr7, bmr8, bmr9, bmr10, bmr11)

bm0 <- rma.mv(yi, vi, data = rom, mod = ~ 1, subset = Biome == "Boreal")
bm1 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_lrr + CNratio, subset = Biome == "Boreal")
bm2 <- rma.mv(yi, vi, data = rom, mod = ~          N_rate + CN_lrr + CNratio, subset = Biome == "Boreal")
bm3 <- rma.mv(yi, vi, data = rom, mod = ~ N_year          + CN_lrr + CNratio, subset = Biome == "Boreal")
bm4 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate          + CNratio, subset = Biome == "Boreal")
bm5 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_lrr          , subset = Biome == "Boreal")
AICc(bm0, bm1, bm2, bm3, bm4, bm5)

bm6 <- rma.mv(yi, vi, data = rom, mod = ~                   CN_lrr + CNratio, subset = Biome == "Boreal")
bm7 <- rma.mv(yi, vi, data = rom, mod = ~ N_year                   + CNratio, subset = Biome == "Boreal")
bm8 <- rma.mv(yi, vi, data = rom, mod = ~ N_year          + CN_lrr          , subset = Biome == "Boreal")
AICc(bm0, bm3, bm6, bm7, bm8)

bm9 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio, subset = Biome == "Boreal")
bm10 <- rma.mv(yi, vi, data = rom, mod = ~ N_year, subset = Biome == "Boreal")
AICc(bm0, bm7, bm9, bm10)




tmr1  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_lrr,         subset = Biome == "Temperate")
tmr2  <- rma.mv(yi, vi, data = rom, mod = ~ N_added,                 subset = Biome == "Temperate")
tmr3  <- rma.mv(yi, vi, data = rom, mod = ~ CN_lrr,                   subset = Biome == "Temperate")
tmr4  <- rma.mv(yi, vi, data = rom,                                  subset = Biome == "Temperate")
tmr5  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_lrr, subset = Biome == "Temperate")
tmr6  <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + CN_lrr,          subset = Biome == "Temperate")
tmr7  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + CN_lrr,          subset = Biome == "Temperate")
tmr8  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,         subset = Biome == "Temperate")
tmr9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year,                  subset = Biome == "Temperate")
tmr10 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate,                  subset = Biome == "Temperate")
AICc(tmr1, tmr2, tmr3, tmr4, tmr5, tmr6, tmr7, tmr8, tmr9, tmr10, tmr11)
tmr11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt,                subset = Biome == "Temperate")


rom_us <- filter(rom, Biome == "Temperate")
names(rom_us)
tmr_mm0 <- rma.mv(yi, vi, data = rom_us)
forest(tmr_mm0, slab = rom_us$Site)
tmr_mm1 <- rma.mv(yi, vi, data = rom_us, mod = ~ MAT)
tmr_mm2 <- rma.mv(yi, vi, data = rom_us, mod = ~ N_deposition)
tmr_mm3 <- rma.mv(yi, vi, data = rom_us, mod = ~ N_added)
tmr_mm4 <- rma.mv(yi, vi, data = rom_us, mod = ~ CNratio)
tmr_mm5 <- rma.mv(yi, vi, data = rom_us, mod = ~ CN_cnt)
tmr_mm6 <- rma.mv(yi, vi, data = rom_us, mod = ~ Vegetation)
tmr_mm7 <- rma.mv(yi, vi, data = rom_us, mod = ~ CN_lrr)
tmr_mm8 <- rma.mv(yi, vi, data = rom_us, mod = ~ MAP)
tmr_mm9 <- rma.mv(yi, vi, data = rom_us, mod = ~ Latitude)
tmr_mm10 <- rma.mv(yi, vi, data = rom_us, mod = ~ N_year)
tmr_mm11 <- rma.mv(yi, vi, data = rom_us, mod = ~ N_rate)
tmr_mm12 <- rma.mv(yi, vi, data = rom_us, mod = ~ Longitude)
tmr_mm13 <- rma.mv(yi, vi, data = rom_us, mod = ~ Fertiliser)
AICc(tmr_mm0, tmr_mm1, tmr_mm2, tmr_mm3, tmr_mm4, tmr_mm5, tmr_mm6, tmr_mm7, 
     tmr_mm8, tmr_mm9, tmr_mm10, tmr_mm11, tmr_mm12, tmr_mm13)
