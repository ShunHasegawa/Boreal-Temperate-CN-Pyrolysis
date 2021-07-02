
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
  summarise_at(.vars = vars(CNratio), .funs = funs(mean)) %>% 
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
         wi = 1/sqrt(rom$vi),
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
AICc(mr22, mr25, mr5)
summary(mr5)

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
summary(mm3)

mm6 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_year)
mm7 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + CN_lrr)
mm8 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + CN_lrr)
AICc(mm3, mm6, mm7, mm8)
summary(mm6)

mm9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year)
mm10 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt)
AICc(mm6, mm9, mm10)
summary(mm6)

# . CN_cnt, N_added, CN_lrr 
mm11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_added + CN_lrr)
mm12 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + N_added)
mm13 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt + CN_lrr)
mm14 <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_lrr)
AICc(mm0, mm11, mm12, mm13, mm14)
summary(mm12)

mm15 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt)
mm16 <- rma.mv(yi, vi, data = rom, mod = ~ N_added)
AICc(mm12, mm15, mm16)
summary(mm15)

AICc(mm6, mm15)
summary(mm6)

# figure with predicted values
range(rom$CN_cnt)
range(rom$N_year)
mm6_pred_cnc <- predict(mm6, cbind(seq(16, 48, length.out = 100), mean(rom$N_year)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi = pred,
         CN_cnt = X.CN_cnt)
mm6_pred_nyr <- predict(mm6, cbind(mean(rom$CN_cnt), seq(7, 33, length.out = 100)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi = pred,
         N_year = X.N_year)

mm6_CNcntrl_p <- ggplot(rom, aes(x = CN_cnt, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm6_pred_cnc, col = "blue")+
  geom_line(data = mm6_pred_cnc, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm6_pred_cnc, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = FALSE)+
  labs(x = "C:N ratio at control", y = "Log RR of lignin:carbohydrate")+
  theme(legend.position   = c(.2, .85),
        legend.title      = element_blank(),
        legend.background = element_blank())

mm6_Nyr_p <- ggplot(rom, aes(x = N_year, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm6_pred_nyr, col = "blue")+
  geom_line(data = mm6_pred_nyr, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm6_pred_nyr, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = FALSE)+
  theme(legend.position = "none")+
  labs(x = "N-added period (year)", y = NULL)

mm6_p <- cbind(ggplotGrob(mm6_CNcntrl_p), 
               ggplotGrob(mm6_Nyr_p))  
# grid.newpage()
# grid.draw(mm6_p)
ggsavePP("Output/Figs/metaanalysis_model1", mm6_p, 6.5, 3.5)


# Alternative model
mm17 <- rma.mv(yi, vi, data = rom, mod = ~ CNratio + N_added)
AICc(mm6, mm17)
range(rom$CNratio)
range(rom$N_added)
mm17_pred_cnr <- predict(mm17, cbind(seq(18, 49, length.out = 100), mean(rom$N_added)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi      = pred,
         CNratio = X.CNratio)
mm17_pred_nad <- predict(mm17, cbind(mean(rom$CNratio), seq(48, 2000, length.out = 100)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi      = pred,
         N_added = X.N_added)

mm17_cnr_p <- ggplot(rom, aes(x = CNratio, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm17_pred_cnr, col = "blue")+
  geom_line(data = mm17_pred_cnr, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm17_pred_cnr, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = FALSE)+
  labs(x = "C:N ratio", y = "Log RR of lignin:carbohydrate")+
  theme(legend.position   = c(.2, .85),
        legend.title      = element_blank(),
        legend.background = element_blank())

mm17_nad_p <- ggplot(rom, aes(x = N_added, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm17_pred_nad, col = "blue")+
  geom_line(data = mm17_pred_nad, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm17_pred_nad, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = FALSE)+
  theme(legend.position = "none")+
  labs(x = "Total added N (kg ha-1)", y = NULL)

mm17_p <- cbind(ggplotGrob(mm17_cnr_p), 
                ggplotGrob(mm17_nad_p))  
# grid.newpage()
# grid.draw(mm17_p)
ggsavePP("Output/Figs/metaanalysis_model2", mm17_p, 6.5, 3.5)




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
rom %>% 
  select(Site, Trt, f_nyr, yi) %>% 
  arrange(f_nyr)

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
rom %>% 
  select(Site, Trt, f_nrt) %>% 
  distinct() %>% 
  arrange(f_nrt)
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
  select(Site, Trt, f_nad) %>% 
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
rom %>% 
  select(Site, Trt, f_cnr, yi, CN_lrr) %>% 
  arrange(f_cnr)
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
rom %>% 
  select(Site, Trt, f_rcn, yi, CN_lrr) %>% 
  arrange(f_rcn)
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



# Subset analysis ---------------------------------------------------------

bmr1  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + CN_lrr,         random = ~1|Site, subset = Biome == "Boreal")
bmr2  <- rma.mv(yi, vi, data = rom, mod = ~ N_added,                 random = ~1|Site, subset = Biome == "Boreal")
bmr3  <- rma.mv(yi, vi, data = rom, mod = ~ CN_lrr,                   random = ~1|Site, subset = Biome == "Boreal")
bmr4  <- rma.mv(yi, vi, data = rom,                                  random = ~1|Site, subset = Biome == "Boreal")
bmr5  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + CN_lrr, random = ~1|Site, subset = Biome == "Boreal")
bmr6  <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + CN_lrr,          random = ~1|Site, subset = Biome == "Boreal")
bmr7  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + CN_lrr,          random = ~1|Site, subset = Biome == "Boreal")
bmr8  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,         random = ~1|Site, subset = Biome == "Boreal")
bmr9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year,                  random = ~1|Site, subset = Biome == "Boreal")
bmr10 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate,                  random = ~1|Site, subset = Biome == "Boreal")
bmr11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt,                random = ~1|Site, subset = Biome == "Boreal")
AICc(bmr1, bmr2, bmr3, bmr4, bmr5, bmr6, bmr7, bmr8, bmr9, bmr10, bmr11)


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
tmr11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cnt,                subset = Biome == "Temperate")
AICc(tmr1, tmr2, tmr3, tmr4, tmr5, tmr6, tmr7, tmr8, tmr9, tmr10, tmr11)
