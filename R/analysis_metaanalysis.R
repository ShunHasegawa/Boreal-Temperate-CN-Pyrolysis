library(metafor)
names(spect_tb_prop)
some(site_dd)
some(trt_dd_tb)
dim(trt_dd_tb)
dim(spect_tb_prop)
plot(spect_tb_prop$lcratio)
plot(spect_tb_prop$carbohydrate)
which.max(spect_tb_prop$lcratio)
which.min(spect_tb_prop$carbohydrate)
spect_tb_prop[23, ]
d1 <- spect_tb_prop[-23, ] %>% 
  select(Site, lcratio, Trt) %>% 
  group_by(Site, Trt) %>% 
  summarise_at(.vars = vars(lcratio), .funs = list(M = mean, N = get_n, SD = sd)) %>% 
  ungroup()

d_c <- d1 %>% 
  filter(Trt == "Control") %>% 
  select(Site, M, N, SD) %>% 
  rename(c_M = M, c_N = N, c_SD = SD)
sid_order <- c("Aheden_N3kg", "Aheden_N6kg", "Aheden_N12kg", "Aheden_N50kg", 
               "Flakaliden_Fertilised", "Rosinedal_Fertilised", 
               "Svartberget_N1", "Svartberget_N2", 
               "CA_Fertilised", "FE_Fertilised", "HF_Fertilised", "ME_Fertilised", "NH_Fertilised")
d2 <- d1 %>% 
  filter(Trt != "Control") %>% 
  rename(f_M = M, f_N = N, f_SD = SD) %>% 
  left_join(d_c)  %>% 
  mutate(sid = paste(Site, Trt, sep = "_"),
         sid = factor(sid, levels = sid_order)) %>% 
  arrange(sid)

trt_frtl_dd_mean <- trt_frtl_dd %>%  
  group_by(Site, Trt) %>%
  summarise_at(.vars = vars(CNratio, CN_cntrl), .funs = funs(mean)) %>% 
  ungroup() %>% 
  mutate(rr_CN = log(CNratio/CN_cntrl))
  
  
# 1. Log Ratio of Means (ROM)
rom <- escalc(measure = "ROM",
              n1i     = f_N, 
              n2i     = c_N,
              m1i     = f_M,
              m2i     = c_M,
              sd1i    = f_SD,
              sd2i    = c_SD,
              slab    = sid,
              data    = d2) %>% 
  left_join(site_fert_dd) %>% 
  left_join(trt_frtl_dd_mean) %>% 
  mutate(lci = yi - 1.96 * sqrt(vi), uci = yi + 1.96 * sqrt(vi))

rom %>% 
  select(yi, vi) %>% 
  mutate(tt = yi - 1.96 * sqrt(vi))
res <- rma(yi, vi, data = rom)
model_performance(res)
forest(res, main = "Lignin:Carbohrate")
t1 <- tree::tree(yi ~ N_added, data = rom)
plot(t1)
text(t1)
boxplot(N_added ~  Biome, rom)
boxplot(N_year ~  Biome, rom)
boxplot(N_rate ~  Biome, rom)
boxplot(rr_CN ~  Biome, rom)
abline(h = 775)

# turn numeric into simple factors
rom <- rom %>% 
  mutate(f_nad = quantileCut(N_added, 3),
         f_nrt = quantileCut(N_rate , 3),
         f_nyr = quantileCut(N_year , 4),
         f_rcn = quantileCut(rr_CN , 3),
         sdn = 1:n())
  


# no moderator
res0 <- rma.mv(yi, vi, data = rom, random = ~1|Site)
forest(res0)
summary(res0)
model_performance(res0)

# with moderator and random factors
mr1 <- rma.mv(yi, vi, data = rom, mod = ~ N_added + rr_CN, random = ~1|Site)
mr2 <- rma.mv(yi, vi, data = rom, mod = ~ N_added, random = ~1|Site)
mr3 <- rma.mv(yi, vi, data = rom, mod = ~ rr_CN, random = ~1|Site)
mr4 <- rma.mv(yi, vi, data = rom, random = ~1|Site)
mr5 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + rr_CN, random = ~1|Site)
mr6 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + rr_CN, random = ~1|Site)
mr7 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + rr_CN, random = ~1|Site)
mr8 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate, random = ~1|Site)
mr9 <- rma.mv(yi, vi, data = rom, mod = ~ N_year, random = ~1|Site)
mr10 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate, random = ~1|Site)
mr11 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_added -1, random = ~1|Site)
mr12 <- rma.mv(yi, vi, data = rom, mod = ~ Biome + N_year, random = ~1|Site)
mr13 <- rma.mv(yi, vi, data = rom, mod = ~ Biome - 1, random = ~1|Site)
AICc(mr1, mr2, mr3, mr4, mr5, mr6, mr7, mr8, mr9, m10, mr11, mr12, mr13)

# mr2
summary(mr2)
model_performance(mr2)
mr2_pred <- predict(mr2, newmods = seq(0, 2000, length.out = 100), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(N_added = X.N_added)
rom <- rom %>% 
  mutate(biom_col = ifelse(Biome == "Boreal", "gray20", "red"),
         wi = 1/sqrt(rom$vi),
         size = .5 + 3 * (wi - min(wi))/(max(wi) - min(wi)))
plot(yi ~ N_added, rom, type = "n", pch = 16, ylim = c(-.5, 1), xlab = "Added N (kg ha-1)", ylab = "log response ratio", main = "Lignin:Carbohydrate")
d_ply(rom, .(Biome), function(x) points(yi ~ N_added, x, col = biom_col, pch = 16, cex = size))
lines(pred ~ N_added, mr2_pred, col = "blue")
lines(ci.lb ~ N_added, mr2_pred, lty = "dotted", col = "blue")
lines(ci.ub ~ N_added, mr2_pred, lty = "dotted", col = "blue")
abline(h = 0, lty = "dashed", col = "gray30")
with(distinct(select(rom, Biome, biom_col)), 
     legend(0, 1, legend = Biome, pch = 16, col = biom_col, bty = "n"))


# mr13
summary(mr13)
forest(mr13)



# without random factor ---------------------------------------------------
mm0 <- rma.mv(yi, vi, data = rom)
mm1 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + N_year + N_rate + rr_CN)
mm2 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + N_year + N_rate)
mm3 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + N_year          + rr_CN)
mm4 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl          + N_rate + rr_CN)
mm5 <- rma.mv(yi, vi, data = rom, mod = ~           N_year + N_rate + rr_CN)
AICc(mm0, mm1, mm2, mm3, mm4, mm5)
mm6 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + N_year)
mm7 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + rr_CN)
mm8 <- rma.mv(yi, vi, data = rom, mod = ~ N_year   + rr_CN)
AICc(mm3, mm6, mm7, mm8)
mm9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year)
mm10 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl)
AICc(mm6, mm9, mm10)
summary(mm6)

mm11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + N_added + rr_CN)
mm12 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + N_added)
mm13 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl + rr_CN)
mm14 <- rma.mv(yi, vi, data = rom, mod = ~ N_added + rr_CN)
AICc(mm0, mm11, mm12, mm13, mm14)
mm15 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl)
mm16 <- rma.mv(yi, vi, data = rom, mod = ~ N_added)
AICc(mm12, mm15, mm16)
AICc(mm6, mm15)
summary(mm6)

range(rom$CN_cntrl)
range(rom$N_year)
mm6_pred_cnc <- predict(mm6, cbind(seq(16, 48, length.out = 100), mean(rom$N_year)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi = pred,
         CN_cntrl = X.CN_cntrl)
mm6_pred_nyr <- predict(mm6, cbind(mean(rom$CN_cntrl), seq(7, 33, length.out = 100)), addx = TRUE) %>% 
  data.frame(.) %>%
  rename(yi = pred,
         N_year = X.N_year)

mm6_CNcntrl_p <- ggplot(rom, aes(x = CN_cntrl, y = yi))+
  geom_hline(yintercept = 0, col = "gray30", linetype = "dashed")+
  geom_point(aes(fill = Biome, size = wi), col = "black", alpha = .7, shape = 21)+
  geom_line(data = mm6_pred_cnc, col = "blue")+
  geom_line(data = mm6_pred_cnc, aes(y = ci.lb), col = "blue", linetype = "dotted")+
  geom_line(data = mm6_pred_cnc, aes(y = ci.ub), col = "blue", linetype = "dotted")+
  scale_fill_manual(values = c("black", "red"))+
  scale_size_continuous(range = c(1, 10), guide = FALSE)+
  labs(x = "C:N ratio at control", y = "Log response ratio")+
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
  labs(x = "N-added period (year)", y = "Log response ratio")

mm6_p <- cbind(ggplotGrob(mm6_CNcntrl_p), 
               ggplotGrob(mm6_Nyr_p))  
# grid.newpage()
# grid.draw(mm6_p)
ggsavePP("Output/Figs/metaanalysis_model1", mm6_p, 6.5, 3.5)


# alternative model
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
  labs(x = "C:N ratio", y = "Log response ratio")+
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
  labs(x = "Total added N (kg ha-1)", y = "Log response ratio")

mm17_p <- cbind(ggplotGrob(mm17_cnr_p), 
                ggplotGrob(mm17_nad_p))  
# grid.newpage()
# grid.draw(mm17_p)
ggsavePP("Output/Figs/metaanalysis_model2", mm17_p, 6.5, 3.5)






plot(yi ~ CNratio, rom)
plot(yi ~ CN_cntrl, rom)
plot(CNratio ~ CN_cntrl, rom)
summary(lm(CNratio ~ CN_cntrl + N_rate, rom))
summary(lm(CNratio ~ CN_cntrl + N_added, rom))








plot(yi ~ N_added, rom, col = "red", pch = 16, ylim = c(-1, 1.5))
with(rom, arrows(N_added, lci, N_added, uci, col = "red", length = .05, code = 3, angle = 90))
lines(pred ~ X.N_added, mm2_pred)
lines(ci.lb ~ X.N_added, mm2_pred, lty = "dotted")
lines(ci.ub ~ X.N_added, mm2_pred, lty = "dotted")
abline(h = 0, lty = "dashed")





res1 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + rr_CN)
res2 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + rr_CN)
res3 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + rr_CN)
res4 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate)
AICc(res1,res2, res3, res4)
res5 <- rma.mv(yi, vi, data = rom, mod = ~ N_year)
res6 <- rma.mv(yi, vi, data = rom, mod = ~ rr_CN)
res7 <- rma.mv(yi, vi, data = rom)
AICc(res5,res6, res7)


res1 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + rr_CN, random= ~ 1|Site)
res2 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + rr_CN, random= ~ 1|Site)
res3 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + rr_CN, random= ~ 1|Site)
res4 <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate, random= ~ 1|Site)
AICc(res1,res2, res3, res4)
res5 <- rma.mv(yi, vi, data = rom, mod = ~ N_year, random= ~ 1|Site)
res6 <- rma.mv(yi, vi, data = rom, mod = ~ rr_CN, random= ~ 1|Site)
res7 <- rma.mv(yi, vi, data = rom, random= ~ 1|Site)
AICc(res5,res6, res7)
model_performance(res7)


# Biome
res_bm <- rma.mv(yi, vi, data = rom, mods = ~ Biome - 1, random = ~ 1|Site)
summary(res_bm)

bm_N <- rom %>% 
  group_by(Biome) %>% 
  summarise(N = n())
res_bm_cfd <- coef(summary(res_bm)) %>% 
  data.frame(.) %>%
  bind_cols(bm_N)
with(res_bm_cfd, forest(x = estimate, sei = se, slab= Biome,
                        ilab=paste0("(", N ,")"), ilab.xpos=-.4))


# Vegetation
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


# Fertiliser
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

# N year
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


# N rate
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


# N add
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

# rr_CN
plot(yi ~ f_rcn, rom)
rom %>% 
  select(Site, Trt, f_rcn, yi, rr_CN) %>% 
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

    



# Subset analysis ---------------------------------------------------------

bmr1  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + rr_CN,         random = ~1|Site, subset = Biome == "Boreal")
bmr2  <- rma.mv(yi, vi, data = rom, mod = ~ N_added,                 random = ~1|Site, subset = Biome == "Boreal")
bmr3  <- rma.mv(yi, vi, data = rom, mod = ~ rr_CN,                   random = ~1|Site, subset = Biome == "Boreal")
bmr4  <- rma.mv(yi, vi, data = rom,                                  random = ~1|Site, subset = Biome == "Boreal")
bmr5  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + rr_CN, random = ~1|Site, subset = Biome == "Boreal")
bmr6  <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + rr_CN,          random = ~1|Site, subset = Biome == "Boreal")
bmr7  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + rr_CN,          random = ~1|Site, subset = Biome == "Boreal")
bmr8  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,         random = ~1|Site, subset = Biome == "Boreal")
bmr9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year,                  random = ~1|Site, subset = Biome == "Boreal")
bmr10 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate,                  random = ~1|Site, subset = Biome == "Boreal")
bmr11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl,                random = ~1|Site, subset = Biome == "Boreal")
AICc(bmr1, bmr2, bmr3, bmr4, bmr5, bmr6, bmr7, bmr8, bmr9, bmr10, bmr11)


tmr1  <- rma.mv(yi, vi, data = rom, mod = ~ N_added + rr_CN,         subset = Biome == "Temperate")
tmr2  <- rma.mv(yi, vi, data = rom, mod = ~ N_added,                 subset = Biome == "Temperate")
tmr3  <- rma.mv(yi, vi, data = rom, mod = ~ rr_CN,                   subset = Biome == "Temperate")
tmr4  <- rma.mv(yi, vi, data = rom,                                  subset = Biome == "Temperate")
tmr5  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate + rr_CN, subset = Biome == "Temperate")
tmr6  <- rma.mv(yi, vi, data = rom, mod = ~ N_rate + rr_CN,          subset = Biome == "Temperate")
tmr7  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + rr_CN,          subset = Biome == "Temperate")
tmr8  <- rma.mv(yi, vi, data = rom, mod = ~ N_year + N_rate,         subset = Biome == "Temperate")
tmr9  <- rma.mv(yi, vi, data = rom, mod = ~ N_year,                  subset = Biome == "Temperate")
tmr10 <- rma.mv(yi, vi, data = rom, mod = ~ N_rate,                  subset = Biome == "Temperate")
tmr11 <- rma.mv(yi, vi, data = rom, mod = ~ CN_cntrl,                subset = Biome == "Temperate")
AICc(tmr1, tmr2, tmr3, tmr4, tmr5, tmr6, tmr7, tmr8, tmr9, tmr10, tmr11)






  
  
# res2 <- rma.mv(yi, vi, data = rom, mods = N_added, random = ~1|Site)
res2 <- rma.mv(yi, vi, data = rom, mods = ~ N_added + CN_cntrl + rr_CN)
res_veg <- rma.mv(yi, vi, data = rom, mods = ~ Vegetation)
res_veg_coefd <- data.frame(coef(summary(res_veg))) %>% 
  mutate(variable = levels(factor(rom$Vegetation)))
veg_n <- rom %>% 
  group_by(Vegetation) %>% 
  summarise(N = n()) %>% 
  .$N

forest(x = res_veg_coefd$estimate, sei = res_veg_coefd$se, slab = res_veg_coefd$variable,
       ilab=paste0("(",veg_n,")"),ilab.xpos=-.6)


res_bio <- rma.mv(yi, vi, data = rom, mods = ~ Biome, random = ~ 1|Site)
res_bio_coefd <- data.frame(coef(summary(res_bio))) %>% 
  mutate(variable = levels(factor(rom$Biome)))
bio_n <- rom %>% 
  group_by(Biome) %>% 
  summarise(N = n()) %>% 
  .$N

with(res_bio_coefd, forest(x = estimate, sei = se, slab = variable,
       ilab=paste0("(",bio_n,")"),ilab.xpos=-1))


res_cn <- rma.mv(yi, vi, data = rom, mods = ~ cnc_f)
res_cn_coefd <- data.frame(coef(summary(res_cn))) %>% 
  mutate(variable = levels(factor(rom$cnc_f)))
cnc_n <- rom %>% 
  group_by(cnc_f) %>% 
  summarise(N = n()) %>% 
  .$N

with(res_cn_coefd, forest(x = estimate, sei = se, slab = variable,
                           ilab=paste0("(",cnc_n,")"),ilab.xpos=-1))



res_rcn <- rma.mv(yi, vi, data = rom, mods = ~ rcn_f, random = ~ 1|Site)
res_rcn_coefd <- data.frame(coef(summary(res_rcn))) %>% 
  mutate(variable = levels(factor(rom$rcn_f)))
rcn_n <- rom %>% 
  group_by(rcn_f) %>% 
  summarise(N = n()) %>% 
  .$N

with(res_rcn_coefd, forest(x = estimate, sei = se, slab = variable,
                          ilab=paste0("(",rcn_n,")"),ilab.xpos=-.2))


res_nad <- rma.mv(yi, vi, data = rom, mods = ~ nad_f, random = ~ 1|Site)
res_nad_coefd <- data.frame(coef(summary(res_nad))) %>% 
  mutate(variable = levels(factor(rom$nad_f)))
rnad_n <- rom %>% 
  group_by(nad_f) %>% 
  summarise(N = n()) %>% 
  .$N
with(res_nad_coefd, forest(x = estimate, sei = se, slab = variable,
                           ilab=paste0("(",rnad_n,")"),ilab.xpos=-.2))



with(res_bio_coefd, forest(x = estimate, sei = se, slab = variable,
                           ilab=paste0("(",bio_n,")"),ilab.xpos=-1))




res3 <- rma.mv(yi, vi, data = rom, mods = ~ N_added + CN_cntrl)
res4 <- rma.mv(yi, vi, data = rom, mods = ~ N_added + rr_CN)
res5 <- rma.mv(yi, vi, data = rom, mods = ~ CN_cntrl + rr_CN)
res6 <- rma.mv(yi, vi, data = rom, mods = ~ N_added)
res7 <- rma.mv(yi, vi, data = rom, mods = ~ CN_cntrl)
forest(res7)
res8 <- rma.mv(yi, vi, data = rom, mods = ~ 1)
AICc(res2, res3, res4, res5, res6, res7, res8)
pred_d <- data.frame(predict(res5, addx = TRUE)) %>% 
  arrange(X.CN_cntrl)
plot(pred ~ X.CN_cntrl, pred_d, ylim = c(-.3, 1.1), type = "l")
lines(ci.lb ~ X.CN_cntrl, pred_d, lty = "dashed")
lines(ci.ub ~ X.CN_cntrl, pred_d, lty = "dashed")
points(yi ~ CN_cntrl, col = "red", rom)
with(rom, arrows(x0 = CN_cntrl, y0 = lci, x1 = CN_cntrl, y1 = uci, angle = 90, code = 3, length = .1, col = "red"))
abline(h = 0, lty = "dotted")


res2
anova(res2)
### forest plot with extra annotations

forest(res2)



dat <- dat.bcg

### calculate log risk ratios and corresponding sampling variances (and use
### the 'slab' argument to store study labels as part of the data frame)
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
              slab=paste(author, year, sep=", "))

### fit random-effects model
res <- rma(yi, vi, data=dat)

### forest plot with extra annotations
forest(res, atransf=exp, at=log(c(.05, .25, 1, 4)), xlim=c(-16,6),
       ilab=cbind(dat.bcg$tpos, dat.bcg$tneg, dat.bcg$cpos, dat.bcg$cneg),
       ilab.xpos=c(-9.5,-8,-6,-4.5), cex=.75, header="Author(s) and Year",
       mlab="")
forest(res)



