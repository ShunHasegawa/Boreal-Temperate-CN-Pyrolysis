
# Prepare df --------------------------------------------------------------

# Convert lignin:carbonhydrate ratios with a log before computing means (i.e.,
# computed means are geometric means)

# mean, sd and N
lcmean_d1 <- spect_tb_d %>% 
  select(Site, lcratio, Treatment, TrtID) %>% 
  group_by(Site, Treatment, TrtID) %>%
  mutate(loglcratio = log(lcratio)) %>% # get log-scale
  summarise_at(.vars = vars(loglcratio), .funs = list(M = mean, N = get_n, SD = sd)) %>% 
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
  mutate(CN_lrr = log(CNratio/CN_cnt),
         obs = as.character(1:n()))
  
  
# Log Ratio of Means (ROM)
rom <- escalc(measure = "MD", # Means are log-scale. so use mean-difference (MD). it is actually log(lcratio1)-log(lcratio2)=log(lcratio1/lcratio2)
              n1i     = f_N, 
              n2i     = c_N,
              m1i     = f_M,
              m2i     = c_M,
              sd1i    = f_SD,
              sd2i    = c_SD,
              data    = lcmean_f) %>% 
  left_join(trt_bt_smmry_f) %>% 
  mutate(lci = yi - 1.96 * sqrt(vi), uci = yi + 1.96 * sqrt(vi)) %>% 
  # turn numeric into factor
  mutate(f_nad = quantileCut(N_added, 3),
         f_nrt = quantileCut(N_rate , 3),
         f_nyr = quantileCut(N_year , 3),
         f_rcn = quantileCut(CN_lrr , 3),
         f_cnr = quantileCut(CNratio, 3),
         sdn = 1:n()) %>% 
  # size will be use for plotting; the larger, the smaller variance
  mutate(biom_col = ifelse(Biome == "Boreal", "gray20", "red"),
         wi = 1/sqrt(vi),
         size = .5 + 3 * (wi - min(wi))/(max(wi) - min(wi))) %>% 
  group_by(Site) %>% 
  mutate(wi2 = 1/n()) %>%  # Weight based on the number of values per site. Each site (and not each measurement) has the equal weight
  ungroup() %>% 
  mutate(Site = factor(Site,levels = site_order))



# Meta-analysis -----------------------------------------------------------


# . Glabal effects ---------------------------------------------------------

mt0 <- rma.mv(yi, vi, random = ~ 1|Site, data = rom)
summary(mt0)
forest(mt0, slab = TrtID, main = "log(Lignin:Carbohrate)", cex = .7)
lrr2perc(coef(mt0))




# . Biome -----------------------------------------------------------------
mt_bio <- rma.mv(yi, vi, mod = ~ Biome - 1, random = ~ 1|Site, data = rom)
bio_N <- rom %>% 
  group_by(Biome) %>% 
  summarise(N = n())
bio_d <- coef(summary(mt_bio)) %>% 
  mutate(type = "Biome",
         grp = levels(factor(rom$Biome)),
         size = bio_N$N)
lrr2perc(coef(mt_bio))

# .  N added period -------------------------------------------------------
mt_nyr <- rma.mv(yi, vi, mod = ~ f_nyr - 1, random = ~ 1|Site, data = rom)
nyr_N <- rom %>% 
  group_by(f_nyr) %>% 
  summarise(N = n())
nyr_d <- coef(summary(mt_nyr)) %>% 
  mutate(type = "N added years",
         grp = levels(rom$f_nyr),
         size = nyr_N$N)


# . N addition rates ------------------------------------------------------
mt_nrt <- rma.mv(yi, vi, mod = ~ f_nrt - 1, random = ~ 1|Site, data = rom)
nrt_N <- rom %>% 
  group_by(f_nrt) %>% 
  summarise(N = n())
nrt_d <- coef(summary(mt_nrt)) %>% 
  mutate(type = "N added rates",
         grp = levels(rom$f_nrt),
         size = nrt_N$N)



# . N added ---------------------------------------------------------------
mt_Nadd <- rma.mv(yi, vi, data = rom, mod = ~ N_added + I(N_added^2), random = ~1|Site/obs)
N_added_val <- seq(min(rom$N_added), max(rom$N_added), .1)
newmod <- cbind(N_added = N_added_val, N_added2 = N_added_val^2)
pred_d <- predict(mt_Nadd, newmods = newmod) %>% 
  data.frame(.) %>% 
  bind_cols(newmod) %>% 
  rename(yi = pred)
Nadd_p <- ggplot(pred_d, aes(x = N_added, y =  yi))+
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "gray80", alpha = .5)+
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1)+
  geom_line()+
  geom_point(data = rom, aes(y = yi, fill = Biome), 
             size = weights(mt_Nadd)/2, alpha = .7, shape = 21, col = "black")+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = expression(Total~added~N~(kg~ha^'-1')), y = expression(LRR~of~lignin:carbohydrate))+
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 2)))+
  theme(legend.position = c(.25, .8),
        legend.title = element_blank(),
        legend.key.height = unit(.8, units = "lines"))
ggsavePP(filename = "Output/Figs/Metaregression_Nadded", width = 3, height = 3, Nadd_p)


# . Forest plot -----------------------------------------------------------

meta_d <- bind_rows(bio_d, nyr_d, nrt_d)

create_forestplot <- function(){
  par(mar=c(4, 4, 1, 2))
  with(meta_d, forest(x = estimate, sei = se, slab = grp, annotate = FALSE,
                      xlim = c(-1, 1), ylim = c(-2, 17),
                      ilab = paste0("(", size, ")"), 
                      ilab.xpos = -.4, psize = 1.5, xlab = "LRR of lignin:carbohydrate",
                      rows = c(13:12, 9:7, 4:2), cex = .75 
  )
  )
  text(-1, c(14, 10, 5), pos = 4, cex = .75, 
       c(expression(bold(Biome)), 
         expression(bold(N~added~duration~(year))),
         expression(bold(N~added~rate~(kg~ha^'-1'~year^'-1')))))
  addpoly(mt0, row = -1, cex = 1, mlab="", annotate = FALSE)
  text(-1, -1, pos = 4, cex = .75, expression(bold(Summary)))
  
}


png(file = "Output/Figs/Summary_Forest_plot.png", width = 3, height = 5, res = 600, units = "in")
create_forestplot()
dev.off()


pdf(file = "Output/Figs/Summary_Forest_plot.pdf", width = 3, height = 5)
create_forestplot
dev.off()
