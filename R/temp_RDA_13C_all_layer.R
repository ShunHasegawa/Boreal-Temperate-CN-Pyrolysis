dd <- rbind(spect_litter_prop, spect_humus_prop) %>% 
  mutate(VL = paste(Vegetation, Layer, sep = "_"))
ggplot(dd, aes(x = Layer, y = d13C, col = Site))+
  geom_boxplot()

# analyse control

dd_cntr <- filter(dd, treatment == "Control")
dd_cntr_sp <-decostand(select(dd_cntr, aromatic:s_lignin), method = "hellinger")
dd_cntr_m1 <- rda(dd_cntr_sp ~ d13C + wC + wN + Vegetation + Layer, dd_cntr)
anova(dd_cntr_m1)
summary(dd_cntr_m1)
anova(dd_cntr_m1, by = "margin")
anova(dd_cntr_m1, by = "axis")
par(mfrow = c(1, 2))
plot(dd_cntr_m1)
text(dd_cntr_m1, display = "species", col = "red")
ordispider(dd_cntr_m1, dd_cntr$VL, label = TRUE, cex = .5, col = "gray70")
plot(dd_cntr_m1, choice = c(2, 3))
text(dd_cntr_m1, display = "species", col = "red", choice = c(2, 3))
ordispider(dd_cntr_m1, dd_cntr$VL, label = TRUE, cex = .5, col = "gray70", choices = c(2, 3))

# Picea abies
dd_PA <- filter(dd, Vegetation == "Picea abies")
dd_PA_sp <-decostand(select(dd_PA, aromatic:s_lignin), method = "hellinger")
dd_pyr_rda_13c_pa <- rda(dd_PA_sp ~ d13C * Layer, dd_PA)
anova(dd_pyr_rda_13c_pa)
anova(dd_pyr_rda_13c_pa, by = "margin")
anova(dd_pyr_rda_13c_pa, by = "axis")
plot(dd_pyr_rda_13c_pa)


# Pynus sylvestris
dd_PS <- filter(dd, Vegetation == "Pynus sylvestris")
dd_PA_sp <-decostand(select(dd_PS, aromatic:s_lignin), method = "hellinger")
dd_pyr_rda_13c_ps <- rda(dd_PA_sp ~ d13C * Layer, dd_PS)
anova(dd_pyr_rda_13c_ps)
anova(dd_pyr_rda_13c_ps, by = "axis")
anova(dd_pyr_rda_13c_ps, by = "margin")


dd_pyr_rda_13c_ps_m2 <- rda(dd_PA_sp ~ d13C + Layer, dd_PS)
anova(dd_pyr_rda_13c_ps_m2)
anova(dd_pyr_rda_13c_ps_m2, by = "axis")
anova(dd_pyr_rda_13c_ps_m2, by = "margin")
plot(dd_pyr_rda_13c_ps_m2)
text(dd_pyr_rda_13c_ps_m2, display = "species", col = "red")

# NMDS
dd_nmds <- metaMDS(dd_sp, distance = "bray")
env_nmds <- envfit(dd_nmds ~ Vegetation + Layer + d13C + wN + wC, dd)

plot(dd_nmds)
plot(env_nmds)
# ordisurf(dd_nmds, dd$d13C, type = "n")
text(dd_nmds, display = "species", col = "red")
# points(dd_nmds, display = "sites", 
#        col = factor(dd$Vegetation),
#        pch = as.numeric(as.factor(dd$Layer)))

abline(h = 0, lty = 2)
abline(v = 0, lty = 2)



# mvabund
# library(mvabund)
# mv_spp <- mvabund(select(dd, aromatic:s_lignin))
# # mv_spp <- mvabund(decostand(select(dd, aromatic:s_lignin), method = "hellinger"))
# par(mar=c(2,10,2,2)) # adjusts the margins
# boxplot(select(dd, aromatic:s_lignin),horizontal = TRUE,las=2)
# meanvar.plot(mv_spp)
# mod1 <- manyglm(mv_spp ~ d13C, data = dd, family = "gamma")
# plot(mod1)
# plot(mod1, which = 2)
# anova(mod1)
# anova(mod1, p.uni="adjusted")
# 

# lignin:carbohydrate ratios



# G, S, Phenol lignin analysis --------------------------------------------

# G lignin
ggplot(dd, aes(x = log(FC), y = g_lignin, col = Site))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_grid(Vegetation ~ Layer)

glig_m1 <- lm(g_lignin ~ log(FC) * Vegetation * Layer, dd)
summary(glig_m1)


glig_ps_m1 <- lm(g_lignin ~ log(FC) * Layer, filter(dd, Vegetation == "Pynus sylvestris"))
Anova(glig_ps_m1)
summary(glig_ps_m1)


glig_pa_m1 <- lm(g_lignin ~ log(FC) * Layer, filter(dd, Vegetation == "Picea abies"))
summary(glig_pa_m1)


# S lignin
ggplot(dd, aes(x = log(FC), y = log(s_lignin), col = Site))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_grid(Vegetation ~ Layer)

slig_m1 <- lm(log(s_lignin) ~ log(FC) * Vegetation * Layer, dd)
summary(slig_m1)
Anova(slig_m1)
visreg(slig_m1, xvar = "FC")


# Phenol
ggplot(dd, aes(x = log(FC), y = Phenol, col = Site))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_grid(Vegetation ~ Layer)

plig_m1 <- lm(Phenol ~ log(FC) * Vegetation * Layer, dd)
summary(plig_m1)
Anova(plig_m1)

plig_ps_m1 <- lm(Phenol ~ log(FC) * Layer, filter(dd, Vegetation == "Pynus sylvestris"))
summary(plig_ps_m1)
Anova(plig_ps_m1)
visreg(plig_ps_m1, xvar = "FC", by = "Layer", overlay = TRUE)

plig_pa_m1 <- lm(Phenol ~ log(FC) * Layer, filter(dd, Vegetation == "Picea abies"))
summary(plig_pa_m1)
Anova(plig_pa_m1)
visreg(plig_pa_m1, xvar = "FC", by = "Layer", overlay = TRUE)
