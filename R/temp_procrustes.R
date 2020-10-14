install.packages("bigutilsr")
library(bigutilsr)
A <- matrix(rnorm(200), ncol = 20)
B <- matrix(rnorm(length(A)), nrow = nrow(A))
proc <- procrustes(B, A)
str(proc)
plot(B, predict(proc, A)); abline(0, 1, col = "red")


protest(B, A)
ff1 <- vegan::procrustes(B, A)

plot(B, predict(proc, A)); abline(0, 1, col = "red")
plot(predict(proc, A), predict(ff1, A)); abline(0, 1, col = "red")
summary(vegan::procrustes(B, A))



data(varespec)
vare.dist <- vegdist(wisconsin(varespec))
mds.null <- monoMDS(vare.dist, y = cmdscale(vare.dist))
mds.alt <- monoMDS(vare.dist)
vare.proc <- vegan::procrustes(mds.alt, mds.null)
vare.proc
protest(vare.proc)
summary(vare.proc)
plot(vare.proc)
plot(vare.proc, kind=2)
residuals(vare.proc)



# dd_cntr <- filter(dd, treatment == "Control")
dd_cntr <- dd
dd_cntr_sp <- decostand(select(dd_cntr, aromatic:s_lignin), method = "hellinger")
dd_cntr_pca <- rda(dd_cntr_sp ~ 1)
summary(dd_cntr_pca)
dd_cntr_env <- dd_cntr %>% 
  select(wN, wC, FC)
  # select(wN, wC, FC, Vegetation, Layer)
# %>%
  # mutate(Vegetation = ifelse(Vegetation == "Pynus sylvestris", 1, 0),
         # Layer = ifelse(Layer == "Litter", 1, 0))
dd_cntr_env_tr <- decostand(dd_cntr_env, method = "hellinger")
  
a1 <- scores(metaMDS(dd_cntr_sp, distance = "bray"), display = "sites")
a2 <- scores(metaMDS(dd_cntr_env_tr, distance = "bray") , display = "sites")

a1 <- monoMDS(vegdist(dd_cntr_sp, "bray"))
a2 <- monoMDS(vegdist(dd_cntr_env_tr, "bray"))

pp <- vegan::procrustes(a1, a2)
summary(pp)
plot(pp)
plot(pp, kind = 2)
protest(a1, a2)

newd <- cbind(dd, resid = residuals(pp))
ggplot(newd, aes(x = Trt, y = resid, col = Site))+
  geom_boxplot()+
  facet_grid(Layer ~ Site, scale = "free_x", space = "free")


proc <- bigutilsr::procrustes(a1, a2)
proc
plot(a1, predict(proc, a2)); abline(0, 1, col = "red")
