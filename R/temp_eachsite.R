
# litter ------------------------------------------------------------------


# Aheden
litterall_ah <- filter(litterall_raw, Site == "Aheden")
litterall_ah <- filter(litterall_raw, Site == "Flakaliden")
litterall_ah <- filter(litterall_raw, Site == "Rosinedal")
litterall_ah <- filter(litterall_raw, Site == "Svartberget")


litterall_ah_sp <- decostand(select(litterall_ah, starts_with("Win")), method = "hellinger")
pcoa_litter_ah <- rda(litterall_ah_sp ~ CN, litterall_ah)
anova(pcoa_litter_ah, by = "margin")
summary(pcoa_litter_ah)
plot(pcoa_litter_ah)


# humus -------------------------------------------------------------------


# Aheden
humusall_ah <- filter(humusall_raw, Site == "Aheden")
humusall_ah <- filter(humusall_raw, Site == "Flakaliden")

humusall_ah <- filter(humusall_raw, Site == "Rosinedal")
humusall_ah <- filter(humusall_raw, Site == "Svartberget")


humusall_ah_sp <- decostand(select(humusall_ah, starts_with("Win")), method = "hellinger")
pcoa_humus_ah <- rda(humusall_ah_sp ~ CN, humusall_ah)
anova(pcoa_humus_ah, by = "margin")
summary(pcoa_humus_ah)
plot(pcoa_humus_ah)

