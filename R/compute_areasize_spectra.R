
# Here, I computed the relative abundance of area size of raw spectra. This will
# be helpful when assigning compounds using NIST.

source("R/packages.R")

# Area size of raw spectrum for NIST --------------------------------------

humus_spec_area_raw <- llply(list('Aheden'      = "Data/Spectra_Humus_Aheden.csv", 
                                  'Svartberget' = "Data/Spectra_Humus_Svartberget.csv",
                                  'Flakaliden'  = "Data/Spectra_Humus_Flakaliden.csv",
                                  'RosinedalOF_v2' = "Data/Spectra_Humus_RosinedalOF_v2.csv"),
                             function(x){
                               d0 <- read.csv(x) 
                               d1 <- d0 %>% 
                                 filter(Window == "Win001_C01") %>% 
                                 mutate(Ave_Area_Prop = 0) %>% 
                                 select(Window, RT_s, Ave_Area_Prop)
                               d2 <- d0 %>% 
                                 transmute(Window, RT_s, tot_area = rowSums(select(., starts_with("X")))) %>% 
                                 filter(Window != "Win001_C01") %>%  # remove CO2
                                 mutate(Ave_Area_Prop = tot_area * 100 / sum(tot_area)) %>% 
                                 select(Window, RT_s, Ave_Area_Prop)
                               return(rbind(d1, d2))
                             })
l_ply(names(humus_spec_area_raw),
      function(x) write.csv(humus_spec_area_raw[[x]], 
                            paste0("Output/Data/", x, "_Spectrum_area_humus.csv"),
                            row.names = FALSE))




# US sample ---------------------------------------------------------------

us_spec_area_raw <- llply(list('BR'      = "Data/US_Soil/Spectra/Spectra_BR.csv", 
                               'CA'      = "Data/US_Soil/Spectra/Spectra_CA.csv", 
                               'FE'      = "Data/US_Soil/Spectra/Spectra_FE.csv", 
                               'HF'      = "Data/US_Soil/Spectra/Spectra_HF.csv", 
                               'KA'      = "Data/US_Soil/Spectra/Spectra_KA.csv", 
                               'ME'      = "Data/US_Soil/Spectra/Spectra_ME.csv", 
                               'NH'      = "Data/US_Soil/Spectra/Spectra_NH.csv"
                               ),
                             function(x){
                               d0 <- read.csv(x) 
                               d1 <- d0 %>% 
                                 filter(Window == "Win001_C01") %>% 
                                 mutate(Ave_Area_Prop = 0) %>% 
                                 select(Window, RT_s, Ave_Area_Prop)
                               d2 <- d0 %>% 
                                 transmute(Window, RT_s, tot_area = rowSums(.[, -c(1:3)])) %>%  
                                 filter(Window != "Win001_C01") %>%  # remove CO2
                                 mutate(Ave_Area_Prop = tot_area * 100 / sum(tot_area)) %>% 
                                 select(Window, RT_s, Ave_Area_Prop)
                               return(rbind(d1, d2))
                             })
l_ply(names(us_spec_area_raw),
      function(x) write.csv(us_spec_area_raw[[x]], 
                            paste0("Output/Data/US_soil/", x, "_Spectrum_area.csv"),
                            row.names = TRUE))

