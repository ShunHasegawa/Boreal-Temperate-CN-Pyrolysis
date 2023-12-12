
# Here, I align the lists of compounds by retention time identified by NIST from
# multiple sampling sites to check if assigned compounds are reasonalble relative
# to the other sites

source("R/packages.R")

comp_fls <- dir("Data/NIST_compound_comparison/")

all_comps <- llply(comp_fls, function(x){
  site <- strsplit(x, "_|[.]")[[1]][3]
  type <- strsplit(x, "_|[.]")[[1]][2]
  d <- read.csv(paste0("Data/NIST_compound_comparison/", x)) %>% 
    mutate(RT_s = round(RT_s, 0)) %>%  # round RT_s
    select(Window, RT_s, Comp) %>% 
    group_by(RT_s) %>%
    mutate(rep = 1:n())                # Some RT_s have multiple values so add rep number for each 
  colnames(d) <- c(paste(type, site, "Window", sep = "_"),
                    "RT_s",
                    paste(type, site, "Comp", sep = "_"),
                   "rep")
  return(d)
  
})

merged_comp <- Reduce(function(...) merge(..., all = TRUE), all_comps) %>% 
  select(-contains("rep")) %>% 
  mutate_all(.funs = list(~ifelse(is.na(.), "", as.character(.))))

# check if there is no redundant rows
sapply(all_comps, nrow)
merged_comp %>% 
  select(contains("Window")) %>% 
  apply(., 2, function(x) sum(x != ""))

# save
write.csv(merged_comp, "Output/Tables/Compare_NIST_assigned_compds.csv", row.names = FALSE)
