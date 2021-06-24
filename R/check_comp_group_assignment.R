# Here, I will check compounds and assigned groups for the US samples

source("R/packages.R")
library(utils)

# read files
comp_US_fls <- dir("Data/NIST_compound_comparison/", "Compound_US.*")
us_comps <- llply(comp_US_fls, function(x){
  d <- read.csv(paste0("Data/NIST_compound_comparison/", x)) %>% 
    select(Comp, Group) %>% 
    filter(Group != "unk")
  site <- strsplit(x, "_|[.]")[[1]][3]
  colnames(d)[2] <- paste0(site, "_Group")
  return(d)
})

# merge by compound names
us_merged <- Reduce(function(x, y) merge(x, y, by = "Comp", all = TRUE), us_comps)

# found the rows with discrepancy in the assigned groups between sites
dscrpcy <- apply(us_merged, 1, function(x){
  d1 <- as.matrix(x)[-1]  # turn into a vector
  d2 <- d1[!is.na(d1)]    # remove NA
  l <- ifelse(length(d2) == 1, TRUE, all(combn(d2, 2, function(x) identical(x[1], x[2])))) # check if assigned groups are identical for all combination between the sites
  return(l)
})
us_dscrpcy_d <- us_merged[!dscrpcy, ]


