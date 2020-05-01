# create axis labels from RDA obj for PCA plot ----------------------------

get_PCA_axislab <- function(pcaobj){
  d        <- summary(pcaobj)$cont$importance  # get eigen values
  pro_expl <- round(d[2, ] * 100, 1)
  axes     <- paste0(names(pro_expl), " (", pro_expl, "%)")
  return(axes)
}

