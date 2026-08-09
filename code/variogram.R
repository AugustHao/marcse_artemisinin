setwd("~/Desktop/MARCSE/k13_seafrica")

source("code/setup.R")

library(gstat)



for (marker in nice_name_lookup_all){
  
}


marker = "k13_marcse"
preds <- rast(paste0("output/", marker, "/bb_gne/preds_medians.tif"))
preds_df <- as.data.frame(preds, xy = TRUE, na.rm = TRUE)
# bug: gstat::variogram doesn't function with column names beginning with numbers
names(preds_df) <- gsub("_50", "", names(preds_df))
preds_df$tmp <- preds_df$`2000_50`

# this may be a job for the cluster as it takes a while ...
v <- gstat::variogram(tmp ~ 1, locations = ~ x + y, data = preds_df)

plot(v)
# show annual variograms over same axes, subpanelled into models?