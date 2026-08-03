
extract_pfpr <- function(df,
                          covs,
                          max_year = 2024,
                          buffer = 0){
  # it's time to throw in the towel chatty g can do this way better than me
  
  yrs_covs <- str_extract(names(covs), "\\d{4}")
  
  # get predictions for each row in `mut_data`
  df$pfpr <- NA
  yrs_to_extract <- unique(df$year)
  for (yr in yrs_to_extract){
    message(yr)
    idx <- which(df$year == yr)
    
    yr <- min(yr, max_year)
    val <- terra::extract(covs[[paste0("pfpr_", yr)]], 
                          df[idx, c("x", "y")],
                          ID = FALSE, 
                          search_radius = buffer) 
    
    df[idx, "pfpr"] <- val[, paste0("pfpr_", yr)]
  }
  
  df
}

pair_of_rasts <- function(ras, year, pfpr_upper_limit = 2024){
  # give me the right pair of pfpr and prediction rasters for a given year
  c(ras[paste0(year, "_50")],
    ras[paste0("pfpr_", min(as.numeric(year), pfpr_upper_limit))]) %>%
    setNames(c("pred", "pfpr"))
}

#' Run-up to response plots
#' Partitioning as it's a big old object and I don't want to run every time
#' I make up the plot
#'
#' @param ras rast of pfprs and preds
#'
#' @returns df of pfprs and preds
#' @export
#'
#' @examples
response_plots_runup <- function(ras){
  years <- ras %>%
    names() %>%
    str_extract(pattern = "\\d{4}")
  
  lapply(years, function(year){
    pair_of_rasts(ras, year) %>%
      as.data.frame() %>%
      mutate(year = as.numeric(year)) %>%
      drop_na()
  }) %>%
    do.call(what = rbind)
}

#' "Response" plots
#' How does prediction vary over covariate ?
#'
#' @param preds rast
#' @param covar rast
#'
#' @returns
#' @export
#'
#' @examples
response_plot <- function(df, dat = NULL, covar = "pfpr", xax_breaks = 100){
  # summarise distribution of all median preds against all covar vals
  # for each value in pfpr, find median and 95% quantile
  # might be more convenient to grab from prediction design matrix ..
  
  xax_bins <- seq(min(df[, covar]), max(df[, covar]), length.out = xax_breaks)
  df$xbinned <- cut(df[, covar], breaks = xax_bins)
  
  to_plot <- df %>% 
    group_by(xbinned) %>%
    summarise(med = median(pred),
              lower2.5 = quantile(pred, probs = c(0.025)),
              lower25 = quantile(pred, probs = c(0.25)),
              upper25 = quantile(pred, probs = c(0.75)),
              upper2.5 = quantile(pred, probs = c(0.975)),
              min = min(pred),
              max = max(pred)) %>%
    mutate(xnume = xax_bins[2:xax_breaks])
  
  p <- ggplot(data = to_plot) +
    geom_line(aes(x = xnume, y = med)) +
    # geom_line(aes(x = xnume, y = min)) +
    # geom_line(aes(x = xnume, y = max)) +
    geom_ribbon(aes(ymin = lower2.5, ymax = upper2.5, x = xnume),
                alpha = 0.2) +
    geom_ribbon(aes(ymin = lower25, ymax = upper25, x = xnume),
                alpha = 0.2) +
    xlab("PFPR (unscaled)") +
    ylab("Median predicted prevalence")

  if (!is.null(dat)){
    # p + geom_point(aes(x = pfpr, y = pred), data = dat)
    p + geom_point(aes(x = pfpr, y = present/tested, size = tested), data = dat,
                   alpha = 0.6, pch = 1) +
      scale_size_continuous(trans = "sqrt", breaks = c(10, 100, 1000, 3000))
  }

}

response_plot(df = response_plot_runup,
              dat = mut_dat_assoc_with_preds$k13_marcse)

ggplot(data = response_plot_runup) +
  geom_density(aes(x = pfpr))

#' Spatial/temporal variogram
#' How does a prediction at one location relate to a prediction at another location ?
#'
#' @param preds 
#'
#' @returns
#' @export
#'
#' @examples
variogram <- function(preds){
  
}

library(terra)

source("code/setup.R")
source("code/build_design_matrix.R")

pfpr_unscaled <- rast("data/pfpr_rasters_afr_2025.tif")
names(pfpr_unscaled) <- paste0("pfpr_", years)

# from validation.R
mut_dat_assoc_with_preds <- lapply(names(nice_name_lookup_all), function(marker){
  extract_preds(data_path = data_path_lookup[[marker]],
                      pred_path = paste0(bb_paths[[marker]], "preds_medians.tif"),
                      buffer = BUFFER) %>%
    extract_pfpr(covs = covariates)
  
}) %>%
  setNames(names(nice_name_lookup_all)) %>%
  suppressMessages()

preds <- rast("output/k13_marcse/bb_gne/preds_medians.tif")
ras <- c(preds, covariates)
response_plot_runup <- response_plots_runup(ras)




