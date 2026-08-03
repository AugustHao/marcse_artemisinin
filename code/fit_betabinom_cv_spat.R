# WRAP MODEL FITTING WITH SPATIALLY STRATIFIED CROSS VALIDATION

library(parallel)
library(caret)
library(parallelly)
library(blockCV)
print(paste("Cores:", detectCores()))
print(paste("Cores:", availableCores()))

suppressMessages(source("code/setup.R"))
suppressMessages(source("code/build_design_matrix.R"))
suppressMessages(source("code/betabinomial_p_rho.R"))
suppressMessages(source("code/wrap_fit.R"))

args <- commandArgs(trailingOnly = TRUE)
marker <- args[1]  # of "k13_marcse", "mdr86", "mdr184", "mdr1246", "crt76" etc
seed <- as.numeric(args[2])
print(paste0("Marker: ", marker))
print(paste0("Seed: ", seed))

set.seed(seed)

out_dir <- paste0(marker, "/bb_gne/")

in_dat <- data_path_lookup[[marker]]


print(paste0("Reading in from: ", in_dat))
print("Enforcing min year for surveyor data - 2000")
print("Fitting betabinom")

mut_data <- setup_mut_data(in_dat, 
                           min_year = MIN_YEAR, 
                           buffer = BUFFER)
write_rds(mut_data, paste0("output/", out_dir, "mut_data.rds"))

NFOLD <- 10

# pop mut_data into sf format
mut_data_sf <- st_as_sf(mut_data, coords = c("x", "y"), crs = st_crs(afr))

# from blockCV docs - seed set above
folds <- cv_spatial(
  x = mut_data_sf,
  k = NFOLD, # number of folds
  size = 500000, # size of the blocks in metres
  selection = "random", # random blocks-to-fold
  iteration = 50, # find evenly dispersed folds
  progress = FALSE
)

# for supp?:
mut_data_sf$folds <- folds$folds_ids

p <- ggplot() +
  geom_sf(data = afr, fill = NA) +
  geom_sf(data = folds$blocks) +
  geom_sf(data = mut_data_sf, pch = 1) +
  facet_wrap(~folds)

ggsave(p, paste("figures/spat_folds_", marker, ".png"))

# write in same format as caret::createFolds
folds <- lapply(folds$folds_list, function(x){
  x[[2]]
})

names(folds) <- paste0("Fold", 1:NFOLD)

write_rds(folds, paste0("output/", out_dir, "cv_folds_spat.rds"))

system.time(mclapply(1:NFOLD, function(x){
  fit_betabinom(mut_data = mut_data,
                covariates = covariates,
                pfpr_years = pfpr_years,
                out_dir = out_dir,
                fold = x,
                folds = folds,
                nchains = 6,
                warmup = 5000,
                nsamples = 30000)
}, mc.cores = NFOLD))