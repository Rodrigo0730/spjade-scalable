# -----------------------------------------------------------------------------
# Spatial BSS MPI Data Generation Framework
# -----------------------------------------------------------------------------
#
# In Slurm choose the setting to run and values of the parameter d. Note that
# the total number of data points is d^2. Add the following after module loads:
#
# SETTING="setting_1/2/3/4"
# DS_VALUES="40,80,120"
#
# Then add to the srun of this R script the following:
#
# --setting=${SETTING} --ds=${DS_VALUES}
#
# See /scripts/example_mpi.sh for and example. This script generates the data
# for a given setting and values of d. The data is saved in 
# /data/setting_X/data_d.rds for setting number "X" and d value "d".
#
# Requirements:
#   - Relies on helper data generation functions sourced from "R/datagen.R".
#
# -----------------------------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(name, default = NULL) {
  val <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(val) == 0) return(default)
  sub(paste0("^--", name, "="), "", val)
}

setting <- parse_arg("setting", default = "setting_1")
ds_input <- parse_arg("ds", default = "5")
ds <- as.integer(strsplit(ds_input, ",")[[1]])

message("Setting: ", setting)
message("Values of parameter d: ", paste(ds, collapse = ", "))

master <- c("doRNG", "Matrix", "doMPI", "spdep", "doParallel", "SpatialBSS")
workers <- c("SpatialBSS","JADE","spGARCH","spdep","sp",
             "dplyr","sf","moments", "Matrix")

#.libPaths(c("/projappl/project_2012081/project_rpackages_440", .libPaths()))

#load packages for master script
for (pkg in master) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg, character.only = TRUE, quietly=TRUE)
  }
}

#load helper functions
source("~/spBSS/R/datagen.R")

message("[", Sys.time(), "] Starting MPI cluster…")

#start MPI cluster
cl <- makeCluster(8)
registerDoParallel(cl)

registerDoRNG(seed=123)

message("[", Sys.time(), "] Cluster started and RNG registered.")

n_reps <- 2000

for (d in ds) {
  
  message("[", Sys.time(), "] Beginning loop for d = ", d)
  
  coords <- gen_field(d)
  nb <- gen_nb(d, type = "rook")
  W1 <- nb2W(nb, style = "W")  
  W2 <- nth_W(nb, k = 2)
  W3 <- nth_W(nb, k = 3)
  
  gen_sources_fun <- switch(setting,
                            setting_1 = gen_sources_setting_1,
                            setting_2 = gen_sources_setting_2,
                            setting_3 = gen_sources_setting_3,
                            setting_4 = gen_sources_setting_4,
                            stop("Unknown setting: ", setting)
  )
  
  message("[", Sys.time(), "] Starting foreach for ", n_reps, " reps…")
  

  results_d <- foreach(
    i = 1:n_reps,
    .combine = "c",
    .packages = workers
  ) %dopar% {
    list(gen_sources_fun(d, W1, W2, W3))
  }

  message("[", Sys.time(), "] Saving results to data_", d, ".rds …")
  
  save_filename <- sprintf("~/spBSS/data/%s/data_%d.rds", setting, d)
  
  saveRDS(results_d, file = save_filename, compress = "gzip")
  
  message("[", Sys.time(), "] Saved results. Clearing memory…")
  
  rm(results_d); gc()
  
  message("[", Sys.time(), "] Cleanup done.")
}


message("[", Sys.time(), "] About to closeCluster(cl)…")
stopCluster(cl)
message("[", Sys.time(), "] closeCluster() returned.")

