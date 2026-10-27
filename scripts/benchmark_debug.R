# Debug: nn_spjade vs spJADE — exact equivalence test (3 reps, d=5)

library(torch)
library(JADE)
library(SpatialBSS)
library(spdep)
library(Matrix)
library(dplyr)

source("R/datagen.R")
source("R/analysis.R")

d <- 5
n_reps <- 3
bd <- c(0, 1, 1, 2, 2, 3)

filename <- file.path("data", "setting_3", sprintf("data_%d.rds", d))
data <- readRDS(filename)
field <- gen_field(d)

rings_all <- gen_rings(field, bd,
  row_standardize = FALSE, symmetrize = TRUE, f0 = TRUE
)
kernels        <- rings_all$kernels[-1]
kernels_sparse <- rings_all$kernels_sparse[-1]
Sf_list <- lapply(kernels, function(kmat) sum(kmat) / nrow(kmat))

for (r in 1:n_reps) {
  cat(sprintf("\n========== Rep %d ==========\n", r))

  sources <- as.matrix(data[[r]])
  p <- ncol(sources)

  set.seed(2025 + r)
  A <- matrix(rnorm(p * p), p, p)
  X <- tcrossprod(sources, A)

  t0 <- proc.time()
  md_sp <- spJADE(X, A, kernels, kernels_sparse, Sf_list)$md
  t_sp <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  spJADE     MD=%.6f  t=%.3fs\n", md_sp, t_sp))

  t0 <- proc.time()
  res_nn <- nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
    lr = 0.01, maxiter = 500, lambda_orth = 100, seed = 2025 + r
  )
  t_nn <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  nn_spjade  MD=%.6f  t=%.3fs\n", res_nn$md, t_nn))
  cat(sprintf("             loss: start=%.4f  end=%.4f  min=%.4f\n",
    res_nn$loss_history[1],
    tail(res_nn$loss_history, 1),
    min(res_nn$loss_history)))

  rel_change <- (tail(res_nn$loss_history, 1) - res_nn$loss_history[1]) / abs(res_nn$loss_history[1])
  cat(sprintf("             rel change: %.4f%%\n", rel_change * 100))
  cat(sprintf("  nn vs sp:  ΔMD=%.6f\n", res_nn$md - md_sp))
}

cat("\nDone.\n")
