# Speed benchmark: spJADE vs riemann_nn_spjade at d=80
# n = d^2 = 6400 — where O(n^2) vs O(n) scaling should matter

library(torch); library(JADE); library(SpatialBSS); library(spdep); library(Matrix)
source("R/datagen.R"); source("R/analysis.R")

d <- 80; n_reps <- 10; bd <- c(0, 1, 1, 2, 2, 3)

cat(sprintf("Loading d=%d (n=%d) data...\n", d, d^2))
data  <- readRDS(sprintf("data/setting_3/data_%d.rds", d))
field <- gen_field(d)
rings <- gen_rings(field, bd, row_standardize = FALSE, symmetrize = TRUE, f0 = TRUE)
kernels        <- rings$kernels[-1]
kernels_sparse <- rings$kernels_sparse[-1]
Sf_list        <- lapply(kernels, function(kmat) sum(kmat) / nrow(kmat))

times_sp <- numeric(n_reps)
times_ri <- numeric(n_reps)

for (r in 1:n_reps) {
  sources <- as.matrix(data[[r]])
  p <- ncol(sources)
  set.seed(2025 + r)
  A <- matrix(rnorm(p * p), p, p)
  X <- tcrossprod(sources, A)

  t0 <- proc.time()
  res_sp <- spJADE(X, A, kernels, kernels_sparse, Sf_list)$md
  t_sp  <- (proc.time() - t0)["elapsed"]

  t0 <- proc.time()
  res_ri <- riemann_nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
    lr = 0.01, maxiter = 100, seed = 2025 + r)$md
  t_ri  <- (proc.time() - t0)["elapsed"]

  times_sp[r] <- t_sp
  times_ri[r] <- t_ri

  cat(sprintf("  rep %2d | spJADE MD=%.4f  t=%.2fs | riemann MD=%.4f  t=%.2fs | ratio=%.1fx\n",
    r, res_sp, t_sp, res_ri, t_ri, t_ri / t_sp))
}

cat(sprintf("\n--- Summary (%d reps) ---\n", n_reps))
cat(sprintf("  spJADE:    med=%.2fs  mean=%.2fs  min=%.2fs  max=%.2fs\n",
  median(times_sp), mean(times_sp), min(times_sp), max(times_sp)))
cat(sprintf("  riemann:   med=%.2fs  mean=%.2fs  min=%.2fs  max=%.2fs\n",
  median(times_ri), mean(times_ri), min(times_ri), max(times_ri)))
cat(sprintf("  riemann/spJADE ratio: %.1fx\n", median(times_ri) / median(times_sp)))
