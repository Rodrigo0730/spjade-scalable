# Benchmark: iterative kernel refinement — geographic vs data-driven
# Methods:
#   - classical spJADE (geographic ring, α=1, no iteration)
#   - iterative α=0.5 (balanced blend)
#   - iterative α=0.2 (mostly source-driven)
#   - iterative α=0.0 (pure source-space kernel, no geography)
# Setting 3, d = 5, 10, 20, 10 reps each, max 5 iterations

library(Matrix); library(spdep); library(JADE); library(SpatialBSS); library(sp); library(sf)
source("R/datagen.R"); source("R/analysis.R")

SETTING   <- "setting_3"
DS_VALUES <- c(5, 10, 20)
N_REPS    <- 20
BD        <- c(0, 1, 1, 2, 2, 3)
SEED_BASE <- 2025
MAX_ITER  <- 5
ALPHAS    <- c(1.0, 0.5, 0.2, 0.0)

cat("Legend: α=1 (geo only), α=0.5 (blended), α=0.2 (mostly src), α=0 (src only)\n")
cat(sprintf("%-6s %-3s %8s %8s %8s %8s %s\n", "d", "rep", "cl-geo", "α=0.5", "α=0.2", "α=0.0", "iters"))

all_results <- list()

for (d in DS_VALUES) {
  cat(sprintf("\n=== d = %d ===\n", d))
  filename <- file.path("data", SETTING, sprintf("data_%d.rds", d))
  data   <- readRDS(filename)
  field  <- gen_field(d)
  rings  <- gen_rings(field, BD, row_standardize = FALSE, symmetrize = TRUE, f0 = TRUE)

  kernels        <- rings$kernels[-1]
  kernels_sparse <- rings$kernels_sparse[-1]
  Sf_list        <- lapply(kernels, function(kmat) sum(kmat) / nrow(kmat))

  results <- vector("list", N_REPS)

  for (r in 1:N_REPS) {
    sources <- as.matrix(data[[r]])
    p <- ncol(sources)
    set.seed(SEED_BASE + r)
    A <- matrix(rnorm(p * p), p, p)
    X <- tcrossprod(sources, A)

    # Classical spJADE (geographic only)
    md_geo <- spJADE(X, A, kernels, kernels_sparse, Sf_list, norm_scheme = "none")$md

    mds <- c(md_geo)
    niters <- c(0)

    for (alpha in ALPHAS[2:4]) {  # α=0.5, 0.2, 0.0
      res <- iterative_spjade(X, A, kernels, kernels_sparse, Sf_list,
        max_iter = MAX_ITER, alpha = alpha, sigma_scale = 1, tol = 1e-5)
      mds    <- c(mds, res$best_md)
      niters <- c(niters, res$n_iter)
    }

    results[[r]] <- data.frame(
      d = d, rep = r,
      md_geo = mds[1], md_alpha05 = mds[2],
      md_alpha02 = mds[3], md_alpha00 = mds[4],
      iters_05 = niters[2], iters_02 = niters[3], iters_00 = niters[4],
      stringsAsFactors = FALSE
    )

    cat(sprintf("%-6d %-3d %8.4f %8.4f %8.4f %8.4f %s\n",
      d, r, mds[1], mds[2], mds[3], mds[4],
      paste(niters[2:4], collapse = "/")))
  }

  df <- do.call(rbind, results)

  cat(sprintf("\n  geo (α=1):    mean=%.4f med=%.4f sd=%.4f\n",
    mean(df$md_geo), median(df$md_geo), sd(df$md_geo)))
  cat(sprintf("  blend (α=0.5): mean=%.4f med=%.4f sd=%.4f  Δ=%.4f  iters=%.1f\n",
    mean(df$md_alpha05), median(df$md_alpha05), sd(df$md_alpha05),
    mean(df$md_alpha05 - df$md_geo), mean(df$iters_05)))
  cat(sprintf("  mostly (α=0.2):mean=%.4f med=%.4f sd=%.4f  Δ=%.4f  iters=%.1f\n",
    mean(df$md_alpha02), median(df$md_alpha02), sd(df$md_alpha02),
    mean(df$md_alpha02 - df$md_geo), mean(df$iters_02)))
  cat(sprintf("  src-only (α=0):mean=%.4f med=%.4f sd=%.4f  Δ=%.4f  iters=%.1f\n",
    mean(df$md_alpha00), median(df$md_alpha00), sd(df$md_alpha00),
    mean(df$md_alpha00 - df$md_geo), mean(df$iters_00)))

  all_results[[as.character(d)]] <- df
}

saveRDS(all_results, "results/iterative_benchmark.rds")
cat("\nDone.\n")
