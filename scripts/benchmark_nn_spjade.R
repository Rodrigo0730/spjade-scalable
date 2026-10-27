# Full benchmark: nn_spjade vs spJADE — exact equivalence
# Uses the vectorized Adam JD approach (no C_f^{kl} construction)

pkg_list <- c("Matrix", "spdep", "dplyr", "SpatialBSS", "JADE", "sp", "sf", "torch")
for (pkg in pkg_list) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

source("R/datagen.R")
source("R/analysis.R")

SETTING        <- "setting_3"
DS_VALUES      <- c(5, 10, 20)
N_REPS         <- 100
BD             <- c(0, 1, 1, 2, 2, 3)
SEED_BASE      <- 2025

OUTPUT_DIR <- file.path("results", "nn_benchmark", "raw", SETTING)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (d in DS_VALUES) {

  cat(sprintf("\n========== d = %d | %d reps ==========\n", d, N_REPS))
  cat(sprintf("[%s] Loading data...\n", Sys.time()))

  filename  <- file.path("data", SETTING, sprintf("data_%d.rds", d))
  data      <- readRDS(filename)
  field     <- gen_field(d)

  rings_all <- gen_rings(field, BD,
    row_standardize = FALSE, symmetrize = TRUE, f0 = TRUE
  )

  kernels        <- rings_all$kernels[-1]
  kernels_sparse <- rings_all$kernels_sparse[-1]
  Sf_list        <- lapply(kernels, function(kmat) sum(kmat) / nrow(kmat))

  results <- vector("list", N_REPS)

  for (r in 1:N_REPS) {
    sources <- as.matrix(data[[r]])
    p       <- ncol(sources)

    set.seed(SEED_BASE + r)
    A <- matrix(rnorm(p * p), p, p)
    X <- tcrossprod(sources, A)

    t0_sp  <- proc.time()
    res_sp <- spJADE(X, A, kernels, kernels_sparse, Sf_list, norm_scheme = "none")
    time_sp <- (proc.time() - t0_sp)["elapsed"]

    t0_nn  <- proc.time()
    res_nn <- nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
      lr = 0.01, maxiter = 500, lambda_orth = 100, seed = SEED_BASE + r
    )
    time_nn <- (proc.time() - t0_nn)["elapsed"]

    t0_j  <- proc.time()
    md_j  <- jade_md(X, A)
    time_j <- (proc.time() - t0_j)["elapsed"]

    results[[r]] <- data.frame(
      rep            = r, d = d,
      md_spjade      = res_sp$md,
      md_nn          = res_nn$md,
      md_jade        = md_j,
      time_spjade    = time_sp,
      time_nn        = time_nn,
      time_jade      = time_j,
      nn_final_loss  = tail(res_nn$loss_history, 1),
      nn_loss_min    = min(res_nn$loss_history),
      stringsAsFactors = FALSE
    )

    if (r %% 10 == 0) {
      cat(sprintf("  rep %3d | md_sp=%.4f  md_nn=%.4f | Δ=%.4f | t_sp=%.2fs  t_nn=%.2fs\n",
        r, res_sp$md, res_nn$md, res_nn$md - res_sp$md, time_sp, time_nn))
    }
  }

  df <- do.call(rbind, results)
  saveRDS(df, file.path(OUTPUT_DIR, sprintf("res_%d.rds", d)))
  cat(sprintf("[%s] Saved. ", Sys.time()))

  cat(sprintf(
    "spJADE MD: mean=%.4f med=%.4f sd=%.4f t_med=%.2fs | ",
    mean(df$md_spjade, na.rm = TRUE), median(df$md_spjade, na.rm = TRUE),
    sd(df$md_spjade, na.rm = TRUE), median(df$time_spjade, na.rm = TRUE)
  ))
  cat(sprintf(
    "nn MD: mean=%.4f med=%.4f sd=%.4f t_med=%.2fs | Δmean=%.4f\n",
    mean(df$md_nn, na.rm = TRUE), median(df$md_nn, na.rm = TRUE),
    sd(df$md_nn, na.rm = TRUE), median(df$time_nn, na.rm = TRUE),
    mean(df$md_nn - df$md_spjade, na.rm = TRUE)
  ))
}

cat(sprintf("\n[%s] Done.\n", Sys.time()))
