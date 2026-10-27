# Comprehensive benchmark: all nn_spjade variants vs classical spJADE
# Methods: spJADE, nn_spjade, implicit, riemann, rff, hybrid
# Setting 3 (log-spGARCH), d = 5, 10, 20, 40, 60, 50 reps each

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
DS_VALUES      <- c(5, 10, 20, 40, 60)
BD             <- c(0, 1, 1, 2, 2, 3)
SEED_BASE      <- 2025

OUTPUT_DIR <- file.path("results", "benchmark_all", "raw", SETTING)
DIR_PLOTS  <- file.path("results", "benchmark_all", "plots")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_PLOTS, recursive = TRUE, showWarnings = FALSE)

available_d <- DS_VALUES[sapply(DS_VALUES, function(d) {
  file.exists(file.path("data", SETTING, sprintf("data_%d.rds", d)))
})]
cat("Available d values:", paste(available_d, collapse = ", "), "\n")

rep_counts <- setNames(c(100, 100, 100, 50, 50), c("5", "10", "20", "40", "60"))

cat("\n========== Method legend ==========\n")
cat("  M0: spJADE          (classical C_f^{kl} + frjd)\n")
cat("  M1: nn_spjade       (exact JD, 500 Adam epochs)\n")
cat("  M2: implicit        (diagonal auto-kurtosis, 500 Adam)\n")
cat("  M3: riemann         (Stiefel manifold Adam, 100 epochs)\n")
cat("  M4: rff             (RFF kernel approx + exact JD, 500 Adam)\n")
cat("  M5: hybrid          (implicit 50 → exact JD 200)\n")
cat("==================================\n\n")

all_results <- list()

for (d in available_d) {

  N_REPS <- rep_counts[as.character(d)]
  if (is.na(N_REPS)) N_REPS <- 100

  cat(sprintf("\n########## d = %d | %d reps ##########\n", d, N_REPS))
  cat(sprintf("[%s] Loading data...\n", Sys.time()))

  filename  <- file.path("data", SETTING, sprintf("data_%d.rds", d))
  data      <- readRDS(filename)
  field     <- gen_field(d)

  rings <- gen_rings(field, BD, row_standardize = FALSE, symmetrize = TRUE, f0 = TRUE)
  kernels        <- rings$kernels[-1]
  kernels_sparse <- rings$kernels_sparse[-1]
  Sf_list        <- lapply(kernels, function(kmat) sum(kmat) / nrow(kmat))

  results <- vector("list", N_REPS)

  for (r in 1:N_REPS) {
    sources <- as.matrix(data[[r]])
    p       <- ncol(sources)

    set.seed(SEED_BASE + r)
    A <- matrix(rnorm(p * p), p, p)
    X <- tcrossprod(sources, A)

    # ---- M0: classical spJADE ----
    t0 <- proc.time()
    m0 <- spJADE(X, A, kernels, kernels_sparse, Sf_list, norm_scheme = "none")$md
    t0 <- (proc.time() - t0)["elapsed"]

    # ---- M1: exact JD nn_spjade ----
    t1 <- proc.time()
    m1 <- nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
      lr = 0.01, maxiter = 500, lambda_orth = 100, seed = SEED_BASE + r
    )$md
    t1 <- (proc.time() - t1)["elapsed"]

    # ---- M2: implicit diagonal-only ----
    t2 <- proc.time()
    m2 <- implicit_nn_spjade(X, A, kernels_sparse,
      lr = 0.01, maxiter = 500, lambda_orth = 100, seed = SEED_BASE + r
    )$md
    t2 <- (proc.time() - t2)["elapsed"]

    # ---- M3: riemann Stiefel Adam (100 epochs) ----
    t3 <- proc.time()
    m3 <- riemann_nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
      lr = 0.01, maxiter = 100, seed = SEED_BASE + r
    )$md
    t3 <- (proc.time() - t3)["elapsed"]

    # ---- M4: RFF kernel approximation ----
    t4 <- proc.time()
    m4 <- rff_nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
      D = 50, lr = 0.01, maxiter = 500, lambda_orth = 100, seed = SEED_BASE + r
    )$md
    t4 <- (proc.time() - t4)["elapsed"]

    # ---- M5: hybrid implicit → JD ----
    t5 <- proc.time()
    m5 <- hybrid_nn_spjade(X, A, kernels, kernels_sparse, Sf_list,
      lr = 0.01, maxiter_implicit = 50, maxiter_jd = 200,
      lambda_orth = 100, seed = SEED_BASE + r
    )$md
    t5 <- (proc.time() - t5)["elapsed"]

    results[[r]] <- data.frame(
      rep     = r, d = d,
      md_spjade   = m0, time_spjade   = t0,
      md_nn       = m1, time_nn       = t1,
      md_implicit = m2, time_implicit = t2,
      md_riemann  = m3, time_riemann  = t3,
      md_rff      = m4, time_rff      = t4,
      md_hybrid   = m5, time_hybrid   = t5,
      stringsAsFactors = FALSE
    )

    if (r %% 10 == 0 || r == 1) {
      cat(sprintf("  rep %3d | M0=%.4f M1=%.4f M2=%.4f M3=%.4f M4=%.4f M5=%.4f | t %s\n",
        r, m0, m1, m2, m3, m4, m5,
        paste(sprintf("%.1f", c(t0, t1, t2, t3, t4, t5)), collapse = "/")))
    }
  }

  df <- do.call(rbind, results)
  out_file <- file.path(OUTPUT_DIR, sprintf("res_%d.rds", d))
  saveRDS(df, out_file)

  # Per-method summaries
  for (mname in c("spjade", "nn", "implicit", "riemann", "rff", "hybrid")) {
    md_col <- paste0("md_", mname)
    t_col  <- paste0("time_", mname)
    md_vals <- df[[md_col]]
    t_vals  <- df[[t_col]]

    cat(sprintf("  %-10s MD: mean=%.4f med=%.4f sd=%.4f  ΔM0=%.4f  t_med=%.2fs\n",
      mname,
      mean(md_vals, na.rm = TRUE),
      median(md_vals, na.rm = TRUE),
      sd(md_vals, na.rm = TRUE),
      mean(md_vals - df$md_spjade, na.rm = TRUE),
      median(t_vals, na.rm = TRUE)
    ))
  }

  all_results[[as.character(d)]] <- df
}

saveRDS(all_results, file.path(OUTPUT_DIR, "all_results.rds"))
cat(sprintf("\n[%s] Benchmark complete.\n", Sys.time()))
