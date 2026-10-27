# Stochastic subsampling benchmark for spJADE
# Tests: d={20,40,80}, m_frac={0.05,0.10,0.20,0.30,0.50,1.0}, n_stoch=5

library(Matrix); library(spdep); library(JADE); library(SpatialBSS); library(sp); library(sf)
source("R/datagen.R"); source("R/analysis.R")

ds      <- c(20, 40, 80)
m_fracs <- c(0.05, 0.10, 0.20, 0.30, 0.50, 1.0)
bd      <- c(0, 1, 1, 2, 2, 3)
n_reps  <- 20
n_stoch <- 5
seed_base <- 2025

cat(sprintf("Subsampling benchmark: d=%s, m_frac=%s, %d reps, %d stoch\n\n",
  paste(ds, collapse="/"), paste(m_fracs, collapse="/"), n_reps, n_stoch))

all_res <- list()

for (d in ds) {
  cat(sprintf("\n===== d = %d (n = %d) =====\n", d, d^2))
  filename <- file.path("data", "setting_3", sprintf("data_%d.rds", d))
  data  <- readRDS(filename)
  field <- gen_field(d)
  
  for (r in 1:n_reps) {
    sources <- as.matrix(data[[r]])
    p <- ncol(sources)
    set.seed(seed_base + r)
    A <- matrix(rnorm(p * p), p, p)
    X <- tcrossprod(sources, A)
    
    # Full spJADE as reference
    rings <- gen_rings(field, bd, row_standardize=FALSE, symmetrize=TRUE, f0=FALSE)
    kerns  <- rings$kernels
    kerns_s <- rings$kernels_sparse
    Sfs    <- lapply(kerns, function(k) sum(k)/nrow(k))
    
    t0 <- proc.time()
    md_full <- spJADE(X, A, kerns, kerns_s, Sfs, norm_scheme="none")$md
    t_full <- (proc.time() - t0)["elapsed"]
    
    row <- data.frame(d = d, rep = r, m_frac = 1.0, 
                      md_sub = md_full, md_sd = 0, time = t_full,
                      stringsAsFactors = FALSE)
    
    for (mf in m_fracs[m_fracs < 1.0]) {
      res <- spJADE_subsample(X, A, field, bd, m_frac = mf, 
                              n_stoch = n_stoch, seed = seed_base + r)
      row <- rbind(row, data.frame(d = d, rep = r, m_frac = mf,
                    md_sub = res$md_mean, md_sd = res$md_sd, 
                    time = res$time_mean,
                    stringsAsFactors = FALSE))
    }
    
    all_res[[paste(d, r, sep="_")]] <- row
    
    if (r %% 5 == 0 || r == 1) {
      sub <- row
      cat(sprintf("  rep %2d | full MD=%.4f t=%.3fs | ", r, md_full, t_full))
      for (i in 1:nrow(sub)) {
        mf <- sub$m_frac[i]
        cat(sprintf("m=%.2f:MD=%.4f ", mf, sub$md_sub[i]))
      }
      cat("\n")
    }
  }
}

df <- do.call(rbind, all_res)
saveRDS(df, "results/subsample_benchmark.rds")

# Summary
cat("\n\n========== SUMMARY ==========\n\n")
cat(sprintf("%-4s %8s %10s %10s %10s %10s %10s\n", 
  "d", "m_frac", "mean_MD", "sd_MD", "ΔMD", "mean_t", "speedup"))
cat(strrep("-", 68), "\n")

for (d in ds) {
  sub <- df[df$d == d, ]
  t_full <- mean(sub$time[sub$m_frac == 1.0])
  md_full <- mean(sub$md_sub[sub$m_frac == 1.0])
  
  for (mf in m_fracs) {
    r <- sub[sub$m_frac == mf, ]
    cat(sprintf("%-4d %8.2f %10.4f %10.4f %10.4f %10.3f %10.1f\n",
      d, mf, mean(r$md_sub), sd(r$md_sub),
      mean(r$md_sub) - md_full,
      mean(r$time),
      if (mf < 1.0) t_full / mean(r$time) else 1.0))
  }
  cat("\n")
}
cat("Done.\n")
