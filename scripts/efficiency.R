library(doParallel)
library(doRNG)

cl <- makeCluster(6)
registerDoParallel(cl)
registerDoRNG(seed = 123)

message("[", Sys.time(), "] Cluster started and RNG registered.")

n_reps <- 2000

norms_to_test <- c("none", "frobenius", "max_abs")

for (d in ds) {
  message("[", Sys.time(), "] Starting loop for d = ", d, ".")

  filename <- sprintf("~/Research/spBSS/data/%s/data_%d.rds", setting, d)

  data <- readRDS(filename)
  field <- gen_field(d)

  bd <- c(0, 1, 1, 2, 2, 3)

  rings_all <- gen_rings(
    field,
    bd,
    row_standardize = FALSE,
    symmetrize = TRUE,
    f0 = TRUE
  )

  kernels_sparse_f0 <- rings_all$kernels_sparse
  kernels_f0 <- rings_all$kernels

  kernels_sparse <- kernels_sparse_f0[-1]
  kernels <- kernels_f0[-1]

  Sf_list_standard <- lapply(kernels, function(kmat) sum(kmat) / nrow(kmat))
  Sf_list_f0 <- lapply(kernels_f0, function(kmat) sum(kmat) / nrow(kmat))

  message("[", Sys.time(), "] Starting dorng...")

  df <- foreach(r = 1:n_reps, .combine = rbind, .packages = workers) %dopar%
    {
      # THE CRITICAL FIX: Force coercion to base matrix
      sources <- as.matrix(data[[r]])

      p <- ncol(sources)

      A <- matrix(rnorm(p * p), p, p)
      X <- tcrossprod(sources, A)

      spJADE_md <- spJADE(X, A, kernels, kernels_sparse, Sf_list_standard)$md
      spJADE_md_f0 <- spJADE(X, A, kernels_f0, kernels_sparse_f0, Sf_list_f0)$md
      JADE_md <- JADE::MD(JADE::JADE(X)$W, A)

      # compile Results
      rep_res <- data.frame(
        spJADE = spJADE_md,
        spJADE_f0 = spJADE_md_f0,
        JADE = JADE_md,
        stringsAsFactors = FALSE
      )

      for (n_scheme in norms_to_test) {
        jade_md <- spJADE(
          X,
          A,
          kernels,
          kernels_sparse,
          Sf_list_standard,
          norm_scheme = n_scheme
        )$md
        rep_res[[paste0("spJADE_", n_scheme)]] <- jade_md

        jade_md_f0 <- spJADE(
          X,
          A,
          kernels_f0,
          kernels_sparse_f0,
          Sf_list_f0,
          norm_scheme = n_scheme
        )$md
        rep_res[[paste0("spJADE_f0_", n_scheme)]] <- jade_md_f0
      }

      rep_res
    }

  output_dir <- sprintf("~/Research/spBSS/results_norm_test/raw/%s", setting)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  output_file <- file.path(output_dir, sprintf("res_%d.rds", d))

  message("[", Sys.time(), "] Saving results to ", output_file)

  # Because df is already fully combined by foreach, we can save it directly
  saveRDS(df, output_file)

  message("[", Sys.time(), "] Finished d = ", d)

  rm(
    df,
    rings_all,
    kernels,
    kernels_sparse,
    kernels_f0,
    kernels_sparse_f0,
    data,
    field,
    Sf_list_standard,
    Sf_list_f0
  )
  gc()
}

message("[", Sys.time(), "] Closing cluster...")

stopCluster(cl)

message("[", Sys.time(), "] Closed successfully.")
