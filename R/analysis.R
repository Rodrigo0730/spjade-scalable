# -----------------------------------------------------------------------------
# Spatial BSS Helper Functions for MPI Simulation Framework
# -----------------------------------------------------------------------------
#
# Implemented Methods:
#   - spFOBI: Fobi-based Spatial Independent Component Estimation
#   - spJADE: Jade-based Spatial Independent Component Estimation
#   - SBSS: Spatial BSS based on local covariance matrices.
#   - FOBI: Fourth Order Blind Identification.
#   - JADE: Joint Approximate Diagonalization of Eigen-matrices.
#
#
# Requirements:
#   - Relies on helper data generation functions sourced from "datagen.R".
#
# -----------------------------------------------------------------------------

whiten_data <- function(X, na.action = na.fail) {
  # Whitening exactly as in from JADE:::FOBI

  X <- na.action(X)
  X <- as.matrix(X)

  X_centered <- sweep(X, 2, colMeans(X), "-")

  COV <- crossprod(X_centered) / nrow(X_centered)

  EVD <- eigen(COV, symmetric = TRUE)
  U <- EVD$vectors
  Lambda <- EVD$values

  COV_inv_sqrt <- U %*% diag(1 / sqrt(Lambda)) %*% t(U)

  X_whitened <- X_centered %*% COV_inv_sqrt

  return(list(
    Y = X_whitened,
    W0 = COV_inv_sqrt
  ))
}

compute_bkl <- function(k, l, kmat, Y) {
  # For n distinct locations, computes the fourth order cross-moment
  # matrices B_f^{kl}. The computations are vectorized and kmat is
  # a matrix such that the ij-th element corresponds to f(u_i - u_j).

  n <- nrow(Y)
  vkl <- Y[, k] * Y[, l]
  w <- kmat %*% vkl
  return(as.matrix(crossprod(Y * w, Y) / n))
}

compute_ckl <- function(k, l, kmat, lcov, Y, Sf) {
  Bkl <- compute_bkl(k, l, kmat, Y)

  u <- lcov[, k, drop = FALSE]
  v <- lcov[, l, drop = FALSE]

  # Compute the penalty division immediately to save intermediate allocations
  L_Sf <- if (k == l) {
    (2 / Sf) * tcrossprod(u)
  } else {
    (1 / Sf) * (tcrossprod(u, v) + tcrossprod(v, u))
  }

  # Apply the Lcov-based Gaussian penalty
  Ckl <- Bkl - L_Sf

  # Subtract the instantaneous variance penalty (diag(p) * Sf) strictly in-place
  # This prevents R from generating and subtracting a full p x p matrix
  if (k == l) {
    diag(Ckl) <- diag(Ckl) - Sf
  }

  return(Ckl)
}

spFOBI <- function(field, X, A, kernels_sparse, eps = 1e-06, maxiter = 100) {
  # Main function to apply spFOBI to simulated sources. The function first creates
  # a mixing matrix A and computes the mixed sources as X = sources * A^T. Then, X is
  # whitened using whiten_data() and proceeds to compute the B_f matrices for joint
  # diagonalization utilizing JADE::frjd. It returns the array of matrices which underwent
  # joint diagonalization and the Minimum Distance Index value from JADE::MD of the
  # estimated unmixing matrix and the mixing matrix A.
  #
  # Here the kernels matrices must be sparse, generated, for example, using gen_rings()
  # from the datagen.R helper file.
  #
  # Possible future implementation: compare performance with weighted diagonalization.

  if (!is.numeric(field)) {
    stop("invalid field")
  }

  K <- length(kernels_sparse)
  n <- nrow(X)
  p <- ncol(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  B_array <- array(0, dim = c(p, p, K))

  for (i in 1:K) {
    kmat <- kernels_sparse[[i]]
    B_sum <- matrix(0, p, p)
    for (k in 1:p) {
      B <- compute_bkl(k, k, kmat, Y)
      B_sum <- B_sum + B
    }

    B_array[,, i] <- B_sum
  }

  # Weighted Diagonalization could be implemented as follows
  # w <- c(1, 1, 1, 1)
  # B_array <- sweep(B_array, 3, w, `*`)

  jd_res <- tryCatch(
    suppressWarnings(frjd(B_array, maxiter = maxiter, eps = eps)),
    error = function(e) {
      message("frjd failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(jd_res)) {
    md <- NA
  } else {
    W_est <- crossprod(jd_res$V, W0)
    md <- MD(W_est, A)
  }

  return(list(
    md = md,
    B_array = B_array
  ))
}

spJADE <- function(
  X,
  A,
  kernels,
  kernels_sparse,
  Sf_list,
  norm_scheme = c("none", "frobenius", "max_abs", "l1", "infinity", "spectral"),
  eps = 1e-06,
  maxiter = 100
) {
  # Main function to apply spJADE to simulated sources.

  # Match the chosen normalization scheme (defaults to "none")
  norm_scheme <- match.arg(norm_scheme)

  K <- length(kernels_sparse)
  p <- ncol(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  n <- nrow(Y)

  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels, lcov = "lcov")

  C_array <- array(0, dim = c(p, p, p * p * K))
  idx <- 1

  for (i in 1:K) {
    kmat <- kernels_sparse[[i]]
    lcov <- lcovs[[i]]
    Sf <- Sf_list[[i]]
    for (k in 1:p) {
      for (l in 1:p) {
        C_array[,, idx] <- compute_ckl(k, l, kmat, lcov, Y, Sf)
        idx <- idx + 1
      }
    }
  }

  if (norm_scheme != "none") {
    for (m in 1:dim(C_array)[3]) {
      mat <- C_array[,, m]

      norm_val <- switch(
        norm_scheme,
        "frobenius" = norm(mat, type = "F"),
        "max_abs" = max(abs(mat)),
        "l1" = norm(mat, type = "O"), # "O" computes the One norm (max col sum)
        "infinity" = norm(mat, type = "I"), # "I" computes the Infinity norm (max row sum)
        "spectral" = norm(mat, type = "2") # "2" computes the Spectral norm
      )

      if (!is.na(norm_val) && norm_val > 1e-12) {
        C_array[,, m] <- mat / norm_val
      }
    }
  }

  jd_res <- tryCatch(
    suppressWarnings(frjd(C_array, maxiter = maxiter, eps = eps)),
    error = function(e) {
      message("frjd failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(jd_res)) {
    md <- NA
    W_est <- NA
  } else {
    W_est <- crossprod(jd_res$V, W0)
    md <- MD(W_est, A)
  }

  return(list(
    md = md,
    unmixing_matrix = W_est,
    C_array = C_array
  ))
}

lcovbss_fobi_jade <- function(X, A, kernels) {
  # This function performs the same procedure as for spFOBI and spJADE but for
  # sbss (SpatialBSS::sbss), classic FOBI and JADE (JADE::FOBI; JADE::JADE). Returns
  # the Minimum Distance Index values (JADE::MD) of the estimated unmixing matrix and
  # mixing matrix A.

  K <- length(kernels)
  n <- nrow(X)
  p <- ncol(X)

  W_est_sbss <- tryCatch(
    sbss(x = X, kernel_list = kernels, lcov = "lcov")$w,
    error = function(e) NA
  )

  W_est_fobi <- tryCatch(
    FOBI(X = X)$W,
    error = function(e) NA
  )

  W_est_jade <- tryCatch(
    JADE(X = X, n.comp = p)$W,
    error = function(e) NA
  )

  md_sbss <- tryCatch(
    MD(W_est_sbss, A),
    error = function(e) NA
  )
  md_fobi <- if (!is.matrix(W_est_fobi)) NA else MD(W_est_fobi, A)
  md_jade <- if (!is.matrix(W_est_jade)) NA else MD(W_est_jade, A)

  return(list(
    md_sbss = md_sbss,
    md_fobi = md_fobi,
    md_jade = md_jade
  ))
}

jade_md <- function(X, A) {
  # This function computes the Minimum Distance Index (JADE::MD) of the estimated
  # unmixing matrix and mixing matrix A using JADE::JADE.

  W_est_jade <- tryCatch(
    JADE(X = X, n.comp = ncol(X))$W,
    error = function(e) NA
  )

  md_jade <- if (!is.matrix(W_est_jade)) NA else MD(W_est_jade, A)

  return(md_jade)
}


nn_spjade <- function(
  X,
  A,
  kernels,
  kernels_sparse,
  Sf_list,
  lr = 0.01,
  maxiter = 500,
  lambda_orth = 100,
  seed = 42
) {
  # Exact neural spJADE: Adam solves the same JD objective as frjd.
  #
  # Precomputes R_f = F @ (Y[:,k] ⊙ Y[:,l]) for all (k,l) — O(p²·K·nnz).
  # Per epoch: batches all (a,b) dot products via Q^T @ R_stack — 1 matmul.
  #
  # Equivalent to classical spJADE — same objective, different solver.

  K <- length(kernels_sparse)
  p <- ncol(X)
  n <- nrow(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels, lcov = "lcov")

  R_stack <- list()
  for (f in 1:K) {
    kmat <- kernels_sparse[[f]]
    R_f <- matrix(0, n, p * p)
    col <- 1
    for (k in 1:p) {
      for (l in 1:p) {
        v_kl <- Y[, k] * Y[, l]
        R_f[, col] <- as.numeric(kmat %*% v_kl)
        col <- col + 1
      }
    }
    R_stack[[f]] <- R_f
  }

  off_pairs <- list()
  off_rows <- integer(0)
  off_cols <- integer(0)
  idx <- 1
  for (a in 1:p) {
    for (b in 1:p) {
      if (a != b) {
        off_pairs[[idx]] <- c(a, b)
        off_rows[idx] <- a
        off_cols[idx] <- b
        idx <- idx + 1
      }
    }
  }
  n_off <- length(off_pairs)

  Sf_vec <- unlist(Sf_list)

  Yt  <- torch_tensor(Y, dtype = torch_float32())
  Rt  <- lapply(R_stack, function(R) torch_tensor(R, dtype = torch_float32()))
  Lt  <- lapply(lcovs, function(L) torch_tensor(L, dtype = torch_float32()))
  Sft <- torch_tensor(Sf_vec, dtype = torch_float32())

  off_rows_t <- torch_tensor(off_rows, dtype = torch_long())
  off_cols_t <- torch_tensor(off_cols, dtype = torch_long())

  pa1 <- sapply(off_pairs, `[`, 1)
  pa2 <- sapply(off_pairs, `[`, 2)
  idx1_t <- torch_tensor(pa1, dtype = torch_long())
  idx2_t <- torch_tensor(pa2, dtype = torch_long())

  set.seed(seed)
  U_init <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  U <- torch_tensor(U_init, dtype = torch_float32(), requires_grad = TRUE)

  opt <- optim_adam(list(U), lr = lr)
  Id <- torch_eye(p, dtype = torch_float32())

  loss_history <- numeric(maxiter)

  for (epoch in 1:maxiter) {
    opt$zero_grad()

    S_hat <- torch_matmul(Yt, U$transpose(1, 2))

    S1 <- S_hat$index_select(2, idx1_t)
    S2 <- S_hat$index_select(2, idx2_t)
    Q <- S1 * S2

    L <- torch_tensor(0, dtype = torch_float32())

    for (f in 1:K) {
      Vals <- torch_matmul(Q$transpose(1, 2), Rt[[f]]) / n

      V_L <- torch_matmul(U, Lt[[f]])
      Sf  <- Sft[f]

      A_pen <- V_L[off_rows_t, ]
      B_pen <- V_L[off_cols_t, ]

      A_exp <- A_pen$view(c(-1, p, 1L))
      B_exp <- B_pen$view(c(-1, 1L, p))
      P <- A_exp * B_exp
      P_sym <- P + P$transpose(2, 3)
      pen_mat <- P_sym / Sf
      pen_vec <- pen_mat$view(c(n_off, p * p))

      L <- L + (Vals - pen_vec)$pow(2)$sum()
    }

    UtU <- torch_matmul(U, U$transpose(1, 2))
    penalty <- lambda_orth * (UtU - Id)$pow(2)$sum()

    total_loss <- L + penalty
    total_loss$backward()
    opt$step()

    loss_history[epoch] <- as.numeric(total_loss)
  }

  svdU <- svd(as.matrix(U))
  Uorth <- svdU$u %*% t(svdU$v)
  W_est <- Uorth %*% W0

  md <- tryCatch(
    MD(W_est, A),
    error = function(e) {
      message("MD computation failed: ", conditionMessage(e))
      NA
    }
  )

  return(list(
    md = md,
    unmixing_matrix = W_est,
    loss_history = loss_history
  ))
}

# ---------------------------------------------------------------------------
# nn_spjade variants: speed-optimized neural approximations
# ---------------------------------------------------------------------------

implicit_nn_spjade <- function(
  X,
  A,
  kernels_sparse,
  lr = 0.01,
  maxiter = 500,
  lambda_orth = 100,
  seed = 42
) {
  # Implicit neural spJADE: maximizes only the spatial auto-kurtosis
  # (diagonal elements of C_f^{ii}), relying on orthogonality + distinct
  # spatial structures to implicitly achieve separation.
  #
  # Loss: -Σ_{c,f} [ (1/n) Σ_{i,j} K_ij S_i²(c) S_j²(c) ]²
  # Cost per epoch: O(K·nnz) — only sparse kernel matvecs, no tensor expansion.

  K <- length(kernels_sparse)
  p <- ncol(X)
  n <- nrow(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  Yt <- torch_tensor(Y, dtype = torch_float32())
  Kt <- lapply(kernels_sparse, function(kmat) {
    torch_tensor(as.matrix(kmat), dtype = torch_float32())
  })

  set.seed(seed)
  U_init <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  U <- torch_tensor(U_init, dtype = torch_float32(), requires_grad = TRUE)

  opt <- optim_adam(list(U), lr = lr)
  Id <- torch_eye(p, dtype = torch_float32())

  loss_history <- numeric(maxiter)

  for (epoch in 1:maxiter) {
    opt$zero_grad()

    S_hat <- torch_matmul(Yt, U$transpose(1, 2))
    S_sq  <- S_hat$pow(2)

    L <- torch_tensor(0, dtype = torch_float32())

    for (f in 1:K) {
      B_f <- (S_sq * torch_matmul(Kt[[f]], S_sq))$sum(dim = 1) / n
      L <- L - B_f$pow(2)$sum()
    }

    UtU <- torch_matmul(U, U$transpose(1, 2))
    penalty <- lambda_orth * (UtU - Id)$pow(2)$sum()

    total_loss <- L + penalty
    total_loss$backward()
    opt$step()

    loss_history[epoch] <- as.numeric(total_loss)
  }

  svdU <- svd(as.matrix(U))
  Uorth <- svdU$u %*% t(svdU$v)
  W_est <- Uorth %*% W0

  md <- tryCatch(MD(W_est, A), error = function(e) NA)

  return(list(
    md = md,
    unmixing_matrix = W_est,
    loss_history = loss_history
  ))
}

riemann_nn_spjade <- function(
  X,
  A,
  kernels,
  kernels_sparse,
  Sf_list,
  lr = 0.01,
  maxiter = 100,
  seed = 42,
  betas = c(0.9, 0.999),
  eps = 1e-8
) {
  # Riemannian Adam on the Stiefel manifold: optimizes the exact JD loss
  # with hard orthogonality at every step (QR retraction), eliminating the
  # penalty term and final SVD. Converges in 50-100 iterations vs 500.
  #
  # Gradient projected to Stiefel tangent space before Adam update.

  K <- length(kernels_sparse)
  p <- ncol(X)
  n <- nrow(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels, lcov = "lcov")

  R_stack <- list()
  for (f in 1:K) {
    kmat <- kernels_sparse[[f]]
    R_f <- matrix(0, n, p * p)
    col <- 1
    for (k in 1:p) {
      for (l in 1:p) {
        v_kl <- Y[, k] * Y[, l]
        R_f[, col] <- as.numeric(kmat %*% v_kl)
        col <- col + 1
      }
    }
    R_stack[[f]] <- R_f
  }

  off_pairs <- list()
  off_rows <- integer(0)
  off_cols <- integer(0)
  idx <- 1
  for (a in 1:p) {
    for (b in 1:p) {
      if (a != b) {
        off_pairs[[idx]] <- c(a, b)
        off_rows[idx] <- a
        off_cols[idx] <- b
        idx <- idx + 1
      }
    }
  }
  n_off <- length(off_pairs)

  pa1 <- sapply(off_pairs, `[`, 1)
  pa2 <- sapply(off_pairs, `[`, 2)

  Sf_vec <- unlist(Sf_list)

  Yt  <- torch_tensor(Y, dtype = torch_float32())
  Rt  <- lapply(R_stack, function(R) torch_tensor(R, dtype = torch_float32()))
  Lt  <- lapply(lcovs, function(L) torch_tensor(L, dtype = torch_float32()))
  Sft <- torch_tensor(Sf_vec, dtype = torch_float32())

  off_rows_t <- torch_tensor(off_rows, dtype = torch_long())
  off_cols_t <- torch_tensor(off_cols, dtype = torch_long())
  idx1_t <- torch_tensor(pa1, dtype = torch_long())
  idx2_t <- torch_tensor(pa2, dtype = torch_long())

  set.seed(seed)
  U_init <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  U <- torch_tensor(U_init, dtype = torch_float32())

  M <- torch_zeros(c(p, p), dtype = torch_float32())
  V <- torch_zeros(c(p, p), dtype = torch_float32())

  loss_history <- numeric(maxiter)

  for (epoch in 1:maxiter) {
    U$requires_grad_(TRUE)

    S_hat <- torch_matmul(Yt, U$transpose(1, 2))
    S1 <- S_hat$index_select(2, idx1_t)
    S2 <- S_hat$index_select(2, idx2_t)
    Q <- S1 * S2

    L <- torch_tensor(0, dtype = torch_float32())

    for (f in 1:K) {
      Vals <- torch_matmul(Q$transpose(1, 2), Rt[[f]]) / n
      V_L <- torch_matmul(U, Lt[[f]])
      Sf  <- Sft[f]

      A_pen <- V_L[off_rows_t, ]
      B_pen <- V_L[off_cols_t, ]
      A_exp <- A_pen$view(c(-1, p, 1L))
      B_exp <- B_pen$view(c(-1, 1L, p))
      P <- A_exp * B_exp
      P_sym <- P + P$transpose(2, 3)
      pen_vec <- (P_sym / Sf)$view(c(n_off, p * p))

      L <- L + (Vals - pen_vec)$pow(2)$sum()
    }

    L$backward()
    G <- U$grad

    U$requires_grad_(FALSE)

    UtG <- torch_matmul(U$transpose(1, 2), G)
    G_tan <- G - 0.5 * torch_matmul(U, UtG + UtG$transpose(1, 2))

    M <- betas[1] * M + (1 - betas[1]) * G_tan
    V <- betas[2] * V + (1 - betas[2]) * G_tan$pow(2)

    m_hat <- M / (1 - betas[1]^epoch)
    v_hat <- V / (1 - betas[2]^epoch)

    U_step <- U - lr * m_hat / (v_hat$sqrt() + eps)

    U_step_r <- as.matrix(U_step)
    qr_res   <- qr(U_step_r)
    U_new    <- qr.Q(qr_res) %*% diag(sign(diag(qr.R(qr_res))), p, p)
    U        <- torch_tensor(U_new, dtype = torch_float32())

    loss_history[epoch] <- as.numeric(L)
  }

  W_est <- as.matrix(U) %*% W0
  md <- tryCatch(MD(W_est, A), error = function(e) NA)

  return(list(
    md = md,
    unmixing_matrix = W_est,
    loss_history = loss_history
  ))
}

rff_nn_spjade <- function(
  X,
  A,
  kernels,
  kernels_sparse,
  Sf_list,
  D = 100,
  lr = 0.01,
  maxiter = 500,
  lambda_orth = 100,
  seed = 42
) {
  # Random Fourier Feature approximation of spatial kernel matvecs.
  # R_f^{k,l} = K @ v_kl ≈ Φ @ (Φ^T @ v_kl) where Φ is n×D.
  # Penalty computation uses exact lcovs (still O(n·p·avg_deg), cheap).
  # Speedup: O(n·D) kernel matvecs vs O(n·avg_deg) for sparse K.
  # Gains appear when D ≪ avg_deg (e.g., large d with many neighbors).

  K <- length(kernels_sparse)
  p <- ncol(X)
  n <- nrow(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels, lcov = "lcov")

  set.seed(seed)
  Omega   <- matrix(rnorm(D, 0, 1), D, 1)
  b_coeff <- runif(D, 0, 2 * pi)
  coords  <- matrix(1:n, n, 1)
  Phi     <- sqrt(2 / D) * cos(coords %*% t(Omega) + matrix(b_coeff, n, D, byrow = TRUE))

  R_fft <- list()
  for (f in 1:K) {
    R_f <- matrix(0, n, p * p)
    col <- 1
    for (k in 1:p) {
      for (l in 1:p) {
        v_kl <- Y[, k] * Y[, l]
        # RFF approximation: K @ v_kl ≈ Φ @ (Φ^T @ v_kl)
        R_f[, col] <- as.numeric(Phi %*% (t(Phi) %*% v_kl))
        col <- col + 1
      }
    }
    R_fft[[f]] <- R_f
  }

  off_pairs <- list()
  off_rows <- integer(0)
  off_cols <- integer(0)
  idx <- 1
  for (a in 1:p) {
    for (b in 1:p) {
      if (a != b) {
        off_pairs[[idx]] <- c(a, b)
        off_rows[idx] <- a
        off_cols[idx] <- b
        idx <- idx + 1
      }
    }
  }
  n_off <- length(off_pairs)
  pa1 <- sapply(off_pairs, `[`, 1)
  pa2 <- sapply(off_pairs, `[`, 2)
  Sf_vec <- unlist(Sf_list)

  Yt  <- torch_tensor(Y, dtype = torch_float32())
  Rt  <- lapply(R_fft, function(R) torch_tensor(R, dtype = torch_float32()))
  Lt  <- lapply(lcovs, function(L) torch_tensor(L, dtype = torch_float32()))
  Sft <- torch_tensor(Sf_vec, dtype = torch_float32())

  off_rows_t <- torch_tensor(off_rows, dtype = torch_long())
  off_cols_t <- torch_tensor(off_cols, dtype = torch_long())
  idx1_t <- torch_tensor(pa1, dtype = torch_long())
  idx2_t <- torch_tensor(pa2, dtype = torch_long())

  set.seed(seed + 1)
  U_init <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  U <- torch_tensor(U_init, dtype = torch_float32(), requires_grad = TRUE)

  opt <- optim_adam(list(U), lr = lr)
  Id <- torch_eye(p, dtype = torch_float32())

  loss_history <- numeric(maxiter)

  for (epoch in 1:maxiter) {
    opt$zero_grad()

    S_hat <- torch_matmul(Yt, U$transpose(1, 2))
    S1 <- S_hat$index_select(2, idx1_t)
    S2 <- S_hat$index_select(2, idx2_t)
    Q <- S1 * S2

    L <- torch_tensor(0, dtype = torch_float32())

    for (f in 1:K) {
      Vals <- torch_matmul(Q$transpose(1, 2), Rt[[f]]) / n
      V_L  <- torch_matmul(U, Lt[[f]])
      Sf   <- Sft[f]

      A_pen <- V_L[off_rows_t, ]
      B_pen <- V_L[off_cols_t, ]
      A_exp <- A_pen$view(c(-1, p, 1L))
      B_exp <- B_pen$view(c(-1, 1L, p))
      P <- A_exp * B_exp
      P_sym <- P + P$transpose(2, 3)
      pen_vec <- (P_sym / Sf)$view(c(n_off, p * p))

      L <- L + (Vals - pen_vec)$pow(2)$sum()
    }

    UtU <- torch_matmul(U, U$transpose(1, 2))
    penalty <- lambda_orth * (UtU - Id)$pow(2)$sum()

    total_loss <- L + penalty
    total_loss$backward()
    opt$step()

    loss_history[epoch] <- as.numeric(total_loss)
  }

  svdU <- svd(as.matrix(U))
  Uorth <- svdU$u %*% t(svdU$v)
  W_est <- Uorth %*% W0

  md <- tryCatch(MD(W_est, A), error = function(e) NA)

  return(list(
    md = md,
    unmixing_matrix = W_est,
    loss_history = loss_history
  ))
}

hybrid_nn_spjade <- function(
  X,
  A,
  kernels,
  kernels_sparse,
  Sf_list,
  lr = 0.01,
  maxiter_implicit = 50,
  maxiter_jd = 200,
  lambda_orth = 100,
  seed = 42
) {
  # Hybrid: implicit (diagonal loss) warm-start → exact JD fine-tuning.
  # Phase 1: implicit_nn_spjade for maxiter_implicit epochs (fast)
  # Phase 2: exact JD loss for maxiter_jd epochs (accurate)

  K <- length(kernels_sparse)
  p <- ncol(X)
  n <- nrow(X)

  whitening <- whiten_data(X)
  Y <- whitening$Y
  W0 <- whitening$W0

  lcovs <- local_covariance_matrix(x = Y, kernel_list = kernels, lcov = "lcov")

  R_stack <- list()
  for (f in 1:K) {
    kmat <- kernels_sparse[[f]]
    R_f <- matrix(0, n, p * p)
    col <- 1
    for (k in 1:p) {
      for (l in 1:p) {
        v_kl <- Y[, k] * Y[, l]
        R_f[, col] <- as.numeric(kmat %*% v_kl)
        col <- col + 1
      }
    }
    R_stack[[f]] <- R_f
  }

  off_pairs <- list()
  off_rows <- integer(0)
  off_cols <- integer(0)
  idx <- 1
  for (a in 1:p) {
    for (b in 1:p) {
      if (a != b) {
        off_pairs[[idx]] <- c(a, b)
        off_rows[idx] <- a
        off_cols[idx] <- b
        idx <- idx + 1
      }
    }
  }
  n_off <- length(off_pairs)
  pa1 <- sapply(off_pairs, `[`, 1)
  pa2 <- sapply(off_pairs, `[`, 2)
  Sf_vec <- unlist(Sf_list)

  Yt  <- torch_tensor(Y, dtype = torch_float32())
  Rt  <- lapply(R_stack, function(R) torch_tensor(R, dtype = torch_float32()))
  Lt  <- lapply(lcovs, function(L) torch_tensor(L, dtype = torch_float32()))
  Sft <- torch_tensor(Sf_vec, dtype = torch_float32())

  Kt <- lapply(kernels_sparse, function(kmat) {
    torch_tensor(as.matrix(kmat), dtype = torch_float32())
  })

  off_rows_t <- torch_tensor(off_rows, dtype = torch_long())
  off_cols_t <- torch_tensor(off_cols, dtype = torch_long())
  idx1_t <- torch_tensor(pa1, dtype = torch_long())
  idx2_t <- torch_tensor(pa2, dtype = torch_long())

  set.seed(seed)
  U_init <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  U <- torch_tensor(U_init, dtype = torch_float32(), requires_grad = TRUE)

  opt <- optim_adam(list(U), lr = lr)
  Id <- torch_eye(p, dtype = torch_float32())

  # --- Phase 1: Implicit (diagonal-only) warm-start ---
  for (epoch in 1:maxiter_implicit) {
    opt$zero_grad()

    S_hat <- torch_matmul(Yt, U$transpose(1, 2))
    S_sq  <- S_hat$pow(2)

    L <- torch_tensor(0, dtype = torch_float32())
    for (f in 1:K) {
      B_f <- (S_sq * torch_matmul(Kt[[f]], S_sq))$sum(dim = 1) / n
      L <- L - B_f$pow(2)$sum()
    }

    UtU <- torch_matmul(U, U$transpose(1, 2))
    penalty <- lambda_orth * (UtU - Id)$pow(2)$sum()

    (L + penalty)$backward()
    opt$step()
  }

  # --- Phase 2: Exact JD fine-tuning (fresh optimizer) ---
  opt <- optim_adam(list(U), lr = lr)
  loss_history <- numeric(maxiter_jd)
  for (epoch in 1:maxiter_jd) {
    opt$zero_grad()

    S_hat <- torch_matmul(Yt, U$transpose(1, 2))
    S1 <- S_hat$index_select(2, idx1_t)
    S2 <- S_hat$index_select(2, idx2_t)
    Q <- S1 * S2

    L <- torch_tensor(0, dtype = torch_float32())

    for (f in 1:K) {
      Vals <- torch_matmul(Q$transpose(1, 2), Rt[[f]]) / n
      V_L <- torch_matmul(U, Lt[[f]])
      Sf  <- Sft[f]

      A_pen <- V_L[off_rows_t, ]
      B_pen <- V_L[off_cols_t, ]
      A_exp <- A_pen$view(c(-1, p, 1L))
      B_exp <- B_pen$view(c(-1, 1L, p))
      P <- A_exp * B_exp
      P_sym <- P + P$transpose(2, 3)
      pen_vec <- (P_sym / Sf)$view(c(n_off, p * p))

      L <- L + (Vals - pen_vec)$pow(2)$sum()
    }

    UtU <- torch_matmul(U, U$transpose(1, 2))
    penalty <- lambda_orth * (UtU - Id)$pow(2)$sum()

    total_loss <- L + penalty
    total_loss$backward()
    opt$step()

    loss_history[epoch] <- as.numeric(total_loss)
  }

  svdU <- svd(as.matrix(U))
  Uorth <- svdU$u %*% t(svdU$v)
  W_est <- Uorth %*% W0

  md <- tryCatch(MD(W_est, A), error = function(e) NA)

  return(list(
    md = md,
    unmixing_matrix = W_est,
    loss_history = loss_history
  ))
}

# ---------------------------------------------------------------------------
# Iterative kernel refinement: jointly learn kernel and unmixing
# ---------------------------------------------------------------------------

iterative_spjade <- function(
  X,
  A,
  kernels,
  kernels_sparse,
  Sf_list,
  max_iter = 5,
  alpha = 0.5,
  sigma_scale = 1,
  tol = 1e-4
) {
  n <- nrow(X)
  K <- length(kernels)
  K_dense  <- kernels
  K_sparse <- kernels_sparse
  Sf       <- Sf_list
  geo_kernels <- kernels

  history_md   <- numeric(max_iter + 1)
  history_W    <- vector("list", max_iter + 1)

  res0 <- spJADE(X, A, K_dense, K_sparse, Sf, norm_scheme = "none")
  history_md[1] <- res0$md
  history_W[[1]] <- res0$unmixing_matrix

  if (max_iter == 0) {
    return(list(md_history = history_md, W_history = history_W,
                best_md = res0$md, best_W = res0$unmixing_matrix,
                n_iter = 0, converged = TRUE))
  }

  for (iter in 1:max_iter) {
    S_est <- X %*% t(history_W[[iter]])
    S_scaled <- scale(S_est)
    dist_S <- as.matrix(dist(S_scaled))
    sigma_s <- sigma_scale * median(dist_S[upper.tri(dist_S)])
    K_src <- exp(-dist_S^2 / (2 * sigma_s^2))

    for (f in 1:K) {
      K_blend <- alpha * geo_kernels[[f]] + (1 - alpha) * K_src
      K_dense[[f]]  <- K_blend
      K_sparse[[f]] <- as(K_blend, "dgCMatrix")
      Sf[[f]] <- sum(K_blend) / n
    }

    res <- spJADE(X, A, K_dense, K_sparse, Sf, norm_scheme = "none")
    history_md[iter + 1]  <- res$md
    history_W[[iter + 1]] <- res$unmixing_matrix

    delta_md <- abs(history_md[iter + 1] - history_md[iter])
    if (!is.na(delta_md) && delta_md < tol) {
      return(list(md_history = history_md[1:(iter + 1)],
                  W_history  = history_W[1:(iter + 1)],
                  best_md = res$md, best_W = res$unmixing_matrix,
                  n_iter = iter, converged = TRUE))
    }
  }

  best_idx <- which.min(history_md)
  return(list(md_history = history_md, W_history = history_W,
              best_md = history_md[best_idx],
              best_W = history_W[[best_idx]],
              n_iter = max_iter, converged = FALSE))
}

# ---------------------------------------------------------------------------
# Stochastic subsampling spJADE
# ---------------------------------------------------------------------------

spJADE_subsample <- function(
  X, A, field, bd, m_frac, n_stoch = 5, seed = NULL
) {
  # Subsamples m = floor(m_frac * n) locations uniformly without replacement.
  # Rebuilds spatial kernels on the subsampled grid.
  # Repeats n_stoch times; returns mean MD and individual results.
  #
  # The key insight: spatial information is redundant — the 4x4 mixing
  # matrix can be recovered from far fewer than all n locations.
  
  n <- nrow(X)
  m <- max(4, floor(m_frac * n))  # at least p locations
  if (m >= n) m <- n
  
  mds <- numeric(n_stoch)
  times <- numeric(n_stoch)
  
  for (i in 1:n_stoch) {
    if (!is.null(seed)) set.seed(seed + i * 1000 + 1)
    idx <- sort(sample(n, m))
    
    X_sub  <- X[idx, , drop = FALSE]
    A_sub  <- A  # same mixing (p x p)
    field_sub <- field[idx, , drop = FALSE]
    
    rings_sub <- gen_rings(field_sub, bd,
      row_standardize = FALSE, symmetrize = TRUE, f0 = FALSE
    )
    
    kernels_sub        <- rings_sub$kernels
    kernels_sparse_sub <- rings_sub$kernels_sparse
    Sf_list_sub        <- lapply(kernels_sub, function(kmat) sum(kmat) / nrow(kmat))
    
    t0 <- proc.time()
    res <- spJADE(X_sub, A_sub, kernels_sub, kernels_sparse_sub, Sf_list_sub, 
                  norm_scheme = "none")
    times[i] <- (proc.time() - t0)["elapsed"]
    mds[i]   <- res$md
  }
  
  list(
    md_mean  = mean(mds, na.rm = TRUE),
    md_sd    = sd(mds, na.rm = TRUE),
    md_vals  = mds,
    time_mean = mean(times),
    time_vals = times,
    m        = m,
    m_frac   = m_frac,
    n_stoch  = n_stoch
  )
}
