# -------------------------------------------------------------------------
# Helper functions for data generation MPI framework
#
# This script provides a collection of helper functions for spatial grid
# generation, neighborhood construction, spatial weights matrix handling,
# and source generation (log-spARCH and log-spGARCH processes). These 
# functions are designed to support distributed simulations or computations 
# in an MPI environment.
#
# Functions included:
#   - gen_field: Creates a grid of coordinates for a d x d field.
#   - gen_rings: Builds ring-shaped spatial kernel matrices (dense and sparse).
#   - gen_nb: Generates rook/queen adjacency-based neighbors.
#   - nb2W: Converts neighbors lists to sparse weights matrices.
#   - nth_W: Constructs k-th order spatial weight matrices.
#   - gen_sources_setting_n: simulate sources
# 


gen_field <- function(d) {
  
  # Generates a d x d grid of coordinates (field) as a matrix with
  # columns x and y. Returns a matrix with d^2 rows, where each row
  # represents a grid point.

  xs <- rep(1:d, each = d)
  ys <- rep(1:d, times = d)
  coords <- cbind(x = xs, y = ys)
  return(coords)
}

gen_rings <- function(field, bandwidths, row_standardize = FALSE, symmetrize = TRUE, f0 = FALSE) {
  
  # Generates a series of ring-shaped spatial kernel matrices (both dense and 
  # sparse formats) for a given field and a set of bandwidths. 
  # Includes toggles for row-standardization, symmetrization, and appending 
  # an identity matrix to represent the kernel f0 (no spatial information).
  
  n <- nrow(field)
  kernels <- spatial_kernel_matrix(field, "ring", bandwidths)
  
  kernels <- lapply(kernels, function(W) {
    if (row_standardize) {
      rs <- rowSums(W)
      rs[rs == 0] <- 1 
      W <- sweep(W, 1, rs, "/") 
    }
    
    if (symmetrize) {
      W <- (W + t(W)) / 2
    }
    
    return(W)
  })
  
  kernels_sparse <- lapply(kernels, function(W) { as(W, "dgCMatrix") })
  
  if (f0) {
    kernels <- c(list(diag(1, n, n)), kernels)
    kernels_sparse <- c(list(Matrix::Diagonal(n, x = 1)), kernels_sparse)
  }
  
  return(list(kernels_sparse = kernels_sparse, kernels = kernels))
}

gen_nb <- function(d, type = c("rook","queen")) {
  
  # Generates a neighbors list (nb) for a d x d grid based on the specified 
  # neighborhood definition ("rook" or "queen").
  
  type <- match.arg(type)
  nb <- cell2nb(nrow = d, ncol = d, type = type, torus = TRUE)
  
  return(nb)
}

nb2W <- function(nb, style = "W") {
  
  # Converts a neighbors list (nb) into a sparse spatial weights
  # matrix of class "dgCMatrix". The style parameter defines how the weights 
  # are standardized (e.g., "W" for row-standardized weights).

  lw  <- nb2listw(nb, style = style)
  W   <- as(listw2mat(lw), "dgCMatrix")
  return(W)
}

nth_W <- function(nb, k = 2, style = "W") {
  
  # Creates a spatial weights matrix of the k-th order neighbors from a given 
  # neighbors list (nb). It computes neighbors up to lag k, extracts the 
  # k-th order neighbor list, and converts it into a spatial weights matrix 
  # using the specified style (default "W" for row-standardized weights).
  
  lagged_l <- suppressWarnings(nblag(nb, maxlag = k))
  nb_k <- lagged_l[[k]]
  Wk <- nb2W(nb_k, style = style)
  return(Wk)
}

sim_SAR <- function(coords, W, rho, alpha = 0, error = "normal", df = 8, seed = NULL) {
  
  # Simulates data from a spatial Autoregressive model as described
  # in the paper. 
  
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(coords)
  I <- diag(n)
  eps <- if (error == "normal") rnorm(n) else rt(n, df = df)
  mu <- rep(alpha, n)
  X <- solve(I - rho * as.matrix(W), mu + eps)
  as.numeric(X)
}

sim_CAR <- function(coords, W, rho, alpha = 0, seed = NULL) {
  
  # Simulates data from a conditional Autoregressive model as described in
  # the paper.
  
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(coords)
  A <- diag(n) - rho * as.matrix(W)
  Sigma <- solve(A)
  mu <- rep(alpha, n)
  X <- MASS::mvrnorm(1, mu = mu, Sigma = Sigma)
  as.numeric(X)
}

sim_SEM <- function(coords, W, rho = 0.5, alpha = 0, error = "normal", df = 7, seed = NULL) {
  
  # Simulates data from a spatial Error model as described in 
  # the paper.
  
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(coords)
  I <- diag(n)
  eps <- if (error == "normal") rnorm(n) else rt(n, df=df)
  u <- as.numeric(solve(I - rho * as.matrix(W), eps))
  X <- rep(alpha, n) + u
  as.numeric(X)
}

gen_sources_setting_1 <- function(d, W1, W2, W3, seed = NULL) {
  
  # helper function to generate the sources for the third setting
  # Rho is the spatial dependence parameter and alpha represents the 
  # baseline volatility or unconditional variance. Innovations are Gaussian.

  if (!is.null(seed)) set.seed(seed)
  
  coords <- gen_field(d)
  
  sar1 <- sim_SAR(coords, W=0.9*W1, rho = 1, alpha = 0.5)
  car <- sim_CAR(coords, W=0.3*W1 + 0.4*W2, rho=1, alpha=0.5)
  sar2 <- sim_SAR(coords, W=0.8*W1, rho=1, alpha = 0.5)
  sem <- sim_SEM(coords, W= 0.2*W1 + 0.2*W2 + 0.3*W3, rho = 1, alpha=1)
  

  return(cbind(sar1, car, sem, sar2))
}

gen_sources_setting_2 <- function(d, W1, W2, W3, seed = NULL) {
  
  # helper function to generate the sources for the third setting
  # Rho is the spatial dependence parameter and alpha represents the 
  # baseline volatility or unconditional variance. Innovations are Student-t
  # expect for one source.
  
  if (!is.null(seed)) set.seed(seed)
  
  coords <- gen_field(d)
  
  sar1 <- sim_SAR(coords, W=0.9*W1, rho = 1, alpha = 0.5, error="t", df=7)
  car <- sim_CAR(coords, W=0.3*W1 + 0.4*W2, rho=1, alpha=0.5)
  sar2 <- sim_SAR(coords, W=0.8*W1, rho=1, alpha = 0.5, error="t", df=8)
  sem <- sim_SEM(coords, W= 0.2*W1 + 0.2*W2 + 0.3*W3, rho = 1, alpha=1, error="t", df=8)
  
  
  return(cbind(sar1, car, sem, sar2))
}

gen_sources_setting_3 <- function(d, W1, W2, W3, seed = NULL) {
  
  # helper function to generate the sources for the third setting
  # Rho and lambda are the
  # spatial dependence parameters for W1 (ARCH component) and
  # W2 (GARCH component) respectively, and alpha represents the 
  # baseline volatility or unconditional variance.
  
  if (!is.null(seed)) set.seed(seed)
  
  params <- list(
    list(rho = 1, lambda = 0.0, alpha = 0.5, W1 = 0.9*W1, W2 = W1),
    list(rho = 1, lambda = 0.0, alpha = 0.5, W1 = 0.8*W1, W2 = W1),
    list(rho = 1, lambda = 1, alpha = 1, W1 = 0.4*W1 + 0.3*W2 + 0.3*W3, W2 = 0.5*W1),
    list(rho = 1, lambda = 1, alpha = 0.5, W1 = 0.2*W1 + 0.7*W2, W2 = 0.2*W1 + 0.2*W2 + 0.4*W3)
  )
  
  sources_list <- lapply(seq_along(params), function(i) {
    par <- params[[i]]
    sim.spGARCH(
      n       = d^2,
      rho     = par$rho,
      alpha   = par$alpha,
      lambda  = par$lambda,
      W1      = par$W1,
      W2      = par$W2,
      type    = "log-spGARCH",
    )

  })
  
  return(do.call(cbind, sources_list))
}

gen_sources_setting_4 <- function(d, W1, W2, W3, seed = NULL) {
  
  # helper function to generate the sources for the fourth setting
  # Rho and lambda are the
  # spatial dependence parameters for W1 (ARCH component) and
  # W2 (GARCH component) respectively, and alpha represents the 
  # baseline volatility or unconditional variance. Innovations are
  # Student-t expect for one source.
  
  if (!is.null(seed)) set.seed(seed)
  
  params <- list(
    list(rho = 1, lambda = 0.0, alpha = 0.5, W1 = 0.9*W1, W2 = W1, df=7),
    list(rho = 1, lambda = 0.0, alpha = 0.5, W1 = 0.8*W1, W2 = W1, df=7),
    list(rho = 1, lambda = 1, alpha = 1, W1 = 0.4*W1 + 0.3*W2 + 0.3*W3, W2 = 0.5*W1, df=8),
    list(rho = 1, lambda = 1, alpha = 0.5, W1 = 0.2*W1 + 0.7*W2, W2 = 0.2*W1 + 0.2*W2 + 0.4*W3, df=8)
  )
  
  sources_list <- lapply(seq_along(params), function(i) {
    par <- params[[i]]
    sim.spGARCH(
      n       = d^2,
      rho     = par$rho,
      alpha   = par$alpha,
      lambda  = par$lambda,
      W1      = par$W1,
      W2      = par$W2,
      type    = "log-spGARCH",
      error   = "student_t",
      df      = par$df,
    )
    
  })
  
  return(do.call(cbind, sources_list))
}























