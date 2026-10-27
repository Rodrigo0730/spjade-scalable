library(ggplot2)
library(dplyr)
library(tidyr)

# Operation count formulas for spJADE and its variants
# n = d^2, avg_deg ≈ 8, K = 3, frjd_sweeps = 10, deflation_iters = 10, adam_epochs = 100

ops_classical <- function(d, p) {
  3 * p^2 * d^2 * (8 + p^2) + 30 * p^5
}

ops_deflation <- function(d, p) {
  30 * p * d^2 * (8 + p)
}

ops_neural <- function(d, p) {
  # preprocessing = classical preprocessing
  3 * p^2 * d^2 * (8 + p^2) + 300 * p^4 * d^2
}

ops_fastica <- function(d, p) {
  # True FastICA: O(p * n * p * iters) per component = O(p^2 * n * iters)
  # But this is a different problem (classical ICA, not spatial)
  p^2 * d^2 * 10  # ballpark for classical FastICA
}

# Grid
ds <- c(5, 10, 20, 40, 80, 100, 200)
ps <- c(4, 10, 20, 50)

df <- expand.grid(d = ds, p = ps)
df$classical  <- with(df, ops_classical(d, p))
df$deflation  <- with(df, ops_deflation(d, p))
df$neural     <- with(df, ops_neural(d, p))
df$n           <- df$d^2

# Winner (fewest ops)
df$winner <- apply(df[, c("classical", "deflation", "neural")], 1, function(x) {
  c("classical", "deflation", "neural")[which.min(x)]
})

# Ratio
df$defl_vs_cl <- df$deflation / df$classical
df$nn_vs_cl   <- df$neural / df$classical

cat("========== Operation counts (millions) ==========\n\n")
cat(sprintf("%6s %3s %8s %12s %12s %12s %8s %10s %10s\n",
  "d", "p", "n", "classical", "deflation", "neural", "winner", "defl/cl", "nn/cl"))
cat(strrep("-", 90), "\n")

for (i in 1:nrow(df)) {
  row <- df[i, ]
  cat(sprintf("%6d %3d %8d %12.2f %12.2f %12.2f %8s %10.2f %10.2f\n",
    row$d, row$p, row$n,
    row$classical/1e6, row$deflation/1e6, row$neural/1e6,
    row$winner,
    row$defl_vs_cl, row$nn_vs_cl))
}

cat("\n========== When does deflationary win? ==========\n")
cat("Deflationary ~ 30·p·d²·(8+p),  Classical ~ 3·p²·d²·(8+p²)\n")
cat("Deflationary/Classical ≈ 10·(8+p) / [p·(8+p²)]\n")
cat("Ratio < 1 when: 10(8+p) < p(8+p²) → 80+10p < 8p+p³ → p³ > 80+2p\n")
cat("→ p > ~4.4  (crossover at p ≈ 5)\n\n")

cat("========== When does neural win? ==========\n")
cat("Neural preprocessing = classical preprocessing\n")
cat("Neural ONLY wins when preprocessing << neural solver\n")
cat("→ Happens at very large p (p > 50) where p⁴ dominates\n\n")

cat("========== Key takeaway ==========\n")
cat("For p=4: classical is ALWAYS optimal at all d.\n")
cat("For p≥5: deflationary surpasses classical (asymptotically).\n")
cat("Neural: never wins for sparse ring kernels at any (d,p).\n")
