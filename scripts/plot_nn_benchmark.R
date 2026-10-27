# -----------------------------------------------------------------------------
# Analysis: nn_spjade vs spJADE benchmark results
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

INPUT_DIR  <- "results/nn_benchmark/raw/setting_3"
OUTPUT_DIR <- "results/nn_benchmark/plots"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

ds <- c(5, 10, 20)

all_dfs <- lapply(ds, function(d) {
  f <- file.path(INPUT_DIR, sprintf("res_%d.rds", d))
  if (!file.exists(f)) return(NULL)
  df <- readRDS(f)
  df$d <- as.factor(paste0("d = ", df$d))
  df
})

df_all <- do.call(rbind, all_dfs)

# ---------------------------------------------------------------------------
# 1. Summary statistics table
# ---------------------------------------------------------------------------
cat("\n========== Summary Statistics ==========\n\n")

summary_tab <- df_all %>%
  group_by(d) %>%
  summarise(
    n = n(),
    spJADE_mean_MD = mean(md_spjade, na.rm = TRUE),
    spJADE_med_MD  = median(md_spjade, na.rm = TRUE),
    spJADE_sd_MD   = sd(md_spjade, na.rm = TRUE),
    nn_mean_MD     = mean(md_nn, na.rm = TRUE),
    nn_med_MD      = median(md_nn, na.rm = TRUE),
    nn_sd_MD       = sd(md_nn, na.rm = TRUE),
    ΔMD_mean       = mean(md_nn - md_spjade, na.rm = TRUE),
    ΔMD_sd         = sd(md_nn - md_spjade, na.rm = TRUE),
    spJADE_t_med   = median(time_spjade, na.rm = TRUE),
    nn_t_med       = median(time_nn, na.rm = TRUE),
    jade_mean_MD   = mean(md_jade, na.rm = TRUE),
    jade_t_med     = median(time_jade, na.rm = TRUE),
    .groups = "drop"
  )

print(as.data.frame(summary_tab), digits = 4)

# ---------------------------------------------------------------------------
# 2. Paired t-test / Wilcoxon test
# ---------------------------------------------------------------------------
cat("\n========== MD Equality Tests ==========\n\n")

for (d_val in ds) {
  sub <- df_all[df_all$d == paste0("d = ", d_val), ]
  wt <- wilcox.test(sub$md_spjade, sub$md_nn, paired = TRUE)
  tt <- t.test(sub$md_spjade, sub$md_nn, paired = TRUE)
  cat(sprintf(
    "d = %d:  Wilcoxon p = %.2e  |  t-test p = %.2e  |  mean ΔMD = %.4f\n",
    d_val, wt$p.value, tt$p.value,
    mean(sub$md_nn - sub$md_spjade, na.rm = TRUE)
  ))
}

# ---------------------------------------------------------------------------
# 3. Correlation (nn MD vs spJADE MD)
# ---------------------------------------------------------------------------
cat("\n========== MD Correlation (nn vs spJADE) ==========\n\n")
for (d_val in ds) {
  sub <- df_all[df_all$d == paste0("d = ", d_val), ]
  cc <- cor(sub$md_spjade, sub$md_nn, use = "complete.obs")
  cat(sprintf("d = %d:  Pearson r = %.4f\n", d_val, cc))
}

# ---------------------------------------------------------------------------
# 4. Boxplot: MD by method and d
# ---------------------------------------------------------------------------
df_long <- df_all %>%
  select(d, md_spjade, md_nn, md_jade) %>%
  pivot_longer(-d, names_to = "Method", values_to = "MD") %>%
  mutate(
    Method = recode(Method,
      md_spjade = "spJADE",
      md_nn     = "nn_spjade",
      md_jade   = "JADE"
    ),
    Method = factor(Method, levels = c("spJADE", "nn_spjade", "JADE"))
  )

p_box <- ggplot(df_long, aes(x = Method, y = MD, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~ d, nrow = 1) +
  scale_fill_manual(values = c("spJADE" = "#2c7bb6", "nn_spjade" = "#d7191c", "JADE" = "#fdae61")) +
  labs(title = "MD Distribution: spJADE vs nn_spjade vs JADE (Setting 3)", y = "Minimum Distance Index") +
  theme_bw(base_size = 13) +
  theme(axis.title.x = element_blank())

ggsave(file.path(OUTPUT_DIR, "md_boxplot.pdf"), p_box, width = 10, height = 5)
cat(sprintf("\nSaved %s\n", file.path(OUTPUT_DIR, "md_boxplot.pdf")))

# ---------------------------------------------------------------------------
# 5. Scatter: nn MD vs spJADE MD
# ---------------------------------------------------------------------------
p_scatter <- ggplot(df_all, aes(x = md_spjade, y = md_nn)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ d, nrow = 1) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "nn_spjade MD vs spJADE MD",
    x = "spJADE MD", y = "nn_spjade MD"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(OUTPUT_DIR, "md_scatter.pdf"), p_scatter, width = 10, height = 4)
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "md_scatter.pdf")))

# ---------------------------------------------------------------------------
# 6. Difference distribution
# ---------------------------------------------------------------------------
df_all$delta_md <- df_all$md_nn - df_all$md_spjade

p_delta <- ggplot(df_all, aes(x = delta_md, fill = d)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "MD Difference (nn_spjade - spJADE)",
    x = expression(Delta * MD), y = "Density"
  ) +
  theme_bw(base_size = 13)

ggsave(file.path(OUTPUT_DIR, "md_delta_density.pdf"), p_delta, width = 7, height = 5)
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "md_delta_density.pdf")))

# ---------------------------------------------------------------------------
# 7. Timing comparison
# ---------------------------------------------------------------------------
df_time <- df_all %>%
  select(d, time_spjade, time_nn, time_jade) %>%
  pivot_longer(-d, names_to = "Method", values_to = "Time") %>%
  mutate(
    Method = recode(Method,
      time_spjade = "spJADE",
      time_nn     = "nn_spjade",
      time_jade   = "JADE"
    ),
    Method = factor(Method, levels = c("spJADE", "nn_spjade", "JADE"))
  )

p_time <- ggplot(df_time, aes(x = d, y = Time, fill = Method)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = c("spJADE" = "#2c7bb6", "nn_spjade" = "#d7191c", "JADE" = "#fdae61")) +
  labs(
    title = "Wall-clock Time per Repetition (Setting 3)",
    y = "Time (seconds)"
  ) +
  theme_bw(base_size = 13) +
  scale_y_log10()

ggsave(file.path(OUTPUT_DIR, "timing_boxplot.pdf"), p_time, width = 8, height = 5)
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "timing_boxplot.pdf")))

# ---------------------------------------------------------------------------
# 8. Per-repetition MD comparison (connected dots for paired comparison)
# ---------------------------------------------------------------------------
# Sample first 30 reps per d for clarity
df_sample <- df_all %>%
  group_by(d) %>%
  slice(1:30) %>%
  ungroup() %>%
  mutate(pair_id = paste(d, rep, sep = "_"))

df_paired <- bind_rows(
  df_sample %>% select(d, rep, md = md_spjade) %>% mutate(Method = "spJADE"),
  df_sample %>% select(d, rep, md = md_nn)     %>% mutate(Method = "nn_spjade")
) %>%
  mutate(Method = factor(Method, levels = c("spJADE", "nn_spjade")))

p_paired <- ggplot(df_paired, aes(x = Method, y = md, group = interaction(d, rep))) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  geom_point(aes(color = Method), size = 1, alpha = 0.5) +
  facet_wrap(~ d, nrow = 1) +
  scale_color_manual(values = c("spJADE" = "#2c7bb6", "nn_spjade" = "#d7191c")) +
  labs(title = "Paired MD: nn_spjade vs spJADE (30 reps per d)", y = "MD") +
  theme_bw(base_size = 13) +
  theme(axis.title.x = element_blank())

ggsave(file.path(OUTPUT_DIR, "md_paired.pdf"), p_paired, width = 10, height = 5)
cat(sprintf("Saved %s\n", file.path(OUTPUT_DIR, "md_paired.pdf")))

cat("\nAnalysis complete.\n")
