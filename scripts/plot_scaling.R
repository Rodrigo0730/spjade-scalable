# Figures for computation time against field size and the subsampling
# trade-off. Reads the outputs of benchmark_all.R and benchmark_subsample.R
# and writes figures/time_vs_n.pdf and figures/subsample_tradeoff.pdf.

library(ggplot2)
library(dplyr)
library(tidyr)

SETTING <- "setting_3"
INPUT_DIR <- file.path("results", "benchmark_all", "raw", SETTING)
OUTPUT_DIR <- "figures"
dir.create(OUTPUT_DIR, showWarnings = FALSE)

method_labels <- c(
  spjade = "spJADE",
  nn = "nn",
  riemann = "riemann",
  hybrid = "hybrid",
  implicit = "implicit",
  rff = "rff"
)
method_colors <- c(
  spJADE = "#2a78d6",
  nn = "#eb6834",
  riemann = "#1baf7a",
  hybrid = "#eda100",
  implicit = "#e87ba4",
  rff = "#008300"
)
method_shapes <- c(
  spJADE = 16, nn = 17, riemann = 15, hybrid = 18, implicit = 25, rff = 4
)

theme_bench <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.3),
    axis.text = element_text(colour = "#52514e"),
    axis.title = element_text(colour = "#0b0b0b"),
    plot.title = element_text(face = "bold", colour = "#0b0b0b"),
    plot.subtitle = element_text(colour = "#52514e"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# ---- Time against n ---------------------------------------------------------

timing <- list.files(INPUT_DIR, pattern = "^res_\\d+\\.rds$", full.names = TRUE) |>
  lapply(function(f) as.data.frame(readRDS(f))) |>
  bind_rows() |>
  select(d, starts_with("time_")) |>
  pivot_longer(-d, names_to = "method", names_prefix = "time_", values_to = "time") |>
  mutate(n = d^2, method = factor(method_labels[method], levels = method_labels)) |>
  group_by(method, n) |>
  summarise(
    median = median(time, na.rm = TRUE),
    q25 = quantile(time, 0.25, na.rm = TRUE),
    q75 = quantile(time, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Methods whose medians nearly coincide at the largest n share one end label
end_labels <- timing |>
  filter(n == max(n)) |>
  arrange(median) |>
  mutate(group = cumsum(c(TRUE, diff(log10(median)) > 0.05))) |>
  group_by(group) |>
  summarise(
    n = first(n),
    median = exp(mean(log(median))),
    label = paste(method, collapse = ", "),
    .groups = "drop"
  )

p_time <- ggplot(timing, aes(n, median, colour = method, shape = method)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.04, linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.2, fill = "white") +
  geom_text(
    data = end_labels,
    aes(n, median, label = label),
    inherit.aes = FALSE,
    hjust = 0, nudge_x = 0.06, size = 3, colour = "#52514e"
  ) +
  scale_x_log10(
    breaks = sort(unique(timing$n)),
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  scale_y_log10(labels = scales::label_number(accuracy = 0.01)) +
  scale_colour_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(
    title = "Computation time against field size",
    subtitle = "Setting 3, median over replications with interquartile range",
    x = "Number of locations n (log scale)",
    y = "Seconds per fit (log scale)"
  ) +
  theme_bench

ggsave(file.path(OUTPUT_DIR, "time_vs_n.pdf"), p_time, width = 7, height = 4.5)

# ---- Subsampling trade-off --------------------------------------------------

sub <- readRDS(file.path("results", "subsample_benchmark.rds")) |>
  group_by(d, m_frac) |>
  summarise(
    md = median(md_sub, na.rm = TRUE),
    time = median(time, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(n = factor(paste0("n = ", d^2), levels = paste0("n = ", sort(unique(d))^2))) |>
  arrange(d, m_frac)

n_colors <- setNames(c("#2a78d6", "#eb6834", "#1baf7a"), levels(sub$n))

p_sub <- ggplot(sub, aes(time, md, colour = n)) +
  geom_path(linewidth = 0.6) +
  geom_point(size = 2.2) +
  geom_text(
    aes(label = scales::percent(m_frac, accuracy = 1)),
    hjust = -0.3, size = 2.6, colour = "#52514e", show.legend = FALSE
  ) +
  facet_wrap(~n, scales = "free_x") +
  scale_x_continuous(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.08, 0.3))
  ) +
  scale_colour_manual(values = n_colors, guide = "none") +
  labs(
    title = "Subsampling: time against accuracy",
    subtitle = "Setting 3, median over replications; labels give the fraction of locations used",
    x = "Seconds per fit",
    y = "Minimum distance index (lower is better)"
  ) +
  theme_bench

ggsave(file.path(OUTPUT_DIR, "subsample_tradeoff.pdf"), p_sub, width = 9, height = 3.8)
