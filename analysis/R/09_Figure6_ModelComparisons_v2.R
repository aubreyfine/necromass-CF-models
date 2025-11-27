
# ============================================================
# Figure 6 – Microbial necromass as % of SOC across methods
# Input : data/external/pnas.sd04.soilMurAGlcN.csv
# Output: outputs/figures/Fig6.tiff
#         outputs/tables/SourceData_Fig6_long.csv
#         outputs/tables/SourceData_Fig6_summary.csv
# ============================================================

options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
  library(ggplot2)
})


here::i_am("scripts/R/09_Figure6_ModelComparisons_v2.R")
message("Project root: ", here::here())

# ---------- I/O ----------
infile     <- here::here("data", "soil.data.csv")
out_dirF   <- here::here("outputs", "figures")
out_dirT   <- here::here("outputs", "tables")
outfile_fig  <- here::here(out_dirF, "Fig6.tiff")
outfile_long <- here::here(out_dirT, "SourceData_Fig6_long.csv")
outfile_sum  <- here::here(out_dirT, "SourceData_Fig6_summary.csv")

dir.create(out_dirF, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dirT, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(infile)) {
  stop("Input file not found at: ", infile,
       "\nPlease add the soils dataset to data/external or update `infile`.")
}

# ---------- helpers ----------
pct_lab <- function(x) paste0(x, "%")
theme_small <- function() {
  theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")
}
add_median_labels <- TRUE  # toggle to FALSE for a cleaner plot

# ---------- data ----------
soil.data <- readr::read_csv(infile, show_col_types = FALSE)

dat <- soil.data %>% filter(MurA <= 412, GlcN <= 7800)

dat3 <- dat %>%
  mutate(
    obs_id   = row_number(),
    soilMurA = MurA / 1000,
    soilGlcN = GlcN / 1000
  ) %>%
  transmute(
    obs_id, SOC,
    CF_std         = 100 * ((soilMurA * 45)   + (soilGlcN * 9))  / SOC,
    CF_Hu2024      = 100 * ((soilMurA * 31.3) + pmax(soilGlcN - 1.16 * soilMurA, 0) * 10.8) / SOC,
    Cmass_MurAGlcN = 100 * ((soilMurA * 0.42) + (soilGlcN * 0.39)) / SOC,
    Cmass_GlcN     = 100 * (soilGlcN * 0.39) / SOC
  ) %>%
  pivot_longer(c(CF_std, CF_Hu2024, Cmass_MurAGlcN, Cmass_GlcN),
               names_to = "method", values_to = "f_necC") %>%
  mutate(method = factor(method, levels = c("CF_std","CF_Hu2024","Cmass_MurAGlcN","Cmass_GlcN")))

# ---------- summary ----------
sum_F6 <- dat3 %>%
  group_by(method) %>%
  summarise(
    N            = sum(is.finite(f_necC)),
    median       = median(f_necC, na.rm = TRUE),
    q25          = quantile(f_necC, 0.25, na.rm = TRUE),
    q75          = quantile(f_necC, 0.75, na.rm = TRUE),
    over100_n    = sum(f_necC > 100, na.rm = TRUE),
    over100_pct  = 100 * mean(f_necC > 100, na.rm = TRUE),
    ymax         = suppressWarnings(max(f_necC, na.rm = TRUE)),
    .groups = "drop"
  )

# ---------- plot ----------
method_labels_expr <- c(
  expression(CF[std]),
  expression(CF[Hu2024]),
  expression(C[mass][MurA+GlcN]),
  expression(C[mass][GlcN])
)
mako_cols <- c("#312142ff", "#414388ff", "#3482a4ff", "#43BBAD")
y_floor <- 0.0005

p6 <- ggplot(dat3, aes(x = method, y = f_necC, fill = method)) +
  geom_violin(width = 0.9, alpha = 0.9, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = 16, outlier.size = 0.9,
               fill = "white", color = "grey25", fatten = 1, outlier.alpha = 0.5) +
  geom_hline(yintercept = 100, linetype = "dashed", linewidth = 0.4, color = "grey40") +
  scale_fill_manual(values = mako_cols) +
  scale_x_discrete(labels = method_labels_expr) +
  scale_y_continuous(
    trans   = "sqrt",
    breaks  = c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500),
    limits  = c(0, NA),
    labels  = pct_lab,
    expand  = expansion(mult = c(0.12, 0.03))
  ) +
  theme_small() +
  labs(x = NULL, y = expression(italic(f[necC])~"%")) +
  # % >100% labels at bottom
  geom_text(
    data = sum_F6,
    aes(x = method, y = y_floor, label = paste0(round(over100_pct, 1), "% >100%")),
    vjust = 1.2, size = 2.6, fontface = "italic",
    inherit.aes = FALSE
  ) +
  coord_cartesian(clip = "off")

if (isTRUE(add_median_labels)) {
  p6 <- p6 +
    geom_point(data = sum_F6, aes(x = method, y = median),
               inherit.aes = FALSE, size = 1.6, color = "black") +
    geom_text(data = sum_F6, aes(x = method, y = median,
                                 label = paste0(round(median, 1), "%")),
              vjust = -0.6, size = 3, inherit.aes = FALSE)
}

# ---------- save ----------
readr::write_csv(dat3, outfile_long)
readr::write_csv(sum_F6, outfile_sum)

ggsave(filename = outfile_fig, plot = p6,
       width = 183, height = 100, units = "mm",
       dpi = 600, compression = "lzw")

message("Saved figure:  ", outfile_fig)
message("Saved source:  ", outfile_long)
message("Saved summary: ", outfile_sum)


# --- Standardized repo outputs (Fig 6) ---
dir.create(here::here("analysis","figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("manuscript","figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = here::here("analysis","figures","Fig6_panel.png"),
  plot = p_fig, width = 7, height = 5, dpi = 300
)

ggsave(
  filename = here::here("manuscript","figures","Fig6.pdf"),
  plot = p_fig, width = 7, height = 5, device = grDevices::cairo_pdf
)

message("✓ Saved Figure 6: analysis/figures/Fig6_panel.png & manuscript/figures/Fig6.pdf")
