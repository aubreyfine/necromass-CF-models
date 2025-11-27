# ============================================================
# Figure 5 – Diagnostic ratio R and structural feasibility
# Input : outputs/tables/diagnostics_input.csv
# Output: outputs/figures/Fig5_Diagnostics.tiff
# ============================================================

options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(viridis)
})

here::i_am("scripts/R/08_Figure5_Diagnostics.R")
message("Project root: ", here::here())

infile  <- here::here("outputs", "tables", "diagnostics_input.csv")
outdir  <- here::here("outputs", "figures")
outfile <- here::here(outdir, "Fig5_Diagnostics.tiff")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---------- helpers ----------
theme_small <- function() {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.tag.position = c(0, 1),
      plot.tag = element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
}

# ---------- data ----------
raw <- readr::read_csv(infile, show_col_types = FALSE)

# Standardize names; keep optional fnecC columns if present
dat <- raw %>%
  dplyr::rename(
    MS      = tidyselect::any_of(c("MS","soil_MurA")),
    GS      = tidyselect::any_of(c("GS","soil_GlcN")),
    S       = tidyselect::any_of(c("S","SOC")),
    rB      = tidyselect::any_of(c("rB","r_B")),
    CFF     = tidyselect::any_of(c("CFF","CF_fungal")),
    fnecC_int = tidyselect::any_of(c("fnecC_int", "fnecC_integral", "fnecC")),
    fnecC_agg = tidyselect::any_of(c("fnecC_agg", "fnecC_aggregate"))
  ) %>%
  dplyr::select(dplyr::any_of(c("MS","GS","S","rB","CFF","fnecC_int","fnecC_agg")))

# Compute R if missing; keep a clipped version for plotting
if (!"R" %in% names(dat)) {
  dat <- dat %>%
    dplyr::mutate(R = 100 * (rB * MS) / (CFF * GS))
}
dat <- dat %>%
  dplyr::mutate(
    R_plot = pmin(R, 100),                       # cap for plotting
    rB_bin = cut(
      rB,
      breaks = c(-Inf, 0, 1, 2, 3, 4, Inf),
      labels = c("0","0–1","1–2","2–3","3–4",">4"),
      right = TRUE
    )
  ) %>%
  dplyr::filter(!is.na(rB_bin))

# Labels for Panel A
labA <- dat %>%
  dplyr::group_by(rB_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    trunc_prop = mean(GS - rB*MS <= 0, na.rm = TRUE),
    .groups = "drop"
  )

# ---------- Panel A: R by rB bins ----------
pA <- ggplot(dat, aes(x = rB_bin, y = R_plot)) +
  geom_violin(fill = "#43BBADFF", alpha = 0.6, color = NA, width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.14, outlier.size = 0.6, outlier.alpha = 0.4, alpha = 0.55) +
  stat_summary(fun = median, geom = "point", size = 1.6) +
  geom_text(
    data = labA,
    aes(label = paste0("N=", scales::comma(n),
                       "\ntrunc=", scales::percent(trunc_prop, accuracy = 0.1)),
        y = 30),
    vjust = 1, size = 2.2, lineheight = 0.95
  ) +
  coord_cartesian(ylim = c(0, 120)) +
  annotate("text", x = 6, y = 112, label = ">5% truncated", hjust = 1, size = 2.5) +
  labs(
    x = expression(italic(r[B]) ~ "bin"),
    y = expression(italic(R)~"(% of fungal GlcN subtracted)")
  ) +
  theme_small() +
  labs(tag = "A")

# ---------- Panel B: Feasibility heatmap (preferred) or S–GS fallback ----------
has_int <- "fnecC_int" %in% names(dat) && any(is.finite(dat$fnecC_int))
has_agg <- "fnecC_agg" %in% names(dat) && any(is.finite(dat$fnecC_agg))

make_heat <- function(df, yvar, ylab) {
  # % infeasible
  pct_infeasible <- mean(df[[yvar]] > 100, na.rm = TRUE)
  
  ggplot(df, aes(x = S, y = .data[[yvar]])) +
    stat_summary_2d(aes(z = GS), bins = 60, fun = median) +
    scale_fill_viridis_c(
      option = "mako", direction = 1, begin = 0.2, end = 1,
      name = expression(italic(G)[S]~"(mg~GlcN~g~soil"^{-1}*")")
    ) +
    geom_hline(yintercept = 100, linetype = "dashed", alpha = 0.7) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 100, ymax = Inf),
              fill = "#dfab9e", inherit.aes = FALSE, alpha = 0.6) +
    annotate(
      "text", x = -Inf, y = 118, hjust = -0.05, vjust = 1, size = 2.2,
      label = paste0("Values >120% clipped; ",
                     scales::percent(pct_infeasible, accuracy = 0.1),
                     " >100%")
    ) +
    annotate(
      "text", x = 20, y = 105, hjust = 1.05, vjust = -0.2, size = 2.2,
      label = "paste('Infeasible ( > ', italic(S) / italic(CF)[italic(F)], ' )')",
      parse = TRUE
    ) +
    scale_x_continuous(
      trans = "log1p",
      breaks = c(0, 1, 10, 100, 400),
      labels = c("0","1","10","100","400")
    ) +
    scale_y_continuous(
      limits = c(0, 120), breaks = seq(0, 120, 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = expression(italic(S)~"(mg~C~g~soil"^{-1}*")"),
      y = ylab
    ) +
    theme_small()
}

if (has_int) {
  pB <- make_heat(dat, "fnecC_int", expression(italic(f[necC])~"(integral), %"))
} else if (has_agg) {
  pB <- make_heat(dat, "fnecC_agg", expression(italic(f[necC])~"(aggregate), %"))
} else {
  # Fallback: S–GS density with an illustrative feasibility guide
  pB <- ggplot(dat, aes(x = S, y = GS)) +
    geom_bin2d(bins = 40) +
    scale_fill_viridis_c(name = "Count") +
    # Illustrative guide: CFF ≈ 9 → GS ≈ S/9 (units simplified)
    geom_abline(slope = 1/9, intercept = 0, color = "red", linetype = 2) +
    labs(
      x = expression(italic(S)~"(mg~C~g~soil"^{-1}*")"),
      y = expression(italic(G)[S]~"(mg~GlcN~g~soil"^{-1}*")")
    ) +
    theme_small()
}

pB <- pB + labs(tag = "B")

# ---------- Combine & save ----------
fig5 <- pA + pB + plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = outfile,
  plot = fig5,
  dpi = 600,
  width = 90, height = 160, units = "mm",
  compression = "lzw"
)

message("Saved: ", outfile)

# --- Standardized repo outputs (Fig 5) ---
dir.create(here::here("analysis","figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("manuscript","figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = here::here("analysis","figures","Fig5_panel.png"),
  plot = p_fig, width = 7, height = 5, dpi = 300
)

ggsave(
  filename = here::here("manuscript","figures","Fig5.pdf"),
  plot = p_fig, width = 7, height = 5, device = grDevices::cairo_pdf
)

message("✓ Saved Figure 5: analysis/figures/Fig5_panel.png & manuscript/figures/Fig5.pdf")

