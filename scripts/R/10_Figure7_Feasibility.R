# ============================================================
# Figure 7 – CF_std necC (%) vs soil GlcN (Panel C) +
#             Feasibility map (Panel B) with CFF bounds
# Raw Input : data/external/pnas.sd04.soilMurAGlcN.csv
# Outputs   : outputs/figures/Fig7_necC_panels.pdf/.tiff
#             (optional) outputs/tables/SourceData_Fig7_PanelC.csv
#             (optional) outputs/tables/SourceData_Fig7_Feasibility.csv
# ============================================================

options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(geomtextpath)  # for geom_textpath()
  library(grid)          # for unit() in theme
})

# ---- theme (journal small) ----
theme_ncomm_small <- function() {
  theme_minimal(base_size = 6, base_family = "sans") %+replace%
    theme(
      plot.title   = element_text(size = 8, face = "plain", hjust = 0),
      axis.title   = element_text(size = 7),
      axis.text    = element_text(size = 6, colour = "black"),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 6),
      strip.text   = element_text(size = 6),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      axis.ticks   = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = unit(1.5, "pt"),
      legend.key.width  = unit(10, "mm"),
      legend.key.height = unit(3, "mm"),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin  = margin(3, 4, 3, 4),
      panel.spacing = unit(6, "pt"),
      panel.grid.minor = element_blank()
    )
}

# ---- (optional) assert script location if running inside the repo ----
script_rel <- file.path("scripts", "R", "10_Figure7_NecC_Panels.R")
if (file.exists(here::here(script_rel))) {
  here::i_am(script_rel)
} else {
  message("Running without here::i_am() assertion. Root = ", here::here())
}

# ---- I/O ----
infile_raw <- here::here("data", "soil.data.csv")
fig_dir    <- here::here("outputs", "figures")
tab_dir    <- here::here("outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

outfile_pdf   <- here::here(fig_dir, "Fig7_necC_panels.pdf")
outfile_tiff  <- here::here(fig_dir, "Fig7_necC_panels.tiff")

# (optional) write computed source data for reproducibility
write_sources <- TRUE
outfile_panelC_src <- here::here(tab_dir, "SourceData_Fig7_PanelC.csv")
outfile_feas_src   <- here::here(tab_dir, "SourceData_Fig7_Feasibility.csv")

# ---- read raw soil data ----
stopifnot("Raw input not found at data/external/..." = file.exists(infile_raw))
raw <- readr::read_csv(infile_raw, show_col_types = FALSE) %>%
  rename(
    MurA = tidyselect::any_of(c("MurA","soil_MurA","soilMurA","MS","murA")),
    GlcN = tidyselect::any_of(c("GlcN","soil_GlcN","soilGlcN","GS","glcN")),
    SOC  = tidyselect::any_of(c("SOC","S","soc"))
  )

# ---- outlier screen (as in your scripts) ----
dat <- raw %>%
  filter(MurA <= 412, GlcN <= 7800) %>%
  mutate(
    obs_id   = row_number(),
    soilMurA = MurA / 1000,   # convert to g g^-1
    soilGlcN = GlcN / 1000
  )

# ---- compute CF_std necC (%) ----
# Rationale: CF_std = 100 * ((soilMurA * 45) + (soilGlcN * 9)) / SOC
dat_cfstd <- dat %>%
  transmute(
    obs_id, SOC, soilGlcN,
    CF_std = 100 * ((soilMurA * 45) + (soilGlcN * 9)) / SOC
  )

# =========================
# Source tables (computed)
# =========================

# Panel C source (wide)
src_panelC <- dat_cfstd %>%
  select(obs_id, SOC, soilGlcN, CF_std)

# Panel B source (feasibility: flag > 100% SOC under CF_std)
src_feas <- dat_cfstd %>%
  transmute(
    obs_id, SOC, soilGlcN,
    over100 = CF_std > 100
  )

if (isTRUE(write_sources)) {
  readr::write_csv(src_panelC, outfile_panelC_src)
  readr::write_csv(src_feas,   outfile_feas_src)
  message("Wrote source data:\n - ", outfile_panelC_src, "\n - ", outfile_feas_src)
}

# =========================
# Plotting
# =========================

# ---- Panel C: CF_std necC % vs GlcN, by SOC bin ----
pC <- src_panelC |>
  mutate(
    S_bin = factor(if_else(SOC > 100, "> 100", "\u2264 100"),
                   levels = c("> 100", "\u2264 100"))
  ) |>
  ggplot(aes(x = soilGlcN, y = CF_std, shape = S_bin, fill = S_bin)) +
  geom_point(
    position = position_jitter(width = 0.05, height = 0),
    size = 2.4, alpha = 0.7, stroke = 0.35, color = "black"
  ) +
  scale_shape_manual(
    values = c("> 100" = 24, "\u2264 100" = 21),
    name   = expression(italic(S)~"(mg~C~g"^{-1}~"soil)")
  ) +
  scale_fill_manual(
    values = c("> 100" = "#43BBAD", "\u2264 100" = "grey60"),
    name   = expression(italic(S)~"(mg~C~g"^{-1}~"soil)")
  ) +
  labs(
    title = expression(italic(CF)[std]~" model " * italic(f)[necC] * " %"),
    y     = expression(italic(f)[necC] * " %"),
    x     = expression(italic(G)[s]~"(mg~glucosamine~" * g^{-1} * "~soil)")
  ) +
  scale_x_continuous(
    trans  = "log1p",
    breaks = c(0, 0.5, 1, 2, 4, 8),
    labels = c("0","0.5","1","2","4","8")
  ) +
  theme_classic() + theme_ncomm_small() +
  guides(
    fill = guide_legend(
      title.position = "top",
      nrow = 1, byrow = TRUE,
      override.aes = list(shape = c(24, 21), size = 2.4, color = "black", stroke = 0.35)
    ),
    shape = "none"
  ) +
  theme(
    legend.position     = c(0.82, 0.88),
    legend.direction    = "horizontal",
    plot.title.position = "panel"
  )

# ---- Panel B: Feasibility map with CFF bounds ----
feas_dat2 <- src_feas %>%
  mutate(
    over100_lab = if_else(over100, "> 100% SOC", "\u2264 100% SOC"),
    over100_lab = factor(over100_lab, levels = c("> 100% SOC", "\u2264 100% SOC"))
  )

# Feasibility band lines
cff_hi <- 20   # stricter: feasible if y <= S / cff_hi
cff_lo <- 9    # looser:   infeasible if y >  S / cff_lo

soc_range <- range(feas_dat2$SOC, na.rm = TRUE)
glcn_max  <- max(feas_dat2$soilGlcN, na.rm = TRUE)

soc_grid <- tibble(SOC = seq(soc_range[1], soc_range[2], length.out = 400)) |>
  mutate(
    y_hi = SOC / cff_hi,
    y_lo = SOC / cff_lo
  )

# Aesthetics
col_strict <- "grey25"
col_loose  <- "#a94322"
fill_some  <- "#A4C6D6"   # feasible for some CFF
fill_inf   <- "#EAD1B0"   # infeasible
col_teal   <- "#43BBAD"

pB <- ggplot() +
  geom_ribbon(data = soc_grid, aes(SOC, ymin = y_hi, ymax = y_lo),
              fill = fill_some, alpha = 0.15) +
  geom_ribbon(data = soc_grid, aes(SOC, ymin = y_lo, ymax = glcn_max),
              fill = fill_inf,  alpha = 0.12) +
  geom_textpath(data = soc_grid, aes(SOC, y_hi, label = "CFF = 20"),
                color = col_strict, linewidth = 0.5, text_smoothing = 25,
                vjust = -0.6, family = "sans") +
  geom_textpath(data = soc_grid, aes(SOC, y_lo, label = "CFF = 9"),
                color = col_loose, linewidth = 0.5, text_smoothing = 25,
                vjust = -0.6, family = "sans") +
  geom_point(
    data = feas_dat2,
    aes(SOC, soilGlcN, shape = over100_lab, fill = over100_lab),
    size = 2.6, alpha = 0.7, stroke = 0.35, color = "black"
  ) +
  scale_shape_manual(
    values = c("> 100% SOC" = 24, "\u2264 100% SOC" = 21),
    name   = expression(italic(CF)[std])
  ) +
  scale_fill_manual(
    values = c("> 100% SOC" = col_teal, "\u2264 100% SOC" = "grey60"),
    guide  = "none"
  ) +
  guides(
    shape = guide_legend(
      title.position = "top", nrow = 1, byrow = TRUE,
      override.aes = list(size = 2.8, color = "black", stroke = 0.45,
                          fill = c(col_teal, "grey60"))
    )
  ) +
  labs(
    title = expression("Feasibility bound: " * italic(G)[s] <= italic(S) / italic(CF)[F]),
    x     = expression(italic(S)~"(mg~C~" * g^{-1} * "~soil)"),
    y     = expression(italic(G)[s]~"(mg~GlcN~" * g^{-1} * "~soil)")
  ) +
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 7.5)) +
  theme_classic() + theme_ncomm_small() +
  theme(
    legend.position     = c(0.82, 0.88),
    legend.direction    = "horizontal",
    panel.grid.minor    = element_blank(),
    panel.grid.major    = element_line(linewidth = 0.2, colour = "grey90"),
    plot.title.position = "panel"
  )

# ---- Combine & save ----
final_fig <- (pB) / pC + plot_annotation(tag_levels = 'A')

# PDF (cairo preferred for crisp text; fallback if not available)
tryCatch(
  ggsave(outfile_pdf, final_fig, width = 183, height = 120, units = "mm", device = cairo_pdf),
  error = function(e) {
    message("cairo_pdf not available; saving standard PDF.")
    ggsave(outfile_pdf, final_fig, width = 183, height = 120, units = "mm")
  }
)

ggsave(outfile_tiff, final_fig, width = 183, height = 120, units = "mm",
       dpi = 600, compression = "lzw")

message("Saved: ", outfile_pdf)
message("Saved: ", outfile_tiff)


# --- Standardized repo outputs (Fig 7) ---
dir.create(here::here("analysis","figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("manuscript","figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = here::here("analysis","figures","Fig7_panel.png"),
  plot = p_fig, width = 7, height = 5, dpi = 300
)

ggsave(
  filename = here::here("manuscript","figures","Fig7.pdf"),
  plot = p_fig, width = 7, height = 5, device = grDevices::cairo_pdf
)

message("✓ Saved Figure 7: analysis/figures/Fig7_panel.png & manuscript/figures/Fig7.pdf")
