
# ============================================================
# Figure 4 — Trait–function relationships (A–D)
# Inputs :
#   data/simulations/CFB_sim.csv  (CFB, MGP, fGP, cB?, MGN?, …)
#   data/simulations/CFF_sim.csv  (CFF, GF,  cF,  …)
# Outputs:
#   outputs/figures/Fig4_Traits_vs_CF.tiff         (A–B)
#   outputs/figures/Fig4_Traits_vs_CF_EXT.tiff     (A–D, if cB & MGN present)
# ============================================================

options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(here)
  library(mgcv)
  library(grid)
})

# --- Pin repo root (update if you move this file) ---
here::i_am("scripts/R/07_Figure4_Traits_vs_CF.R")
message("Project root: ", here::here())

# --- NatComms small-panel theme ---
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

# ---------- I/O ----------
in_cfb <- here::here("data", "simulations", "CFB_sim.csv")
in_cff <- here::here("data", "simulations", "CFF_sim.csv")
outdir <- here::here("outputs", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile_AB   <- here::here(outdir, "Fig4_Traits_vs_CF.tiff")
outfile_ABCD <- here::here(outdir, "Fig4_Traits_vs_CF_EXT.tiff")

# ---------- Read & normalize (keep all columns) ----------
CFB <- readr::read_csv(in_cfb, show_col_types = FALSE) %>%
  dplyr::rename(
    CFB = tidyselect::any_of(c("CFB","CF_bac","cfb")),
    MGP = tidyselect::any_of(c("MGP","M_GP")),
    fGP = tidyselect::any_of(c("fGP","frac_GP","F_GP")),
    cB  = tidyselect::any_of(c("cB","c_B","c_bac")),
    MGN = tidyselect::any_of(c("MGN","M_GN"))
  ) %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::any_of(c("CFB","MGP","fGP","cB","MGN")),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )

CFF <- readr::read_csv(in_cff, show_col_types = FALSE) %>%
  dplyr::rename(
    CFF = tidyselect::any_of(c("CFF","CF_fun","cff")),
    GF  = tidyselect::any_of(c("GF","G_F")),
    cF  = tidyselect::any_of(c("cF","c_F"))
  ) %>%
  dplyr::mutate(
    dplyr::across(
      tidyselect::any_of(c("CFF","GF","cF")),
      ~ suppressWarnings(as.numeric(as.character(.)))
    )
  )

# Guards
stopifnot(all(c("CFB","MGP","fGP") %in% names(CFB)))
stopifnot(all(c("CFF","GF","cF")  %in% names(CFF)))
has_cB  <- "cB"  %in% names(CFB) && any(is.finite(CFB$cB))
has_MGN <- "MGN" %in% names(CFB) && any(is.finite(CFB$MGN))

# ---------- Robust ranges ----------
ylims_cfb <- quantile(CFB$CFB, c(0.01, 0.99), na.rm = TRUE)
ylims_cff <- quantile(CFF$CFF, c(0.01, 0.99), na.rm = TRUE)
fgp_lim   <- quantile(CFB$fGP, c(0.02, 0.98), na.rm = TRUE)
cF_lim    <- quantile(CFF$cF,  c(0.02, 0.98), na.rm = TRUE)
cB_lim    <- if (has_cB)  quantile(CFB$cB,  c(0.02, 0.98), na.rm = TRUE) else NULL
MGN_lim   <- if (has_MGN) quantile(CFB$MGN, c(0.02, 0.98), na.rm = TRUE) else NULL

# ---------- Helper: tile + GAM (smooth vs log10(x)) ----------
panel_tiles_gam <- function(dat, x, y, z, x_breaks, y_breaks, ylims, fill_name, fill_limits){
  ggplot(dat, aes(x = {{x}}, y = {{y}})) +
    stat_summary_2d(aes(z = {{z}}), bins = 50, fun = median, drop = TRUE, na.rm = TRUE) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(log10(x), k = 6),   # smooth along log10(x)
      color = "black", linewidth = 0.6, se = TRUE, alpha = 0.12
    ) +
    scale_fill_viridis_c(
      option = "mako", direction = 1, begin = 0.2, end = 1,
      limits = fill_limits, oob = squish, name = fill_name
    ) +
    scale_x_continuous(trans = "log10", breaks = x_breaks) +
    scale_y_continuous(trans = "log10", breaks = y_breaks,
                       labels = number_format(accuracy = 1)) +
    coord_cartesian(ylim = ylims) +
    theme_classic() + theme_ncomm_small()
}

library(ggstar)
# ---------- Panels ----------
# A: CFB ~ MGP | fill = fGP
pA <- panel_tiles_gam(
  CFB, MGP, CFB, fGP,
  x_breaks = c(2,5,10,20,40),
  y_breaks = c(3,10,30,100),
  ylims    = ylims_cfb,
  fill_name = expression(italic(f[GP])),
  fill_limits = fgp_lim
) + labs(x = expression(italic(M[GP])~"(mg MurA"~g^{-1}~"bacterial biomass)"),
         y = expression(italic(CF[B]))) +
  geom_star(aes(x = 13.9, y = 45), 
             shape = 21, fill = "white", color = "black", size = 2.8, stroke = 0.6) 
pA

# B: CFF ~ GF | fill = cF
pB <- panel_tiles_gam(
  CFF, GF, CFF, cF,
  x_breaks = c(5,10,20,50,100,200),
  y_breaks = c(3,10,30,100),
  ylims    = ylims_cff,
  fill_name = expression(italic(c[F])),
  fill_limits = cF_lim
) + labs(x = expression(italic(G[F])~"(mg GlcN"~g^{-1}~"fungal biomass)"),
         y = expression(italic(CF[F]))) +
  geom_star(aes(x = 49, y = 9), 
            shape = 21, fill = "white", color = "black", size = 2.8, stroke = 0.6) 

# C: CFB ~ MGP | fill = cB   (if available)
pC <- if (has_cB) {
  panel_tiles_gam(
    CFB, MGP, CFB, cB,
    x_breaks = c(2,5,10,20,40),
    y_breaks = c(3,10,30,100),
    ylims    = ylims_cfb,
    fill_name = expression(italic(c)[B]),
    fill_limits = cB_lim
  ) + labs(x = expression(italic(M[GP])~"(mg MurA"~g^{-1}~"bacterial biomass)"),
           y = expression(italic(CF[B]))) +
    geom_star(aes(x = 13.9, y = 45), 
              shape = 21, fill = "white", color = "black", size = 2.8, stroke = 0.6) 
} else NULL

# D: CFB ~ MGP | fill = MGN (if available)
pD <- if (has_MGN) {
  panel_tiles_gam(
    CFB, MGP, CFB, MGN,
    x_breaks = c(2,5,10,20,40),
    y_breaks = c(3,10,30,100),
    ylims    = ylims_cfb,
    fill_name = expression(italic(M)[GN]),
    fill_limits = MGN_lim
  ) + labs(x = expression(italic(M[GP])~"(mg MurA"~g^{-1}~"bacterial biomass)"),
           y = expression(italic(CF[B]))) +
    geom_star(aes(x = 13.9, y = 45), 
              shape = 21, fill = "white", color = "black", size = 2.8, stroke = 0.6) 
} else NULL

# ---------- Save A–B (main) ----------
fig_AB <- (pA + pB) + plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

#pA + pB + pC + pD + plot_layout(ncol = 2)
#ggsave(filename = outfile_AB, plot = fig_AB, dpi = 600,
 #      width = 185, height = 155, units = "mm", device = "tiff", compression = "lzw")
#message("Saved main: ", normalizePath(outfile_AB, winslash = "\\"))

# ---------- Save A–D (extended) if C & D exist ----------
if (!is.null(pC) && !is.null(pD)) {
  fig_ABCD <- ((pA + pC) / (pB + pD)) +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
  ggsave(filename = outfile_ABCD, plot = fig_ABCD, dpi = 600,
         width = 185, height = 160, units = "mm", device = "tiff", compression = "lzw")
  message("Saved extended: ", normalizePath(outfile_ABCD, winslash = "\\"))
} else {
  if (is.null(pC)) message("Note: cB not found or non-finite — panel C skipped.")
  if (is.null(pD)) message("Note: MGN not found or non-finite — panel D skipped.")
}








# ---------------------------
# 1) Harmonize scales
# ---------------------------
# Shared y for CFB panels (A/C/D)
ylims_cfb <- quantile(CFB$CFB, c(0.01, 0.99), na.rm = TRUE)
ybreaks_cfb <- c(3, 10, 30, 100)

# Fungal panel
ylims_cff <- quantile(CFF$CFF, c(0.01, 0.99), na.rm = TRUE)
ybreaks_cff <- c(3, 10, 30, 100)

fix_axes_cfb <- list(
  scale_y_continuous(trans = "log10", breaks = ybreaks_cfb, labels = scales::number_format(accuracy = 1)),
  coord_cartesian(ylim = ylims_cfb)
)

fix_axes_cff <- list(
  scale_y_continuous(trans = "log10", breaks = ybreaks_cff, labels = scales::number_format(accuracy = 1)),
  coord_cartesian(ylim = ylims_cff)
)

xfix_MGP <- scale_x_continuous(trans = "log10", breaks = c(2, 5, 10, 20, 40))
xfix_GF  <- scale_x_continuous(trans = "log10", breaks = c(5, 10, 20, 50, 100, 200))

pA <- pA + fix_axes_cfb + xfix_MGP
pB <- pB + fix_axes_cff + xfix_GF
if (exists("pC") && !is.null(pC)) pC <- pC + fix_axes_cfb + xfix_MGP
if (exists("pD") && !is.null(pD)) pD <- pD + fix_axes_cfb + xfix_MGP


# ---------------------------
# 3) Legend collection & spacing
# ---------------------------
# Make legends compact and consistent
legend_tweaks <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.key.width  = unit(5, "mm"),
  legend.key.height = unit(3, "mm"),
  legend.box.margin = margin(0, 0, 0, 0)
)

pA <- pA + legend_tweaks
pB <- pB + legend_tweaks
if (exists("pC") && !is.null(pC)) pC <- pC + legend_tweaks
if (exists("pD") && !is.null(pD)) pD <- pD + legend_tweaks

# Panel spacing slightly tighter
panel_space <- theme(panel.spacing = unit(5, "pt"),
                     plot.margin   = margin(3, 4, 3, 4))

# ---------------------------
# 4) Compose
# ---------------------------
library(patchwork)

# Main (A–B)
fig_AB <- (pA + pB) +
  plot_layout(ncol = 2) &
  panel_space &
  theme(legend.position = "bottom")

# Extended (A–D), only if both exist
if (exists("pC") && !is.null(pC) && exists("pD") && !is.null(pD)) {
  fig_CD <- (pC + pD)  +
    plot_annotation(tag_levels = "A") &     # <— patchwork adds the letters
    theme(legend.position = "bottom") &
    panel_space
   # theme(legend.position = "bottom")
}

# ---------------------------
# 5) Save
# ---------------------------
outdir <- here::here("outputs", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggsave(here::here(outdir, "Fig4_Traits_vs_CFs.tiff"),
       plot = fig_AB, dpi = 600, width = 185, height = 105,
       units = "mm", device = "tiff", compression = "lzw")

if (exists("fig_CD")) {
  ggsave(here::here(outdir, "Fig4_Traits_vs_CF_EXT.tiff"),
         plot = fig_CD, dpi = 600, width = 185, height = 105,
         units = "mm", device = "tiff", compression = "lzw")
}


ggsave(here::here(outdir, "Fig4_Traits_vs_CF_EXT.png"),
       plot = fig_CD, dpi = 600, width = 185, height = 105,
       units = "mm")

ggsave(here::here(outdir, "Fig4_Traits_vs_CF.png"),
       plot = fig_AB, dpi = 600, width = 185, height = 105,
       units = "mm")


# --- Standardized repo outputs (Fig 4) ---
dir.create(here::here("analysis","figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("manuscript","figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = here::here("analysis","figures","Fig4_panel.png"),
  plot = p_fig, width = 7, height = 5, dpi = 300
)

ggsave(
  filename = here::here("manuscript","figures","Fig4.pdf"),
  plot = p_fig, width = 7, height = 5, device = grDevices::cairo_pdf
)

message("✓ Saved Figure 4: analysis/figures/Fig4_panel.png & manuscript/figures/Fig4.pdf")

