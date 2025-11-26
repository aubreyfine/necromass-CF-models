# ============================================================
# Figure 3 — Global sensitivity (Sobol) of the CF pipeline
# From-scratch build matching your second script’s style.
# - Reads four CSVs from your Documents path
# - Canonicalizes parameter names (so 'SOC', 'Soil C', etc. → 'S')
# - Ensures numeric columns are numeric
# - Data-driven top-3 highlight per panel
# - Safe expression labels (S always shown)
# - Saves a 600 dpi LZW-compressed TIFF
# ============================================================

options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ---------- INPUTS ----------
in_dir <- "C:/Users/xpf/Documents/CF_update/sobol_outputs"

sobol_cfb          <- read.csv(file.path(in_dir, "sobol_CFB_bootstrap.csv"), check.names = FALSE)
sobol_nf           <- read.csv(file.path(in_dir, "sobol_NF_bootstrap.csv"), check.names = FALSE)
sobol_fnecC_agg    <- read.csv(file.path(in_dir, "sobol_fnecC_agg_bootstrap.csv"), check.names = FALSE)
sobol_fnecC_int    <- read.csv(file.path(in_dir, "sobol_fnecC_int_bootstrap.csv"), check.names = FALSE)

# ---------- HELPERS ----------
canon_param <- function(x){
  x <- trimws(as.character(x))
  x <- gsub("\u200B", "", x)      # strip zero-width spaces
  x <- gsub("\\s+", "", x)        # remove internal spaces
  lx <- tolower(x)
  x[lx %in% c("soc","soilc","totalc","s_percent","soc_pct","soil_c","soilcarbon","s")] <- "S"
  x[lx %in% c("gs","glcn","glucosamine","soilglcn","soilglucosamine")] <- "GS"
  x[lx %in% c("ms","mura","soilmura")] <- "MS"
  x[lx %in% c("cff","cf_fungal","cf_fun")] <- "CFF"
  x[lx %in% c("cfb","cf_bac","cf_bacterial")] <- "CFB"
  x[lx %in% c("cb","c_b")] <- "cB"
  x[lx %in% c("mgp","m_gp")] <- "MGP"
  x[lx %in% c("mgn","m_gn")] <- "MGN"
  x[lx %in% c("fgp","f_gp")] <- "fGP"
  x[lx %in% c("cf","c_f")] <- "cF"
  x[lx %in% c("gf","g_f")] <- "GF"
  x[lx %in% c("rb","r_b")] <- "rB"
  x
}

coerce_numeric_cols <- function(df){
  num_cols <- intersect(names(df), c("S1_median","S1_lo95","S1_hi95","ST_median","ST_lo95","ST_hi95"))
  df[num_cols] <- lapply(df[num_cols], function(z) as.numeric(as.character(z)))
  df
}

label_with_fallback <- function(keys, expr_map) {
  out <- vector("list", length(keys))
  for (i in seq_along(keys)) {
    k <- keys[i]
    out[[i]] <- if (!is.null(expr_map[[k]])) expr_map[[k]] else as.expression(k)
  }
  as.expression(out)
}

# --- Helper: build one panel (keeps your look 1:1) ---
build_panel <- function(df, label_map) {
  df %>%
    arrange(desc(ST_median)) %>%
    mutate(
      param_ord = factor(parameter, levels = parameter),
      # data-driven top3 by ST_median:
      top3 = if_else(parameter %in% head(parameter[order(-ST_median)], 3), "High", "Low")
    ) %>%
    { 
      x_min <- min(.$S1_lo95, .$ST_lo95, 0, na.rm = TRUE)
      x_max <- max(.$S1_hi95, .$ST_hi95, 1, na.rm = TRUE)
      ggplot(., aes(y = param_ord)) +
        geom_errorbarh(aes(xmin = ST_lo95, xmax = ST_hi95),
                       height = 0.55, color = "grey50", linewidth = 0.5, alpha = 0.6) +
        geom_col(aes(x = ST_median, fill = top3), width = 0.6, color = "grey30") +
        scale_fill_manual(values = c("High" = "#43BBADFF", "Low" = "grey80"), guide = "none") +
        geom_errorbarh(aes(xmin = S1_lo95, xmax = S1_hi95),
                       height = 0.22, color = "black", linewidth = 0.6) +
        geom_point(aes(x = S1_median), size = 2, shape = 21, fill = "black") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
        scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = c(0.01, 0.05))) +
        scale_y_discrete(labels = label_map) +
        labs(x = "Sobol index", y = NULL) +
        theme_classic(base_size = 8) +
        theme(
          plot.title.position = "plot",
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "plain", size = 11, hjust = 0)
        )
    }
}

# --- Use it for each panel (your same label sets) ---
labels_cfb <- c(
  cB  = expression(italic(c)[B]),
  MGP = expression(italic(M)[GP]),
  MGN = expression(italic(M)[GN]),
  fGP = expression(italic(f)[GP])
)
labels_nf <- c(
  MS  = expression(italic(M)[S]),
  GS  = expression(italic(G)[S]),
  cF  = expression(italic(c)[F]),
  GF  = expression(italic(G)[F]),
  rB  = expression(italic(r)[B])
)
labels_agg <- c(
  MS  = expression(italic(M)[S]),
  GS  = expression(italic(G)[S]),
  S   = expression(italic(S)),
  CFF = expression(italic(CF)[F]),
  CFB = expression(italic(CF)[B])
)
labels_int <- c(
  MS  = expression(italic(M)[S]),
  GS  = expression(italic(G)[S]),
  S   = expression(italic(S)),
  cB  = expression(italic(c)[B]),
  MGP = expression(italic(M)[GP]),
  MGN = expression(italic(M)[GN]),
  fGP = expression(italic(f)[GP]),
  cF  = expression(italic(c)[F]),
  GF  = expression(italic(G)[F]),
  rB  = expression(italic(r)[B])
)

p_sobol_cfb <- build_panel(sobol_cfb,            labels_cfb)
p_sobol_NF  <- build_panel(sobol_nf,             labels_nf)
p_sobol_agg <- build_panel(sobol_fnecC_agg,      labels_agg)
p_sobol_int <- build_panel(sobol_fnecC_int, labels_int)

design <- "
AADDD
BBDDD
CCDDD
"
plots <- p_sobol_cfb + p_sobol_NF + p_sobol_agg + p_sobol_int +
  plot_layout(design = design) +
  plot_annotation(tag_levels = "A")
plots
ggsave("Figure3_sobolPlots.tiff", dpi = 600,
       width = 185, height = 160, units = "mm",
       compression = "lzw")

# --- Standardized repo outputs (Fig 3) ---
dir.create(here::here("analysis","figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("manuscript","figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = here::here("analysis","figures","Fig3_panel.png"),
  plot = p_fig, width = 7, height = 5, dpi = 300
)

ggsave(
  filename = here::here("manuscript","figures","Fig3.pdf"),
  plot = p_fig, width = 7, height = 5, device = grDevices::cairo_pdf
)

message("✓ Saved Figure 3: analysis/figures/Fig3_panel.png & manuscript/figures/Fig3.pdf")

