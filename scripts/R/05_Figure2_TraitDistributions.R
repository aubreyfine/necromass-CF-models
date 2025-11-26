# ============================================================
# Figure 2 – Variation in microbial amino-sugar content
# Inputs : data/bactMurA.csv, data/fungGlcN.csv
# Output : outputs/figures/Figure2.tiff
# ============================================================

options(scipen = 999)
suppressPackageStartupMessages({
  library(tidyverse)
  library(forcats)
  library(ggridges)
  library(patchwork)
  library(here)
})

# ---- check project root ----
message("Project root: ", here::here())

# ---- paths ----
bact_file <- here::here("data", "bactMurA.csv")
fung_file <- here::here("data", "fungGlcN.csv")
outdir    <- here::here("outputs", "figures")
outfile   <- here::here(outdir, "Figure2.tiff")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- palettes and phyla ----
pal_gram <- c(GN = "#43BBAD", GP = "#414388FF")
gram_map <- c("Pseudomonadota" = "GN", "Firmicutes" = "GP", "Actinomycetota" = "GP")
phyl_A   <- c("Pseudomonadota","Firmicutes","Actinomycetota")
phyl_B   <- c("Ascomycota","Basidiomycota")

# ---- read data ----
bac <- readr::read_csv(bact_file, show_col_types = FALSE) %>%
  dplyr::rename(
    Phylum = tidyselect::any_of(c("Phylum","phylum","BacterialPhylum")),
    MurA   = tidyselect::any_of(c("MurA","murA","MurA_mg_g_biomass"))
  ) %>%
  dplyr::select(Phylum, MurA)

fung <- readr::read_csv(fung_file, show_col_types = FALSE) %>%
  dplyr::rename(
    Phylum = tidyselect::any_of(c("Phylum","phylum","FungalPhylum")),
    GlcN   = tidyselect::any_of(c("GlcN","glcN","GlcN_mg_g_biomass"))
  ) %>%
  dplyr::select(Phylum, GlcN)

# ---- PANEL A — Bacteria ----
bac_dat <- bac %>%
  dplyr::filter(Phylum %in% phyl_A) %>%
  dplyr::transmute(
    Level = factor(Phylum, levels = phyl_A),
    MurA  = as.numeric(MurA),
    Gram  = factor(dplyr::recode(Phylum, !!!gram_map), levels = c("GN","GP"))
  )

nA <- bac_dat %>% dplyr::count(Level, name = "nA")
ylabs_A <- setNames(paste0(nA$Level, " (N=", nA$nA, ")"), nA$Level)

pA <- ggplot(bac_dat, aes(x = MurA, y = Level, fill = Gram)) +
  stat_density_ridges(
    jittered_points = TRUE, point_shape = "|", point_size = 2.5,
    point_alpha = 0.45, point_color = "grey30",
    rel_min_height = 0.001, scale = 0.95,
    quantiles = 4, quantile_lines = TRUE,
    vline_linetype = 1, vline_size = 0.35,
    alpha = 0.65,
    position = position_points_jitter(width = 0.05, height = 0)
  ) +
  scale_fill_manual(values = pal_gram, name = " ",
                    labels = c("Gram-negative", "Gram-positive")) +
  scale_y_discrete(labels = ylabs_A,
                   expand = expansion(mult = c(0.05, 0.04))) +
  labs(tag = "A",
       y = "Bacterial phylum",
       x = expression("MurA (mg" ~ g^{-1} ~ "bacterial biomass)")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 6)),
        axis.title.y = element_text(margin = margin(r = 6)))

# ---- PANEL B — Fungi ----
fung_dat <- fung %>%
  dplyr::filter(Phylum %in% phyl_B) %>%
  dplyr::transmute(
    Level = factor(Phylum, levels = phyl_B),
    GlcN  = as.numeric(GlcN)
  )

nB <- fung_dat %>% dplyr::count(Level, name = "nB")
ylabs_B <- setNames(paste0(nB$Level, " (N=", nB$nB, ")"), nB$Level)

pB <- ggplot(fung_dat, aes(x = GlcN, y = forcats::fct_rev(Level))) +
  stat_density_ridges(
    jittered_points = TRUE, point_shape = "|", point_size = 2.5,
    point_alpha = 0.45, point_color = "grey30",
    rel_min_height = 0.001, scale = 0.95,
    fill = "#3482A4FF", alpha = 0.65,
    quantiles = 4, quantile_lines = TRUE,
    vline_linetype = 1, vline_size = 0.35,
    position = position_points_jitter(width = 0.05, height = 0)
  ) +
  scale_y_discrete(labels = function(l) ylabs_B[as.character(l)],
                   expand  = expansion(mult = c(0.05, 0.04))) +
  labs(tag = "B",
       y = "Fungal phylum",
       x = expression("GlcN (mg" ~ g^{-1} ~ "fungal biomass)")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 6)),
        axis.title.y = element_text(margin = margin(r = 6)))

# ---- Combine & save ----
p_fig <- pA + pB + patchwork::plot_layout(heights = c(3, 2))

ggsave(filename = outfile, plot = p_fig,
       width = 185, height = 185, units = "mm",
       dpi = 500, compression = "lzw")
