# ============================================================
# 01_setup.R — Common R setup for necromass-cf
# ============================================================

cat("\n=== necromass-cf: R setup ===\n")

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  dplyr.summarise.inform = FALSE
)

pkgs <- c(
  "tidyverse",
  "ggridges",
  "patchwork",
  "forcats",
  "here"
)

# Load required packages (assumed installed via renv or manually)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  warning("Missing packages (install via renv or install.packages): ",
          paste(missing, collapse = ", "))
} else {
  cat("Loaded packages:\n  ", paste(pkgs, collapse = ", "), "\n")
}

suppressPackageStartupMessages({
  lapply(pkgs, require, character.only = TRUE)
})

# A consistent base theme for figures
theme_set(ggplot2::theme_minimal(base_size = 9))

cat("R version:\n")
print(R.version.string)
cat("=== setup complete ===\n")
