# geomorph_helpers.R
# Standalone functions for CSV → TPS conversion and geomorph analyses.
# Source this file in any R session or Quarto document.
#
# Usage:
#   source("geomorph_helpers.R")
#   df  <- read_landmark_csv("my_specimen.csv")
#   lm  <- parse_landmarks(df)
#   write_tps(lm, "Genus_species", "output.tps")
#   res <- run_geomorph(lm)

library(readr)
library(dplyr)
library(stringr)

# ── 1. Read CSV ────────────────────────────────────────────────────────────────

#' Read a landmark CSV with automatic delimiter detection
#' @param path Path to CSV file
#' @return data.frame
read_landmark_csv <- function(path) {
  df <- tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) read_delim(path, delim = ";", show_col_types = FALSE)
  )
  df
}

# ── 2. Parse & classify landmarks ─────────────────────────────────────────────

#' Parse a raw landmark data frame into a clean form
#'
#' Your CSV has columns: label, r, a, s (RAS coordinates from 3D Slicer).
#' Semilandmarks are identified by label prefix "OC_" or empty/NA labels.
#'
#' @param df           Raw data frame from read_landmark_csv()
#' @param label_col    Column name for landmark labels  (default "label")
#' @param x_col        Column name for X / R coordinate (default "r")
#' @param y_col        Column name for Y / A coordinate (default "a")
#' @param z_col        Column name for Z / S coordinate (default "s")
#' @param semi_pattern Regex that flags a label as a semilandmark (default "^OC_")
#' @return Tibble with columns: label, x, y, z, type ("fixed"|"semi")
parse_landmarks <- function(df,
                             label_col    = "label",
                             x_col        = "r",
                             y_col        = "a",
                             z_col        = "s",
                             semi_pattern = "^OC_") {
  df <- df |>
    select(label = all_of(label_col),
           x     = all_of(x_col),
           y     = all_of(y_col),
           z     = all_of(z_col)) |>
    mutate(
      label = as.character(label),
      label = if_else(is.na(label) | label == "",
                      paste0("unlabeled_", row_number()), label),
      x = as.numeric(x),
      y = as.numeric(y),
      z = as.numeric(z),
      type = case_when(
        str_detect(label, semi_pattern)  ~ "semi",
        str_detect(label, "^unlabeled_") ~ "semi",
        TRUE                              ~ "fixed"
      )
    ) |>
    filter(!is.na(x), !is.na(y), !is.na(z))
  df
}

# ── 3. Build sliders matrix ────────────────────────────────────────────────────

#' Build a geomorph-compatible sliders matrix from parsed landmarks
#'
#' geomorph expects a 3-column matrix where each row is
#' [before_index, semi_index, after_index].
#' Consecutive semilandmarks are grouped into sliding triplets.
#'
#' @param lm_df  Output of parse_landmarks()
#' @return Matrix (n_semi-2) × 3, or NULL if fewer than 3 semilandmarks
build_sliders <- function(lm_df) {
  semi_idx <- which(lm_df$type == "semi")
  if (length(semi_idx) < 3) {
    message("Fewer than 3 semilandmarks found — no sliders built.")
    return(NULL)
  }
  n <- length(semi_idx)
  sliders <- matrix(NA_integer_, nrow = n - 2, ncol = 3)
  for (i in 2:(n - 1)) {
    sliders[i - 1, ] <- c(semi_idx[i - 1], semi_idx[i], semi_idx[i + 1])
  }
  sliders[complete.cases(sliders), , drop = FALSE]
}

# ── 4. Write TPS ───────────────────────────────────────────────────────────────

#' Write landmarks to TPS format
#'
#' TPS format per specimen:
#'   LM3=<n>
#'   x1\ty1\tz1
#'   ...
#'   ID=<specimen_name>
#'
#' @param lm_df         Output of parse_landmarks()
#' @param specimen_name ID string written to the ID= field (e.g. "Genus_species")
#' @param outfile       Path to write .tps file; if NULL returns lines as a vector
#' @return Invisibly returns the character vector of TPS lines
write_tps <- function(lm_df, specimen_name, outfile = NULL) {
  n <- nrow(lm_df)
  fmt <- function(v) formatC(v, format = "f", digits = 4, width = 12, flag = " ")
  coord_lines <- paste(fmt(lm_df$x), fmt(lm_df$y), fmt(lm_df$z), sep = "  ")
  lines <- c(
    paste0("LM3=", n),
    coord_lines,
    paste0("ID=", specimen_name)
  )
  if (!is.null(outfile)) {
    writeLines(lines, outfile)
    message("TPS written to: ", outfile)
  }
  invisible(lines)
}

#' Write multiple specimens to one TPS file
#'
#' @param specimens_list A named list of parse_landmarks() outputs.
#'   Names are used as specimen IDs. e.g. list(Genus_speciesA = lm1, Genus_speciesB = lm2)
#' @param outfile  Output .tps file path
write_tps_multi <- function(specimens_list, outfile) {
  all_lines <- unlist(mapply(
    FUN           = write_tps,
    lm_df         = specimens_list,
    specimen_name = names(specimens_list),
    outfile       = NULL,
    SIMPLIFY      = FALSE
  ))
  writeLines(all_lines, outfile)
  message("Multi-specimen TPS written to: ", outfile)
}

# ── 5. Run geomorph GPA + PCA ─────────────────────────────────────────────────

#' Run GPA and PCA on a single (or multi-specimen) landmark set
#'
#' @param lm_df         Output of parse_landmarks()
#' @param specimen_name Specimen name (used as dimname)
#' @param use_sliders   Logical; if TRUE and semilandmarks exist, use sliding
#' @param procD         Procrustes distance: "full" or "orthogonal"
#' @return List with elements: gpa, pca, lm_df, sliders
run_geomorph <- function(lm_df,
                          specimen_name = "Specimen_1",
                          use_sliders   = TRUE,
                          procD         = "full") {
  if (!requireNamespace("geomorph", quietly = TRUE))
    stop("Package 'geomorph' is required. Install with: install.packages('geomorph')")
  library(geomorph)

  p <- nrow(lm_df)
  coords_matrix <- as.matrix(lm_df[, c("x", "y", "z")])
  A <- array(coords_matrix, dim = c(p, 3, 1))
  dimnames(A)[[1]] <- lm_df$label
  dimnames(A)[[3]] <- specimen_name

  sliders <- NULL
  if (use_sliders && any(lm_df$type == "semi")) {
    sliders <- build_sliders(lm_df)
    if (!is.null(sliders))
      message("Built sliders matrix: ", nrow(sliders), " sliding triplets")
  }

  message("Running GPA...")
  gpa <- gpagen(A, curves = sliders, ProcD = (procD == "full"),
                print.progress = FALSE)

  message("Running PCA...")
  pca <- gm.prcomp(gpa$coords)

  list(gpa = gpa, pca = pca, lm_df = lm_df, sliders = sliders)
}

#' Quick summary print for run_geomorph() output
print_geomorph_summary <- function(res) {
  cat("=== GPA summary ===\n")
  print(summary(res$gpa))

  eigs    <- res$pca$sdev^2
  var_exp <- round(eigs / sum(eigs) * 100, 2)
  cum_var <- cumsum(var_exp)
  cat("\n=== PCA variance explained ===\n")
  pca_tbl <- data.frame(
    PC             = paste0("PC", seq_along(var_exp)),
    Variance_pct   = var_exp,
    Cumulative_pct = cum_var
  )
  print(pca_tbl)
}