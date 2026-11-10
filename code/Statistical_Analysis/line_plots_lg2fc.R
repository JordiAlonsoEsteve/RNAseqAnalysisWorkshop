# ---------------------------------------------------------------------------
# Gene line plots for a fitted DESeq2 object with a
# treatment * time * temperature design, done with explicit linear algebra.
#
#   cell estimate  =  C %*% beta
#
#   beta : the fitted coefficients of the gene (Intercept, treatment, time,
#          temperature and all interaction terms; log2 scale)
#   C    : one row per treatment x time x temperature cell, one column per
#          coefficient.
#            against = "none"  -> C = X, the model-matrix row of the cell, so
#                                 the estimate is intercept + every effect that
#                                 applies (absolute log2 expression)
#            against = "time" / "treatment" / "temperature"
#                              -> C = X[cell, ] - X[reference cell, ],
#                                 i.e. a log2 fold change
#
#   calc_log2fc()           -> estimates (and SEs) for any genes
#   plot_design_matrices()  -> prints/returns the rows of C used in each plot
#                              panel, restricted to the columns that are used
#   plot_gene_lines()       -> one gene: panel per treatment, line per temperature
#
# Assumptions: dds is fitted, design is ~ treatment * time * temperature,
# temperature levels like "T20", time levels like "t0", treatment any string.
# ---------------------------------------------------------------------------

# ---- small internal helpers -------------------------------------------------
.parse_num <- function(x) {
  suppressWarnings(as.numeric(sub("^[A-Za-z]+", "", as.character(x))))
}

# Sort levels numerically when they carry a number ("t2" < "t10"), else alphabetically
.sort_levels <- function(x) {
  x <- unique(as.character(x))
  n <- .parse_num(x)
  if (anyNA(n)) sort(x) else x[order(n)]
}

# Build X (model matrix of all cells) and C (contrast matrix) ------------------
.contrast_matrix <- function(dds, against, ref,
                             treatment_col, time_col, temperature_col) {
  
  for (pkg in c("DESeq2", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.", call. = FALSE)
  }
  if (!inherits(dds, "DESeqDataSet")) stop("`dds` must be a DESeqDataSet.", call. = FALSE)
  
  mc <- SummarizedExperiment::mcols(SummarizedExperiment::mcols(dds))
  if (is.null(mc$type) || !any(mc$type == "results")) {
    stop("`dds` does not look fitted. Run DESeq2::DESeq(dds) first.", call. = FALSE)
  }
  
  cols <- c(treatment = treatment_col, time = time_col, temperature = temperature_col)
  cd   <- SummarizedExperiment::colData(dds)
  if (any(!cols %in% colnames(cd))) {
    stop("Not in colData(dds): ", paste(setdiff(cols, colnames(cd)), collapse = ", "), call. = FALSE)
  }
  
  # Levels exactly as DESeq2 saw them (alphabetical for character columns);
  # the grid must use these so its model matrix matches the fitted coefficients.
  lv        <- lapply(cols, function(cn) levels(factor(cd[[cn]])))
  lv_sorted <- lapply(lv, .sort_levels)
  
  if (against != "none") {
    if (is.null(ref)) ref <- lv_sorted[[against]][1]
    if (!ref %in% lv[[against]]) {
      stop("`ref` '", ref, "' is not a level of ", cols[[against]], " (",
           paste(lv[[against]], collapse = ", "), ").", call. = FALSE)
    }
  }
  
  # grid of all cells
  grid <- expand.grid(lv$treatment, lv$time, lv$temperature, stringsAsFactors = FALSE)
  names(grid) <- unname(cols[c("treatment", "time", "temperature")])
  for (k in names(cols)) {
    grid[[cols[[k]]]] <- factor(grid[[cols[[k]]]], levels = lv[[k]])
  }
  
  # other design variables (e.g. batch) are held at reference value / mean
  extra <- setdiff(all.vars(DESeq2::design(dds)), cols)
  for (v in extra) {
    x <- cd[[v]]
    grid[[v]] <- if (is.numeric(x)) mean(x) else factor(levels(factor(x))[1], levels = levels(factor(x)))
  }
  if (length(extra) > 0) {
    message("Other design variables held at reference/mean: ", paste(extra, collapse = ", "))
  }
  
  X  <- stats::model.matrix(DESeq2::design(dds), data = grid)
  rn <- DESeq2::resultsNames(dds)
  if (ncol(X) != length(rn)) {
    stop("Design matrix for the cell grid (", ncol(X), " columns) does not match the fitted model (",
         length(rn), " coefficients). Is the design ~ treatment * time * temperature?", call. = FALSE)
  }
  colnames(X) <- rn
  
  cells <- data.frame(
    treatment   = as.character(grid[[cols[["treatment"]]]]),
    time        = as.character(grid[[cols[["time"]]]]),
    temperature = as.character(grid[[cols[["temperature"]]]]),
    stringsAsFactors = FALSE
  )
  labels <- paste(cells$treatment, cells$time, cells$temperature)
  rownames(X) <- labels
  
  if (against == "none") {
    C <- X
  } else {
    ref_cells <- cells
    ref_cells[[against]] <- ref
    j <- match(paste(ref_cells$treatment, ref_cells$time, ref_cells$temperature), labels)
    C <- X - X[j, , drop = FALSE]
  }
  
  list(C = C, X = X, cells = cells, against = against, ref = ref)
}


#' Cell estimates by linear algebra: C %*% beta
#'
#' @param dds     Fitted DESeqDataSet.
#' @param genes   Character vector of rownames(dds).
#' @param against "none" (default): model-estimated log2 expression of each cell
#'                (intercept + all applicable effects). "time", "treatment" or
#'                "temperature": log2 fold change vs the reference level `ref`
#'                of that factor, holding the other two factors fixed.
#' @param ref     Reference level for `against` (default: lowest time / lowest
#'                temperature / first treatment level; set it for treatment,
#'                e.g. "control").
#' @param treatment_col,time_col,temperature_col Column names in colData(dds).
#' @return Long data frame: gene, treatment, time, temperature, time_num,
#'   temperature_num, log2FC, lfcSE. (For against = "none" the `log2FC` column
#'   is the absolute log2 expression estimate.) The contrast matrix is attached
#'   as attr "contrast_matrix", for plot_design_matrices().
#'
#' @examples
#' fc <- calc_log2fc(dds, c("geneA", "geneB"))
#' fc <- calc_log2fc(dds, "geneA", against = "treatment", ref = "control")
calc_log2fc <- function(dds,
                        genes,
                        against         = c("none", "time", "treatment", "temperature"),
                        ref             = NULL,
                        treatment_col   = "treatment",
                        time_col        = "time",
                        temperature_col = "temperature") {
  
  against <- match.arg(against)
  cm <- .contrast_matrix(dds, against, ref, treatment_col, time_col, temperature_col)
  C <- cm$C; cells <- cm$cells
  
  genes <- as.character(genes)
  if (any(!genes %in% rownames(dds))) {
    stop("Genes not found in dds: ", paste(setdiff(genes, rownames(dds)), collapse = ", "), call. = FALSE)
  }
  
  # ---- beta: genes x coefficients (log2 scale, stored in mcols(dds)) ----------
  rn <- colnames(C)
  beta_all <- as.matrix(as.data.frame(SummarizedExperiment::mcols(dds), optional = TRUE)[, rn, drop = FALSE])
  beta <- beta_all[match(genes, rownames(dds)), , drop = FALSE]
  rownames(beta) <- genes
  
  # ---- estimates: (cells x p) %*% (p x genes) = cells x genes ------------------
  est <- C %*% t(beta)
  
  # ---- standard errors of C %*% beta (needs the coefficient covariance, so
  #      taken from DESeq2::results() with the same contrast vector) -------------
  se  <- matrix(NA_real_, nrow(C), length(genes), dimnames = dimnames(est))
  sub <- dds[genes, ]
  for (i in seq_len(nrow(C))) {
    if (all(C[i, ] == 0)) next  # reference cell: estimate is exactly 0
    r <- DESeq2::results(sub, contrast = unname(C[i, ]),
                         independentFiltering = FALSE, cooksCutoff = FALSE)
    se[i, ] <- r$lfcSE
    if (!isTRUE(all.equal(unname(est[i, ]), r$log2FoldChange, tolerance = 1e-6))) {
      warning("C %*% beta disagrees with DESeq2::results() for cell '", rownames(C)[i], "'.",
              call. = FALSE)
    }
  }
  
  n   <- nrow(cells)
  res <- data.frame(
    gene = rep(genes, each = n),
    cells[rep(seq_len(n), times = length(genes)), , drop = FALSE],
    log2FC = as.vector(est),
    lfcSE  = as.vector(se),
    stringsAsFactors = FALSE
  )
  rownames(res) <- NULL
  res$time_num        <- .parse_num(res$time)
  res$temperature_num <- .parse_num(res$temperature)
  res <- res[, c("gene", "treatment", "time", "temperature",
                 "time_num", "temperature_num", "log2FC", "lfcSE")]
  
  attr(res, "against")         <- against
  attr(res, "ref")             <- cm$ref
  attr(res, "contrast_matrix") <- C
  attr(res, "cells")           <- cells
  res
}


#' Design (contrast) matrices behind each plot panel
#'
#' For every treatment panel, returns the rows of C (one per time x temperature
#' cell, ordered by temperature then time) restricted to the columns that are
#' non-zero in that panel, i.e. the coefficients that contribute to it. With
#' against = "none", C is the plain model matrix X.
#'
#' @param x A DESeqDataSet, or the output of calc_log2fc().
#' @param against,ref,treatment_col,time_col,temperature_col As in calc_log2fc()
#'   (only used when `x` is a DESeqDataSet).
#' @param print Print the selected columns and the matrix for each panel.
#' @return Invisibly, a named list (one entry per treatment) of matrices.
#'
#' @examples
#' plot_design_matrices(dds)
#' plot_design_matrices(dds, against = "treatment", ref = "control")
plot_design_matrices <- function(x,
                                 against         = c("none", "time", "treatment", "temperature"),
                                 ref             = NULL,
                                 treatment_col   = "treatment",
                                 time_col        = "time",
                                 temperature_col = "temperature",
                                 print           = TRUE) {
  
  if (inherits(x, "DESeqDataSet")) {
    against <- match.arg(against)
    cm <- .contrast_matrix(x, against, ref, treatment_col, time_col, temperature_col)
  } else {
    cm <- list(C = attr(x, "contrast_matrix"), cells = attr(x, "cells"),
               against = attr(x, "against"), ref = attr(x, "ref"))
    if (is.null(cm$C)) stop("`x` must be a DESeqDataSet or the output of calc_log2fc().", call. = FALSE)
  }
  C <- cm$C; cells <- cm$cells
  
  t_order  <- .sort_levels(cells$temperature)
  tm_order <- .sort_levels(cells$time)
  
  out <- list()
  for (tr in unique(cells$treatment)) {
    idx <- which(cells$treatment == tr)
    o   <- order(match(cells$temperature[idx], t_order), match(cells$time[idx], tm_order))
    M   <- C[idx[o], , drop = FALSE]
    sel <- colSums(abs(M)) > 0
    M   <- M[, sel, drop = FALSE]
    out[[tr]] <- M
    
    if (print) {
      cat("\n== Panel: treatment = ", tr, " ==\n", sep = "")
      cat("Selected columns (", ncol(M), " of ", ncol(C), "): ",
          paste(colnames(M), collapse = ", "), "\n", sep = "")
      print(M)
    }
  }
  invisible(out)
}


#' Line plot of one gene: one panel per treatment, one line per temperature
#'
#' @param dds    Fitted DESeqDataSet (not needed if `fc` is supplied).
#' @param gene   Gene name (a rowname of dds); shown in the title.
#' @param fc     Optional precomputed output of calc_log2fc().
#' @param against,ref Passed to calc_log2fc() (ignored if `fc` is given).
#' @param show_se Draw a +/- 1 SE ribbon around each line.
#' @param print_design Print the selected design-matrix columns for each panel
#'               (see plot_design_matrices()). Set FALSE when looping over genes.
#' @param title  Custom title; default is the gene name.
#' @param file   Optional path; if given the plot is also saved.
#' @param width,height,dpi Size (inches) and resolution for `file`.
#' @param treatment_col,time_col,temperature_col Column names in colData(dds).
#' @return A ggplot object.
#'
#' @examples
#' plot_gene_lines(dds, "geneA")
#' plot_gene_lines(dds, "geneA", against = "treatment", ref = "control", file = "geneA.png")
#'
#' fc <- calc_log2fc(dds, top_genes)
#' for (g in top_genes) plot_gene_lines(gene = g, fc = fc, print_design = FALSE,
#'                                      file = paste0(g, ".png"))
plot_gene_lines <- function(dds,
                            gene,
                            fc              = NULL,
                            against         = c("none", "time", "treatment", "temperature"),
                            ref             = NULL,
                            show_se         = TRUE,
                            print_design    = TRUE,
                            title           = NULL,
                            file            = NULL,
                            width           = 9,
                            height          = 4.5,
                            dpi             = 300,
                            treatment_col   = "treatment",
                            time_col        = "time",
                            temperature_col = "temperature") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  if (length(gene) != 1) stop("`gene` must be a single gene name.", call. = FALSE)
  
  if (is.null(fc)) {
    against <- match.arg(against)
    fc <- calc_log2fc(dds, gene, against = against, ref = ref,
                      treatment_col = treatment_col, time_col = time_col,
                      temperature_col = temperature_col)
  }
  against_used <- attr(fc, "against")
  ref_used     <- attr(fc, "ref")
  
  if (print_design) plot_design_matrices(fc)
  
  d <- fc[fc$gene == gene, , drop = FALSE]
  if (nrow(d) == 0) stop("Gene '", gene, "' not found in `fc`.", call. = FALSE)
  
  d$temperature <- factor(d$temperature, levels = .sort_levels(d$temperature))
  d$se <- ifelse(is.na(d$lfcSE), 0, d$lfcSE)
  
  # x axis: numeric time if the labels carry a number, otherwise ordered index
  tl <- .sort_levels(d$time)
  d$x <- if (!anyNA(d$time_num)) d$time_num else match(d$time, tl)
  brk <- sort(unique(d$x))
  lab <- d$time[match(brk, d$x)]
  
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = log2FC,
                                       colour = temperature, fill = temperature,
                                       group = temperature))
  if (show_se) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = log2FC - se, ymax = log2FC + se),
                                  alpha = 0.15, colour = NA)
  }
  if (!identical(against_used, "none")) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")
  }
  p <- p +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ treatment) +
    ggplot2::scale_x_continuous(breaks = brk, labels = lab) +
    ggplot2::scale_colour_viridis_d(option = "plasma", end = 0.85) +
    ggplot2::scale_fill_viridis_d(option = "plasma", end = 0.85) +
    ggplot2::labs(
      title  = if (is.null(title)) gene else title,
      x      = "Time",
      y      = if (identical(against_used, "none")) {
        "Model-estimated log2 expression"
      } else {
        paste0("log2 fold change vs ", against_used, " = ", ref_used)
      },
      colour = "Temperature", fill = "Temperature"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey92"))
  
  if (!is.null(file)) {
    ggplot2::ggsave(file, p, width = width, height = height, dpi = dpi)
  }
  p
}