#' Generate volcano plots from a fitted DESeq2 object
#'
#' Takes a dds that has already been through DESeq() and writes one volcano
#' plot (PNG + PDF) and one results table (CSV) per comparison into a folder,
#' with the top n genes labelled.
#'
#' @param dds          A fitted DESeqDataSet (DESeq() already run).
#' @param out_dir      Folder to create/write into.
#' @param n_top        Number of top genes to label per plot.
#' @param coefs        Coefficients to plot (character vector from resultsNames(dds)).
#'                     Default: all except "Intercept".
#' @param contrasts    Optional list of contrasts, each c(variable, numerator, denominator).
#'                     If given, these are plotted in addition to `coefs`
#'                     (set coefs = character(0) to plot only contrasts).
#' @param rank_by      "padj" (default; ties broken by |log2FC|) or "lfc"
#'                     (largest |log2FC| among significant genes).
#' @param padj_cutoff  Adjusted p-value threshold.
#' @param lfc_cutoff   |log2 fold change| threshold.
#' @param label_col    Column in rowData(dds) holding gene symbols; NULL = rownames.
#' @param shrink       Apply lfcShrink (type = "ashr") to the log2FC values.
#' @param width,height Plot size in inches.
#' @param dpi          PNG resolution.
#' @return Invisibly, a named list of results data frames.
#'
#' @examples
#' dds <- DESeq2::DESeq(dds)
#' volcano_plots_from_dds(dds, out_dir = "volcano_plots", n_top = 15)
volcano_plots_from_dds <- function(dds,
                                   out_dir     = "volcano_plots",
                                   n_top       = 10,
                                   coefs       = NULL,
                                   contrasts   = NULL,
                                   rank_by     = c("padj", "lfc"),
                                   padj_cutoff = 0.05,
                                   lfc_cutoff  = 1,
                                   label_col   = NULL,
                                   shrink      = FALSE,
                                   width       = 7,
                                   height      = 6,
                                   dpi         = 300) {

  rank_by <- match.arg(rank_by)

  for (pkg in c("DESeq2", "SummarizedExperiment", "ggplot2", "ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }
  if (!inherits(dds, "DESeqDataSet")) stop("`dds` must be a DESeqDataSet.", call. = FALSE)

  mc <- SummarizedExperiment::mcols(SummarizedExperiment::mcols(dds))
  if (is.null(mc$type) || !any(mc$type == "results")) {
    stop("`dds` does not look fitted. Run DESeq2::DESeq(dds) first.", call. = FALSE)
  }

  # ---- Work out what to plot -----------------------------------------------
  if (is.null(coefs)) {
    coefs <- setdiff(DESeq2::resultsNames(dds), "Intercept")
  }
  jobs <- list()
  for (cf in coefs) {
    jobs[[cf]] <- list(name = cf, coef = cf, contrast = NULL, title = cf)
  }
  for (ctr in contrasts) {
    nm <- paste0(ctr[1], "_", ctr[2], "_vs_", ctr[3])
    jobs[[nm]] <- list(name = nm, coef = NULL, contrast = ctr,
                       title = paste0(ctr[1], ": ", ctr[2], " vs ", ctr[3]))
  }
  if (length(jobs) == 0) stop("Nothing to plot: no coefficients or contrasts.", call. = FALSE)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Gene labels -----------------------------------------------------------
  rd <- SummarizedExperiment::rowData(dds)
  gene_labels <- if (!is.null(label_col) && label_col %in% colnames(rd)) {
    lab <- as.character(rd[[label_col]])
    ifelse(is.na(lab) | lab == "", rownames(dds), lab)
  } else {
    rownames(dds)
  }
  names(gene_labels) <- rownames(dds)

  all_results <- list()

  for (job in jobs) {
    message("Processing: ", job$name)

    res <- if (shrink) {
      if (!requireNamespace("ashr", quietly = TRUE)) {
        stop("shrink = TRUE needs the 'ashr' package.", call. = FALSE)
      }
      if (!is.null(job$coef)) {
        DESeq2::lfcShrink(dds, coef = job$coef, type = "ashr", quiet = TRUE)
      } else {
        DESeq2::lfcShrink(dds, contrast = job$contrast, type = "ashr", quiet = TRUE)
      }
    } else {
      if (!is.null(job$coef)) {
        DESeq2::results(dds, name = job$coef)
      } else {
        DESeq2::results(dds, contrast = job$contrast)
      }
    }

    df <- as.data.frame(res)
    df$gene  <- rownames(df)
    df$label <- gene_labels[df$gene]
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), , drop = FALSE]

    # padj can underflow to 0; cap -log10 at the smallest non-zero value
    min_p <- min(df$padj[df$padj > 0], 1e-300)
    df$neglog10padj <- -log10(pmax(df$padj, min_p))

    df$status <- "NS"
    df$status[df$padj < padj_cutoff & df$log2FoldChange >=  lfc_cutoff] <- "Up"
    df$status[df$padj < padj_cutoff & df$log2FoldChange <= -lfc_cutoff] <- "Down"
    df$status <- factor(df$status, levels = c("Up", "Down", "NS"))

    # ---- Top n genes among significant ones -----------------------------------
    sig <- df[df$status != "NS", , drop = FALSE]
    if (nrow(sig) > 0 && n_top > 0) {
      ord <- if (rank_by == "padj") {
        order(sig$padj, -abs(sig$log2FoldChange))
      } else {
        order(-abs(sig$log2FoldChange), sig$padj)
      }
      top_genes <- sig$gene[ord][seq_len(min(n_top, nrow(sig)))]
    } else {
      top_genes <- character(0)
    }
    df$is_top <- df$gene %in% top_genes

    n_up   <- sum(df$status == "Up")
    n_down <- sum(df$status == "Down")

    p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neglog10padj)) +
      ggplot2::geom_point(ggplot2::aes(colour = status), alpha = 0.6, size = 1.3) +
      ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
                          linetype = "dashed", colour = "grey40") +
      ggplot2::geom_hline(yintercept = -log10(padj_cutoff),
                          linetype = "dashed", colour = "grey40") +
      ggrepel::geom_text_repel(
        data = df[df$is_top, , drop = FALSE],
        ggplot2::aes(label = label),
        size = 3, max.overlaps = Inf, box.padding = 0.4,
        segment.colour = "grey30", min.segment.length = 0
      ) +
      ggplot2::scale_colour_manual(
        values = c(Up = "#D55E00", Down = "#0072B2", NS = "grey70"),
        labels = c(Up   = paste0("Up (", n_up, ")"),
                   Down = paste0("Down (", n_down, ")"),
                   NS   = "NS"),
        drop = FALSE
      ) +
      ggplot2::labs(
        title    = job$title,
        subtitle = paste0("padj < ", padj_cutoff, ", |log2FC| \u2265 ", lfc_cutoff,
                          "; top ", length(top_genes), " labelled"),
        x = if (shrink) "log2 fold change (shrunken)" else "log2 fold change",
        y = expression(-log[10]~"adjusted p-value"),
        colour = NULL
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     legend.position = "top")

    ggplot2::ggsave(file.path(out_dir, paste0("volcano_", job$name, ".png")),
                    p, width = width, height = height, dpi = dpi)
    ggplot2::ggsave(file.path(out_dir, paste0("volcano_", job$name, ".pdf")),
                    p, width = width, height = height)

    out_tab <- df[order(df$padj), c("gene", "label", "baseMean", "log2FoldChange",
                                    "padj", "status")]
    utils::write.csv(out_tab, file.path(out_dir, paste0("results_", job$name, ".csv")),
                     row.names = FALSE)

    all_results[[job$name]] <- out_tab
  }

  message("Done. Output written to: ", normalizePath(out_dir))
  invisible(all_results)
}