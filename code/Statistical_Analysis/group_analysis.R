library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)


get_deg_sets <- function(dds, computed_obj,
                         padj_cutoff   = 0.05,
                         lfc_cutoff    = 1,       # |log2FC| >= this (0 = no LFC filter)
                         min_base_mean = 0) {     # optional abundance filter
  
  grid <- computed_obj$grid
  rl   <- computed_obj$results_list
  
  clean <- function(x) sub("^[A-Za-z](?=[-0-9.])", "", x, perl = TRUE)  # "t0.5" -> "0.5"
  keys  <- paste(grid$Temperature, grid$Treatment, grid$Time, sep = "__")
  
  sets    <- list()
  summ    <- list()
  details <- list()
  
  for (i in seq_len(nrow(grid))) {
    res <- rl[[keys[i]]]
    if (is.null(res)) next                       # reference (mock) rows / failed contrasts
    
    res <- res[rownames(res) %in% rownames(dds), ]
    
    sig  <- !is.na(res$padj) & res$padj < padj_cutoff & res$baseMean >= min_base_mean
    up   <- sig & res$log2FoldChange >=  lfc_cutoff
    down <- sig & res$log2FoldChange <= -lfc_cutoff
    
    nm <- sprintf("%s_%sC_%sh", grid$Treatment[i], clean(grid$Temperature[i]), clean(grid$Time[i]))
    
    sets[[nm]] <- list(up   = rownames(res)[up],
                       down = rownames(res)[down])
    
    summ[[nm]] <- data.frame(
      contrast    = nm,
      Treatment   = grid$Treatment[i],
      Temperature = clean(grid$Temperature[i]),
      Time        = clean(grid$Time[i]),
      n_tested    = sum(!is.na(res$padj)),
      n_up        = sum(up),
      n_down      = sum(down),
      stringsAsFactors = FALSE
    )
    
    keep <- up | down
    if (any(keep)) {
      details[[nm]] <- data.frame(
        contrast    = nm,
        Treatment   = grid$Treatment[i],
        Temperature = clean(grid$Temperature[i]),
        Time        = clean(grid$Time[i]),
        gene        = rownames(res)[keep],
        baseMean    = res$baseMean[keep],
        log2FC      = res$log2FoldChange[keep],
        padj        = res$padj[keep],
        direction   = ifelse(res$log2FoldChange[keep] > 0, "up", "down"),
        stringsAsFactors = FALSE
      )
    }
  }
  
  list(
    sets    = sets,                                        # sets[["ABA_28C_2h"]]$up / $down
    summary = do.call(rbind, c(summ, make.row.names = FALSE)),   # one row per contrast
    table   = do.call(rbind, c(details, make.row.names = FALSE)) # long table, one row per DEG
  )
}


run_arabidopsis_go <- function(gene_list, 
                               ont = "BP", 
                               pvalueCutoff = 0.05, 
                               qvalueCutoff = 0.2, 
                               universe = NULL) {
  
  # Validate packages
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.At.tair.db", quietly = TRUE)) {
    stop("Please install 'clusterProfiler' and 'org.At.tair.db' via Bioconductor.")
  }
  
  message("Running GO enrichment for ", length(gene_list), " genes (Ontology: ", ont, ")...")
  
  # Run enrichGO using "TAIR" keyType for AGI identifiers (AT1G..., AT2G..., etc.)
  ego <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = ont,
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    universe      = universe,
    readable      = FALSE
  )
  
  # Handle empty results
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning("No enriched GO terms found matching the thresholds.")
    return(list(enrichResult = ego, table = data.frame(), plot = NULL))
  }
  
  # Extract results table
  ego_df <- as.data.frame(ego)
  
  # Generate standard dotplot
  p <- dotplot(ego, showCategory = 20, title = paste("GO Enrichment Analysis -", ont)) +
    theme_bw(base_size = 12)
  
  # Return bundle
  return(list(
    enrichResult = ego,
    table        = ego_df,
    plot         = p
  ))
}