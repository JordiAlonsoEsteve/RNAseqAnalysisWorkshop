library(clusterProfiler)
library(org.At.tair.db)
library(enrichplot)
library(ggplot2)

#' Perform GO Enrichment Analysis for Arabidopsis AGI Genes
#'
#' @param gene_list Character vector of Arabidopsis gene IDs (e.g. c("AT1G01010", "AT2G02020"))
#' @param ont Character string: "BP" (Biological Process, default), "MF" (Molecular Function), "CC" (Cellular Component), or "ALL"
#' @param pvalueCutoff Numeric: Adjusted p-value cutoff (default 0.05)
#' @param qvalueCutoff Numeric: q-value cutoff (default 0.2)
#' @param universe Optional character vector of background gene IDs
#' @return A list containing the S4 enrichResult object, a data frame table, and a ggplot2 dotplot
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