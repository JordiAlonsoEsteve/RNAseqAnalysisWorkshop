library(shiny)
library(ggplot2)
library(DT)

build_coldata <- function(file,
                          sep = NULL,
                          id_cols = c("gene_id", "gene_name"),
                          drop_treatments = "untreated",
                          ref_treatment = "mock",
                          ref_temperature = 21,
                          ref_time = NULL) {   # NULL = smallest time left after dropping treatments
  
  # only the header is needed
  header <- read.table(file, sep = "\t", header = TRUE, nrows = 1,
                       check.names = FALSE, quote = "\"", comment.char = "")
  samples <- setdiff(colnames(header), id_cols)
  
  # parse <Temperature>_<Treatment>_<Time>_<Replicate>
  pat   <- "^(-?[0-9.]+)_(.+)_(-?[0-9.]+)_([0-9]+)$"
  parts <- regmatches(samples, regexec(pat, samples))
  bad   <- lengths(parts) != 5
  if (any(bad)) {
    stop("Column names not in <Temperature>_<Treatment>_<Time>_<Replicate> form: ",
         paste(head(samples[bad], 10), collapse = ", "))
  }
  
  # the replicate (parts[[i]][5]) is deliberately not stored
  cd <- data.frame(
    temperature = as.numeric(sapply(parts, `[`, 2)),
    treatment   = sapply(parts, `[`, 3),
    time        = as.numeric(sapply(parts, `[`, 4)),
    row.names   = samples,
    stringsAsFactors = FALSE
  )
  
  # drop unwanted treatments (e.g. untreated, which has no time-matched mock)
  if (length(drop_treatments)) cd <- cd[!(cd$treatment %in% drop_treatments), , drop = FALSE]
  
  # resolve and validate the references
  if (is.null(ref_time)) ref_time <- min(cd$time)
  if (!(ref_treatment   %in% cd$treatment))   stop("ref_treatment '", ref_treatment, "' not found. Available: ",
                                                   paste(unique(cd$treatment), collapse = ", "))
  if (!(ref_temperature %in% cd$temperature)) stop("ref_temperature ", ref_temperature, " not found. Available: ",
                                                   paste(sort(unique(cd$temperature)), collapse = ", "))
  if (!(ref_time        %in% cd$time))        stop("ref_time ", ref_time, " not found. Available: ",
                                                   paste(sort(unique(cd$time)), collapse = ", "))
  
  # factors: numeric ordering, with the chosen reference level first
  time_levels <- make.names(paste0("t", sort(unique(cd$time))))
  temp_levels <- make.names(paste0("T", sort(unique(cd$temperature))))
  
  cd$treatment   <- relevel(factor(make.names(cd$treatment)), ref = make.names(ref_treatment))
  cd$time        <- relevel(factor(make.names(paste0("t", cd$time)), levels = time_levels),
                            ref = make.names(paste0("t", ref_time)))
  cd$temperature <- relevel(factor(make.names(paste0("T", cd$temperature)), levels = temp_levels),
                            ref = make.names(paste0("T", ref_temperature)))
  
  # sanity checks
  if (anyDuplicated(rownames(cd))) stop("Duplicated sample names.")
  cell <- interaction(cd$treatment, cd$time, cd$temperature, drop = FALSE)
  n_empty <- sum(table(cell) == 0)
  if (n_empty > 0) {
    warning(n_empty, " treatment x time x temperature cell(s) have no samples; ",
            "the full factorial model will not be full rank.")
  }
  
  message(nrow(cd), " samples | ",
          nlevels(cd$treatment), " treatments | ",
          nlevels(cd$time), " times | ",
          nlevels(cd$temperature), " temperatures")
  message("References: treatment = ", levels(cd$treatment)[1],
          ", time = ", levels(cd$time)[1],
          ", temperature = ", levels(cd$temperature)[1])
  print(ftable(table(cd$treatment, cd$time, cd$temperature,
                     dnn = c("treatment", "time", "temperature"))))
  cd
}

load_counts <- function(file,
                        sep = NULL,
                        id_col = "gene_id",
                        drop_cols = c("gene_id", "gene_name"),
                        samples = NULL) {   # pass rownames(colData) to subset and order
  if (is.null(sep)) sep <- if (grepl("\\.csv$", file, ignore.case = TRUE)) "," else "\t"
  
  df <- read.table(file, sep = sep, header = TRUE, check.names = FALSE,
                   quote = "\"", comment.char = "", stringsAsFactors = FALSE)
  
  counts <- as.matrix(df[, setdiff(colnames(df), drop_cols), drop = FALSE])
  rownames(counts) <- df[[id_col]]
  
  if (!is.null(samples)) {
    missing <- setdiff(samples, colnames(counts))
    if (length(missing)) stop("Samples missing from counts: ", paste(head(missing, 10), collapse = ", "))
    counts <- counts[, samples, drop = FALSE]
  }
  
  counts <- round(counts)          # salmon counts are fractional (Because of uncertainty!)
  storage.mode(counts) <- "integer"
  counts
}



plot_gene_padj_heatmap <- function(dds, gene,
                                   temp_var  = "temperature",
                                   treat_var = "treatment",
                                   time_var  = "time",
                                   alpha = 0.05,
                                   cap   = 5,
                                   get_res = function(nm) results(dds, name = nm)) {
  
  stopifnot(gene %in% rownames(dds))
  
  cd   <- colData(dds)
  vars <- c(temp_var, treat_var, time_var)
  
  # levels per factor; the first level is the reference
  lv  <- lapply(setNames(vars, vars), function(v) levels(factor(cd[[v]])))
  ref <- vapply(lv, `[`, character(1), 1)
  
  # model-matrix columns and resultsNames share the same order
  mm <- colnames(model.matrix(design(dds), cd))
  rn <- resultsNames(dds)
  stopifnot(length(mm) == length(rn))
  
  lookup <- do.call(rbind, lapply(vars, function(v)
    data.frame(var = v, level = lv[[v]][-1], token = paste0(v, lv[[v]][-1]),
               stringsAsFactors = FALSE)))
  
  # one row per coefficient: which cell it sits in and its padj / log2FC for this gene
  stats <- do.call(rbind, lapply(which(mm != "(Intercept)"), function(j) {
    info <- lookup[match(strsplit(mm[j], ":", fixed = TRUE)[[1]], lookup$token), ]
    if (anyNA(info$var)) return(NULL)          # not a plain factor-level term
    cell <- ref
    cell[info$var] <- info$level
    
    res_row <- get_res(rn[j])[gene, ]
    
    data.frame(Temperature = cell[[temp_var]],
               Treatment   = cell[[treat_var]],
               Time        = cell[[time_var]],
               coef = rn[j],
               padj = res_row$padj,
               log2FC = res_row$log2FoldChange,
               stringsAsFactors = FALSE)
  }))
  if (is.null(stats)) stop("No factor-level coefficients found in resultsNames(dds).")
  
  # full grid so cells without a coefficient stay visible (blank)
  grid <- merge(expand.grid(Temperature = lv[[temp_var]],
                            Treatment   = lv[[treat_var]],
                            Time        = lv[[time_var]],
                            stringsAsFactors = FALSE),
                stats, all.x = TRUE)
  
  grid$sig      <- !is.na(grid$padj) & grid$padj < alpha
  grid$nlogp    <- -log10(pmax(grid$padj, 1e-300))
  
  # Diverging fill: positive log2FC -> positive nlogp (green), negative log2FC -> negative nlogp (orange), NA if not significant
  grid$fill_val <- ifelse(grid$sig, 
                          ifelse(grid$log2FC < 0, -pmin(grid$nlogp, cap), pmin(grid$nlogp, cap)), 
                          NA_real_)
  
  is_ref <- grid$Temperature == ref[[temp_var]] &
    grid$Treatment   == ref[[treat_var]] &
    grid$Time        == ref[[time_var]]
  grid$label <- ifelse(is_ref, "ref",
                       ifelse(is.na(grid$padj), "", sprintf("%.1f", grid$nlogp)))
  
  # axis order: numeric for temperature and time; treatments with the reference on top
  num_order <- function(x) x[order(as.numeric(gsub("[^0-9.]", "", x)))]
  clean     <- function(x) sub("^[A-Za-z](?=[-0-9.])", "", x, perl = TRUE)  # "t0.5" -> "0.5"
  t_lv <- num_order(lv[[temp_var]])
  m_lv <- num_order(lv[[time_var]])
  grid$Temperature <- factor(grid$Temperature, levels = t_lv, labels = clean(t_lv))
  grid$Time        <- factor(grid$Time,        levels = m_lv, labels = clean(m_lv))
  grid$Treatment   <- factor(grid$Treatment,   levels = rev(lv[[treat_var]]))
  
  ggplot(grid, aes(x = Time, y = Treatment, fill = fill_val)) +
    geom_tile(color = "grey50", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3.2) +
    facet_wrap(~ Temperature, labeller = label_both) +
    scale_fill_gradient2(
      low = "darkorange", mid = "white", high = "forestgreen",
      midpoint = 0,
      name = expression(-log[10](padj)),
      na.value = "white",
      limits = c(-cap, cap),
      oob = scales::squish,
      labels = function(x) abs(x) # show positive values on color bar axis
    ) +
    theme_minimal() +
    labs(
      title    = paste("Significance Map (-log10 padj) for:", gene),
      subtitle = paste0("Green = Up (LFC > 0) | Orange = Down (LFC < 0) | padj < ", alpha, " (capped at 1e-", cap, ")"),
      x = "Time", y = "Treatment"
    ) +
    theme(
      strip.text    = element_text(face = "bold", size = 11),
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8.5),
      plot.caption  = element_text(hjust = 0.5, color = "grey40", size = 8),
      panel.grid    = element_blank(),
      panel.border  = element_rect(color = "grey75", fill = NA, linewidth = 0.8),
      panel.spacing = unit(1.5, "lines")
    )
}


compute_contrasts <- function(dds, 
                                   ref_treat = "mock",
                                   temp_var  = "temperature", 
                                   treat_var = "treatment", 
                                   time_var  = "time") {
  
  cd <- colData(dds)
  
  temps  <- levels(factor(cd[[temp_var]]))
  treats <- levels(factor(cd[[treat_var]]))
  times  <- levels(factor(cd[[time_var]]))
  
  ref_idx <- which(tolower(treats) == tolower(ref_treat))
  if (length(ref_idx) == 0) stop(sprintf("Reference treatment '%s' not found in %s.", ref_treat, treat_var))
  ref_label <- treats[ref_idx[1]]
  
  design_formula <- design(dds)
  mod_mat <- model.matrix(design_formula, cd)
  
  grid <- expand.grid(
    Temperature = temps, 
    Treatment   = treats, 
    Time        = times, 
    stringsAsFactors = FALSE
  )
  
  contrast_results <- list()
  
  for (i in 1:nrow(grid)) {
    t_val  <- grid$Temperature[i]
    tr_val <- grid$Treatment[i]
    ti_val <- grid$Time[i]
    
    key <- paste(t_val, tr_val, ti_val, sep = "__")
    
    if (tolower(tr_val) == tolower(ref_label)) {
      contrast_results[[key]] <- NULL
      next
    }
    
    test_selector <- (cd[[temp_var]] == t_val) & (cd[[treat_var]] == tr_val) & (cd[[time_var]] == ti_val)
    ref_selector  <- (cd[[temp_var]] == t_val) & (tolower(cd[[treat_var]]) == tolower(ref_label)) & (cd[[time_var]] == ti_val)
    
    if (!any(test_selector) || !any(ref_selector)) {
      contrast_results[[key]] <- NULL
      next
    }
    
    mean_test <- colMeans(mod_mat[test_selector, , drop = FALSE])
    mean_ref  <- colMeans(mod_mat[ref_selector, , drop = FALSE])
    contrast_vec <- mean_test - mean_ref
    
    cat(sprintf("=== Computing Contrast: %s vs %s | Temp: %s | Time: %s ===\n", tr_val, ref_label, t_val, ti_val))
    nz_terms <- contrast_vec[abs(contrast_vec) > 1e-8]
    print(round(nz_terms, 4))
    cat(sprintf("Total active coefficients involved: %d\n\n", length(nz_terms)))
    
    res <- tryCatch({
      results(dds, contrast = contrast_vec)
    }, error = function(e) {
      NULL
    })
    
    contrast_results[[key]] <- res
  }
  
  list(
    grid         = grid,
    results_list = contrast_results,
    ref_label    = ref_label,
    temp_var     = temp_var,
    treat_var    = treat_var,
    time_var     = time_var,
    temps        = temps,
    treats       = treats,
    times        = times
  )
}


