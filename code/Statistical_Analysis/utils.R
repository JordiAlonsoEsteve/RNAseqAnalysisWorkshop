build_coldata <- function(file,
                          sep = NULL,
                          id_cols = c("gene_id", "gene_name"),
                          drop_treatments = "untreated",
                          ref_treatment = "mock",
                          ref_temperature = 21,
                          ref_time = NULL) {   # NULL = smallest time left after dropping treatments
  # delimiter: comma for .csv, tab otherwise (unless given)
  if (is.null(sep)) sep <- if (grepl("\\.csv$", file, ignore.case = TRUE)) "," else "\t"
  
  # only the header is needed
  header <- read.table(file, sep = sep, header = TRUE, nrows = 1,
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