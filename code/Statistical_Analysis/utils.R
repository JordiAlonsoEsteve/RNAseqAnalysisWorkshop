library(shiny)
library(ggplot2)
library(DT)
library(ggtext)

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


get_predicted_mean <- function(dds, gene, temp, treat, time) {
  # 1. Extract the full design matrix
  design_mat <- model.matrix(design(dds), colData(dds))
  sample_names <- rownames(design_mat)
  
  # 2. Find the design vector for the specified condition
  split_names <- strsplit(sample_names, "_")
  matched_idx <- which(sapply(split_names, function(x) {
    match_temp  <- is.null(temp)  || x[1] == temp
    match_treat <- is.null(treat) || x[2] == treat
    match_time  <- is.null(time)  || x[3] == time
    return(match_temp && match_treat && match_time)
  }))
  
  if (length(matched_idx) == 0) {
    stop("No samples found matching the specified criteria.")
  }
  
  # Design vector for this group (using the first matching replicate)
  x_vec <- design_mat[matched_idx[1], ]
  
  # 3. Extract the fitted betas for the specified gene
  beta_mat <- coef(dds)
  if (!(gene %in% rownames(beta_mat))) {
    stop("Gene not found in the coef(dds) table.")
  }
  beta_vec <- beta_mat[gene, ]
  
  # 4. Compute predicted log2 value and transform to natural scale (mu)
  predicted_log2 <- sum(x_vec * beta_vec)
  predicted_mu <- 2^predicted_log2
  
  # 5. Return only the predicted mu value
  return(predicted_mu)
}

plot_gene_trajectory <- function(dds, gene, log2_scale = FALSE) {
  design_mat <- model.matrix(design(dds), colData(dds))
  sample_names <- rownames(design_mat)
  split_names <- strsplit(sample_names, "_")
  
  temps <- unique(sapply(split_names, `[`, 1))
  treats <- unique(sapply(split_names, `[`, 2))
  times <- unique(sapply(split_names, `[`, 3))
  
  temp_nums <- as.numeric(gsub("[^0-9.]", "", temps))
  temps <- temps[order(temp_nums)]
  time_nums <- as.numeric(gsub("[^0-9.]", "", times))
  times <- times[order(time_nums)]
  
  norm_counts <- counts(dds, normalized = TRUE)
  if (!(gene %in% rownames(norm_counts))) {
    stop("Gene not found in counts table.")
  }
  
  actual_data <- data.frame(
    Sample = sample_names,
    Expression = norm_counts[gene, ],
    Temperature = sapply(split_names, `[`, 1),
    Treatment = sapply(split_names, `[`, 2),
    Time = sapply(split_names, `[`, 3),
    stringsAsFactors = FALSE
  )
  
  grid <- expand.grid(
    Temperature = temps, 
    Treatment = treats, 
    Time = times, 
    stringsAsFactors = FALSE
  )
  
  grid$predicted_mu <- sapply(1:nrow(grid), function(i) {
    tryCatch({
      get_predicted_mean(dds, gene = gene, 
                         temp = grid$Temperature[i], 
                         treat = grid$Treatment[i], 
                         time = grid$Time[i])
    }, error = function(e) { NA })
  })
  
  grid <- na.omit(grid)
  
  # Apply log2 scaling if requested
  if (log2_scale) {
    actual_data$Expression <- log2(actual_data$Expression + 1)
    grid$predicted_mu <- log2(grid$predicted_mu + 1)
    y_label  <- "Expression (log2 normalized counts + 1)"
    sub_text <- "Lines = Predicted log₂(μ), Dashed = Mock Control μ, Faint points = Actual Replicates"
  } else {
    y_label  <- "Expression (Normalized Counts)"
    sub_text <- "Lines = Predicted μ, Dashed = Mock Control μ, Faint points = Actual Replicates"
  }
  
  # Case-insensitive factor ordering & mock background reference data
  treats_lower <- tolower(treats)
  mock_idx <- which(treats_lower == "mock")
  
  if (length(mock_idx) > 0) {
    mock_name <- treats[mock_idx[1]]
    treat_levels <- c(mock_name, treats[-mock_idx])
    
    # Extract mock grid and drop 'Treatment' so it repeats across all facet panels
    grid_mock <- subset(grid, Treatment == mock_name)
    grid_mock$Treatment <- NULL
  } else {
    treat_levels <- treats
    grid_mock <- NULL
  }
  
  grid$Time <- factor(grid$Time, levels = times)
  grid$Temperature <- factor(grid$Temperature, levels = temps)
  grid$Treatment <- factor(grid$Treatment, levels = treat_levels)
  
  actual_data$Time <- factor(actual_data$Time, levels = times)
  actual_data$Temperature <- factor(actual_data$Temperature, levels = temps)
  actual_data$Treatment <- factor(actual_data$Treatment, levels = treat_levels)
  
  temp_colors <- ifelse(temps == "28", "red", "dodgerblue3")
  names(temp_colors) <- temps
  
  # Build plot with background mock dashed reference lines
  p <- ggplot()
  
  if (!is.null(grid_mock) && nrow(grid_mock) > 0) {
    p <- p + geom_line(data = grid_mock, 
                       aes(x = Time, y = predicted_mu, color = Temperature, group = Temperature), 
                       linetype = "dashed", linewidth = 0.9, alpha = 0.45)
  }
  
  p <- p +
    geom_point(data = actual_data, 
               aes(x = Time, y = Expression, color = Temperature), 
               alpha = 0.5, size = 1.8, position = position_jitter(width = 0.1, height = 0)) +
    geom_line(data = grid, 
              aes(x = Time, y = predicted_mu, color = Temperature, group = Temperature), 
              linewidth = 1) +
    geom_point(data = grid, 
               aes(x = Time, y = predicted_mu, color = Temperature), 
               size = 2.5) +
    scale_color_manual(values = temp_colors) +
    facet_wrap(~ Treatment, scales = "free_y") +
    theme_minimal() +
    labs(
      title = paste("Trajectory:", gene),
      subtitle = sub_text,
      x = "Time",
      y = y_label,
      color = "Temp"
    ) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(color = "grey75", fill = NA, linewidth = 0.8),
      legend.position = "bottom"
    )
  
  return(p)
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


run_gene_viewer <- function(dds) {
  
  ui <- fluidPage(
    titlePanel("DESeq2 Gene Expression & Significance Explorer"),
    uiOutput("page_content")
  )
  
  server <- function(input, output, session) {
    # App navigation state: "search" or "results"
    app_state <- reactiveVal("search")
    current_gene <- reactiveVal("")
    
    # In-memory cache for fast switching
    plot_cache <- new.env(hash = TRUE)
    
    # Handle search submission
    observeEvent(input$search_btn, {
      req(input$gene_input)
      g <- trimws(input$gene_input)
      
      if (!(g %in% rownames(dds))) {
        showNotification(paste("Gene ID", g, "not found in dds dataset."), type = "error")
        return()
      }
      
      current_gene(g)
      app_state("results")
    })
    
    # Handle back to search button
    observeEvent(input$back_btn, {
      app_state("search")
    })
    
    # Dynamic page layout switcher
    output$page_content <- renderUI({
      if (app_state() == "search") {
        # Landing / Search Page
        fluidPage(
          style = "max-width: 450px; margin-top: 120px; text-align: center;",
          h3("Search Gene"),
          p("Enter a valid gene ID from your dds dataset to inspect its trajectory and significance heatmap."),
          br(),
          textInput("gene_input", NULL, value = "AT5G24110", width = "100%"),
          br(),
          actionButton("search_btn", "Explore Gene", class = "btn-primary btn-lg", width = "100%")
        )
      } else {
        # Results Dashboard Page (Both plots side-by-side with Log2 control in header)
        fluidPage(
          fluidRow(
            column(12,
                   div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px; border-bottom: 1px solid #ddd; padding-bottom: 10px;",
                       h4(paste("Viewing Gene:", current_gene())),
                       div(style = "display: flex; align-items: center; gap: 20px;",
                           checkboxInput("log2_toggle", "Log2 Scale", value = TRUE), # Enabled by default
                           actionButton("back_btn", "← Search Another Gene", class = "btn-default btn-sm")
                       )
                   )
            )
          ),
          fluidRow(
            column(6,
                   tags$div(style = "background: #fff; padding: 10px; border-radius: 5px; border: 1px solid #ddd;",
                            plotOutput("trajectoryPlot", height = "580px")
                   )
            ),
            column(6,
                   tags$div(style = "background: #fff; padding: 10px; border-radius: 5px; border: 1px solid #ddd;",
                            plotOutput("heatmapPlot", height = "580px")
                   )
            )
          )
        )
      }
    })
    
    # Cached Trajectory Plot Generation (depends on log2 toggle state)
    traj_reactive <- reactive({
      req(app_state() == "results", current_gene())
      g <- current_gene()
      use_log2 <- isTRUE(input$log2_toggle)
      
      key <- paste0(g, "_traj_log2_", use_log2)
      if (exists(key, envir = plot_cache)) return(get(key, envir = plot_cache))
      
      p <- plot_gene_trajectory(dds = dds, gene = g, log2_scale = use_log2)
      assign(key, p, envir = plot_cache)
      return(p)
    })
    
    # Cached Heatmap Plot Generation
    heat_reactive <- reactive({
      req(app_state() == "results", current_gene())
      g <- current_gene()
      key <- paste0(g, "_heat")
      if (exists(key, envir = plot_cache)) return(get(key, envir = plot_cache))
      
      p <- plot_gene_padj_heatmap(dds = dds, gene = g)
      assign(key, p, envir = plot_cache)
      return(p)
    })
    
    output$trajectoryPlot <- renderPlot({
      p <- traj_reactive()
      if (!is.null(p)) p
    }, res = 100)
    
    output$heatmapPlot <- renderPlot({
      p <- heat_reactive()
      if (!is.null(p)) p
    }, res = 100)
  }
  
  shinyApp(ui = ui, server = server)
}
