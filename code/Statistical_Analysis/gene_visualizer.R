library(shiny)
library(ggplot2)
library(DT)
library(bslib)

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


# -----------------------------------------------------------------------------
# Trajectory plot
# -----------------------------------------------------------------------------
plot_gene_trajectory <- function(dds, gene,
                                 log2_scale = FALSE,
                                 shared_y   = FALSE,
                                 show_reps  = TRUE,
                                 nrow       = 2) {
  design_mat   <- model.matrix(design(dds), colData(dds))
  sample_names <- rownames(design_mat)
  split_names  <- strsplit(sample_names, "_")
  
  temps  <- unique(sapply(split_names, `[`, 1))
  treats <- unique(sapply(split_names, `[`, 2))
  times  <- unique(sapply(split_names, `[`, 3))
  
  temps <- temps[order(as.numeric(gsub("[^0-9.]", "", temps)))]
  times <- times[order(as.numeric(gsub("[^0-9.]", "", times)))]
  
  norm_counts <- counts(dds, normalized = TRUE)
  if (!(gene %in% rownames(norm_counts))) stop("Gene not found in counts table.")
  
  actual_data <- data.frame(
    Sample      = sample_names,
    Expression  = norm_counts[gene, ],
    Temperature = sapply(split_names, `[`, 1),
    Treatment   = sapply(split_names, `[`, 2),
    Time        = sapply(split_names, `[`, 3),
    stringsAsFactors = FALSE
  )
  
  grid <- expand.grid(Temperature = temps, Treatment = treats, Time = times,
                      stringsAsFactors = FALSE)
  
  grid$predicted_mu <- sapply(seq_len(nrow(grid)), function(i) {
    tryCatch(
      get_predicted_mean(dds, gene = gene,
                         temp  = grid$Temperature[i],
                         treat = grid$Treatment[i],
                         time  = grid$Time[i]),
      error = function(e) NA)
  })
  grid <- na.omit(grid)
  
  if (log2_scale) {
    actual_data$Expression <- log2(actual_data$Expression + 1)
    grid$predicted_mu      <- log2(grid$predicted_mu + 1)
    y_label  <- "Expression (log2 normalized counts + 1)"
    sub_text <- "Solid = predicted log2(mu) | Dashed = Mock | Faint points = replicates"
  } else {
    y_label  <- "Expression (normalized counts)"
    sub_text <- "Solid = predicted mu | Dashed = Mock | Faint points = replicates"
  }
  
  # Mock first, then the rest; mock grid is repeated in every facet
  mock_idx <- which(tolower(treats) == "mock")
  if (length(mock_idx) > 0) {
    mock_name    <- treats[mock_idx[1]]
    treat_levels <- c(mock_name, treats[-mock_idx[1]])
    grid_mock    <- subset(grid, Treatment == mock_name)
    grid_mock$Treatment <- NULL
  } else {
    treat_levels <- treats
    grid_mock    <- NULL
  }
  
  grid$Time        <- factor(grid$Time,        levels = times)
  grid$Temperature <- factor(grid$Temperature, levels = temps)
  grid$Treatment   <- factor(grid$Treatment,   levels = treat_levels)
  actual_data$Time        <- factor(actual_data$Time,        levels = times)
  actual_data$Temperature <- factor(actual_data$Temperature, levels = temps)
  actual_data$Treatment   <- factor(actual_data$Treatment,   levels = treat_levels)
  if (!is.null(grid_mock)) {
    grid_mock$Time        <- factor(grid_mock$Time,        levels = times)
    grid_mock$Temperature <- factor(grid_mock$Temperature, levels = temps)
  }
  
  temp_colors <- ifelse(temps == "28", "red", "dodgerblue3")
  names(temp_colors) <- temps
  
  p <- ggplot()
  
  if (!is.null(grid_mock) && nrow(grid_mock) > 0) {
    p <- p + geom_line(data = grid_mock,
                       aes(x = Time, y = predicted_mu, color = Temperature, group = Temperature),
                       linetype = "dashed", linewidth = 0.9, alpha = 0.45)
  }
  
  if (show_reps) {
    p <- p + geom_point(data = actual_data,
                        aes(x = Time, y = Expression, color = Temperature),
                        alpha = 0.4, size = 1.4,
                        position = position_jitter(width = 0.1, height = 0))
  }
  
  p +
    geom_line(data = grid,
              aes(x = Time, y = predicted_mu, color = Temperature, group = Temperature),
              linewidth = 1) +
    geom_point(data = grid,
               aes(x = Time, y = predicted_mu, color = Temperature),
               size = 2.5) +
    scale_color_manual(values = temp_colors, labels = function(x) paste0(x, " \u00b0C")) +
    facet_wrap(~ Treatment, nrow = nrow, scales = if (shared_y) "fixed" else "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = paste("Trajectory:", gene), subtitle = sub_text,
         x = "Time", y = y_label, color = NULL) +
    theme(
      strip.text       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 9),
      panel.spacing    = unit(1, "lines"),
      panel.border     = element_rect(color = "grey75", fill = NA, linewidth = 0.8),
      legend.position  = "bottom"
    )
}


# -----------------------------------------------------------------------------
# Compact contrast heatmap
# -----------------------------------------------------------------------------
plot_gene_diff <- function(computed_obj, gene,
                           metric = c("significance", "log2fc"),
                           alpha = 0.05, cap_mlog10pval = 10, cap_log2fc = 5,
                           hide_ref = TRUE,          # drop the (empty) Mock row
                           label_sig_only = TRUE) {  # print numbers only where padj < alpha
  
  metric <- match.arg(metric)
  
  grid         <- computed_obj$grid
  results_list <- computed_obj$results_list
  ref_label    <- computed_obj$ref_label
  temps        <- computed_obj$temps
  treats       <- computed_obj$treats
  times        <- computed_obj$times
  
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    t_val  <- grid$Temperature[i]
    tr_val <- grid$Treatment[i]
    ti_val <- grid$Time[i]
    key    <- paste(t_val, tr_val, ti_val, sep = "__")
    
    if (tolower(tr_val) == tolower(ref_label)) {
      return(data.frame(Temperature = t_val, Treatment = tr_val, Time = ti_val,
                        padj = 1.0, log2FC = 0.0, stringsAsFactors = FALSE))
    }
    res <- results_list[[key]]
    if (is.null(res) || !(gene %in% rownames(res))) {
      return(data.frame(Temperature = t_val, Treatment = tr_val, Time = ti_val,
                        padj = NA_real_, log2FC = NA_real_, stringsAsFactors = FALSE))
    }
    res_row <- res[gene, c("padj", "log2FoldChange")]
    data.frame(Temperature = t_val, Treatment = tr_val, Time = ti_val,
               padj = res_row$padj, log2FC = res_row$log2FoldChange,
               stringsAsFactors = FALSE)
  })
  gs <- do.call(rbind, rows)
  
  is_ref <- tolower(gs$Treatment) == tolower(ref_label)
  gs$sig   <- !is.na(gs$padj) & gs$padj < alpha & !is_ref
  gs$nlogp <- -log10(pmax(gs$padj, 1e-300))
  
  if (metric == "significance") {
    gs$fill_val <- ifelse(gs$sig,
                          ifelse(gs$log2FC < 0, -pmin(gs$nlogp, cap_mlog10pval),
                                 pmin(gs$nlogp, cap_mlog10pval)),
                          NA_real_)
    gs$label     <- sprintf("%.1f", gs$nlogp)
    fill_name    <- expression(-log[10](padj))
    title_suffix <- paste0("-log10 padj (padj < ", alpha, ")")
    labels_func  <- function(x) abs(x)
    limits_val   <- c(-cap_mlog10pval, cap_mlog10pval)
  } else {
    gs$fill_val <- ifelse(is_ref, NA_real_,
                          pmax(pmin(gs$log2FC, cap_log2fc), -cap_log2fc))
    gs$label     <- sprintf("%.2f", gs$log2FC)
    fill_name    <- "log2FC"
    title_suffix <- "log2 fold change"
    labels_func  <- waiver()
    limits_val   <- c(-cap_log2fc, cap_log2fc)
  }
  
  gs$label[is.na(gs$padj)] <- ""
  if (label_sig_only) gs$label[!gs$sig] <- ""
  gs$label[is_ref] <- if (hide_ref) "" else ref_label
  
  if (hide_ref) {
    gs     <- gs[!is_ref, , drop = FALSE]
    treats <- treats[tolower(treats) != tolower(ref_label)]
  }
  
  num_order <- function(x) x[order(as.numeric(gsub("[^0-9.]", "", x)))]
  clean     <- function(x) sub("^[A-Za-z](?=[-0-9.])", "", x, perl = TRUE)
  t_lv <- num_order(temps)
  m_lv <- num_order(times)
  
  gs$Temperature <- factor(gs$Temperature, levels = t_lv, labels = clean(t_lv))
  gs$Time        <- factor(gs$Time,        levels = m_lv, labels = clean(m_lv))
  gs$Treatment   <- factor(gs$Treatment,   levels = rev(treats))
  
  ggplot(gs, aes(x = Time, y = Treatment, fill = fill_val)) +
    geom_tile(color = "grey80", linewidth = 0.4) +
    geom_text(aes(label = label), size = 3) +
    facet_wrap(~ Temperature, nrow = 1,
               labeller = labeller(Temperature = function(x) paste0(x, " \u00b0C"))) +
    scale_fill_gradient2(
      low = "darkorange", mid = "white", high = "forestgreen", midpoint = 0,
      name = fill_name, na.value = "white",
      limits = limits_val, oob = scales::squish, labels = labels_func
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title    = paste0("Response vs. ", ref_label, " | ", gene),
      subtitle = paste0(title_suffix, " | green = up, orange = down | blank = n.s."),
      x = "Time", y = NULL
    ) +
    theme(
      strip.text       = element_text(face = "bold", size = 11),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 9),
      panel.grid       = element_blank(),
      panel.border     = element_rect(color = "grey75", fill = NA, linewidth = 0.8),
      panel.spacing    = unit(0.8, "lines"),
      legend.position  = "bottom",
      legend.key.width = unit(1.6, "cm"),
      legend.key.height = unit(0.3, "cm")
    )
}


# -----------------------------------------------------------------------------
# Full per-gene contrast table (feeds the "Table" tab)
# -----------------------------------------------------------------------------
gene_contrast_table <- function(computed_obj, gene) {
  g  <- computed_obj$grid
  rl <- computed_obj$results_list
  clean <- function(x) sub("^[A-Za-z](?=[-0-9.])", "", x, perl = TRUE)
  
  rows <- lapply(seq_len(nrow(g)), function(i) {
    key <- paste(g$Temperature[i], g$Treatment[i], g$Time[i], sep = "__")
    res <- rl[[key]]
    if (is.null(res) || !(gene %in% rownames(res))) return(NULL)
    data.frame(
      Temperature = clean(g$Temperature[i]),
      Treatment   = g$Treatment[i],
      Time        = clean(g$Time[i]),
      baseMean    = res[gene, "baseMean"],
      log2FC      = res[gene, "log2FoldChange"],
      lfcSE       = res[gene, "lfcSE"],
      pvalue      = res[gene, "pvalue"],
      padj        = res[gene, "padj"],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out[order(out$padj, na.last = TRUE), ]
}


# -----------------------------------------------------------------------------
# App
# -----------------------------------------------------------------------------
run_gene_viewer <- function(dds, computed_obj, default_gene = "AT5G24110") {
  
  all_genes <- rownames(dds)
  start_gene <- if (default_gene %in% all_genes) default_gene else all_genes[1]
  
  ui <- bslib::page_sidebar(
    title = "DESeq2 Gene Expression & Significance Explorer",
    theme = bslib::bs_theme(version = 5),
    
    sidebar = bslib::sidebar(
      width = 270,
      selectizeInput("gene", "Gene ID", choices = NULL,
                     options = list(placeholder = "Type a gene ID...", maxOptions = 50)),
      hr(),
      h6("Trajectory"),
      checkboxInput("log2_toggle", "Log2 scale", value = TRUE),
      checkboxInput("shared_y",    "Shared y-axis across panels", value = FALSE),
      checkboxInput("show_reps",   "Show replicate points", value = TRUE),
      hr(),
      h6("Contrast heatmap"),
      selectInput("heatmap_metric", NULL,
                  choices = c("Significance (-log10 padj)" = "significance",
                              "Log2 fold change"           = "log2fc")),
      checkboxInput("hide_ref", "Hide reference row (Mock)", value = TRUE),
      checkboxInput("sig_only", "Label significant cells only", value = TRUE)
    ),
    
    bslib::navset_card_tab(
      full_screen = TRUE,
      
      bslib::nav_panel(
        "Plots",
        bslib::layout_columns(
          col_widths = c(8, 4),
          fill = TRUE,
          bslib::card(
            full_screen = TRUE,
            bslib::card_header(
              class = "d-flex justify-content-between align-items-center",
              "Expression trajectory",
              downloadButton("dl_traj", "PNG", class = "btn-sm btn-outline-secondary")
            ),
            bslib::card_body(padding = 4, plotOutput("trajectoryPlot", height = "100%"))
          ),
          bslib::card(
            full_screen = TRUE,
            bslib::card_header(
              class = "d-flex justify-content-between align-items-center",
              "Treatment vs. Mock",
              downloadButton("dl_heat", "PNG", class = "btn-sm btn-outline-secondary")
            ),
            bslib::card_body(padding = 4, plotOutput("heatmapPlot", height = "100%"))
          )
        )
      ),
      
      bslib::nav_panel(
        "Table",
        DT::DTOutput("contrast_table")
      )
    )
  )
  
  server <- function(input, output, session) {
    
    # server-side selectize: fast even with tens of thousands of genes
    updateSelectizeInput(session, "gene", choices = all_genes,
                         selected = start_gene, server = TRUE)
    
    plot_cache <- new.env(hash = TRUE)
    cached <- function(key, build) {
      if (exists(key, envir = plot_cache, inherits = FALSE)) {
        return(get(key, envir = plot_cache, inherits = FALSE))
      }
      p <- build()
      assign(key, p, envir = plot_cache)
      p
    }
    
    gene <- reactive({
      g <- trimws(input$gene)
      req(nzchar(g), g %in% all_genes)
      g
    })
    
    traj_plot <- reactive({
      g <- gene()
      key <- paste("traj", g, isTRUE(input$log2_toggle), isTRUE(input$shared_y),
                   isTRUE(input$show_reps), sep = "|")
      cached(key, function()
        plot_gene_trajectory(dds, g,
                             log2_scale = isTRUE(input$log2_toggle),
                             shared_y   = isTRUE(input$shared_y),
                             show_reps  = isTRUE(input$show_reps),
                             nrow       = 2))
    })
    
    heat_plot <- reactive({
      g <- gene()
      key <- paste("heat", g, input$heatmap_metric, isTRUE(input$hide_ref),
                   isTRUE(input$sig_only), sep = "|")
      cached(key, function()
        plot_gene_diff(computed_obj, g,
                       metric         = input$heatmap_metric,
                       hide_ref       = isTRUE(input$hide_ref),
                       label_sig_only = isTRUE(input$sig_only)))
    })
    
    output$trajectoryPlot <- renderPlot(traj_plot(), res = 100)
    output$heatmapPlot    <- renderPlot(heat_plot(), res = 100)
    
    output$contrast_table <- DT::renderDT({
      tab <- gene_contrast_table(computed_obj, gene())
      DT::datatable(tab, rownames = FALSE, filter = "top",
                    extensions = "Buttons",
                    options = list(pageLength = 25, dom = "Bfrtip",
                                   buttons = c("copy", "csv"))) |>
        DT::formatRound(c("baseMean", "log2FC", "lfcSE"), 2) |>
        DT::formatSignif(c("pvalue", "padj"), 3)
    })
    
    output$dl_traj <- downloadHandler(
      filename = function() paste0(gene(), "_trajectory.png"),
      content  = function(file) ggsave(file, traj_plot(), width = 12, height = 7, dpi = 200)
    )
    output$dl_heat <- downloadHandler(
      filename = function() paste0(gene(), "_contrasts.png"),
      content  = function(file) ggsave(file, heat_plot(), width = 7, height = 6, dpi = 200)
    )
  }
  
  shinyApp(ui = ui, server = server)
}