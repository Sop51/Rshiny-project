## Author: Sophia Marcotte
## marcotts@bu.edu
## BU BF591
## Final Project

# import the required packages
library(shiny)
library(bslib)
library(ggplot2)
library(gplots)
library(colourpicker) 
library(dplyr)
library(tidyr)
library(DESeq2)
library(pheatmap)
library(igraph)

# increase the shiny upload size option for the counts data
options(shiny.maxRequestSize=30*1024^2) 

# UI FOR THE APPLICATION DEFINED HERE #
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "uploadedSampleFile", label = "Upload a file containing sample information:"),
      fileInput(inputId = "uploadedCountsFile", label = "Upload a file containing counts information:")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Samples",
          tabsetPanel(
            tabPanel("Summary", tableOutput("sampleTable")),
            tabPanel("Table", radioButtons("filterChoice", "Choose Sample Data to View", choices = c("Control", "Huntington's")),
                     tableOutput("sampleFilteredTable")),
            tabPanel("Plots", radioButtons("samplePlotChoice", "Choose Sample Information to View",
                                           choices = c("PMI","age_of_death", "RIN", "age_of_onset",
                                                       "duration", "CAG", "vonsattel_grade", "HV_striatal_score", "HV_cortical_score")),
                     plotOutput("sampleHistogram"))
          )
        ),
        tabPanel(
          "Counts",
          sliderInput("varianceSlider", "Variance Percentile Threshold", min = 0, max = 100, value = 50),
          sliderInput("nonzeroSlider", "Min Non-Zero Samples", min = 0, max = 100, value = 10),
          tabsetPanel(
            tabPanel("Filter Information", tableOutput("countsFilteredTable")),
            tabPanel("Diagnostic Plot", plotOutput("varDiagPlot"), plotOutput("zeroDiagPlot")),
            tabPanel("Heatmap", plotOutput("heatmapDiagPlot")),
            tabPanel("PCA", radioButtons("pc_components", "Select Principal Components:",
                                         choices = c("PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3")),
                     plotOutput("pca_plot"))
          )
        ),
        tabPanel(
          "DE",
          tabsetPanel(
            tabPanel("Differential Expression Results",
                     sidebarPanel(sliderInput("threshold", "Choose a P-Value Threshold", min = 0, max = 0.5, value = 0.05, step = 0.001)),
                     mainPanel(tableOutput("detable"))
            ),
            tabPanel("Differential Expression Plots",
                     p("Adjusted P-Values From DESeq2 Results"),
                     plotOutput("plotPval"),
                     p("Log2FC Values of Genes Significant at Padj Threshold of 0.1"),
                     plotOutput("log2FCplot"),
                     p("Log2foldchange vs -log10(padj), Labeled by Status"),
                     plotOutput("volcanoPlot"))
          )
        ),
        tabPanel(
          "Correlation Network Analysis",
          textAreaInput("geneInput", "Enter Gene Names (one per line):", rows = 10),
          tabsetPanel(
            tabPanel("Clustered Heatmap",
                     actionButton("plotHeatmapButton", "Plot Heatmap"),
                     plotOutput("corrHeatmapPlot")
                     ),
            tabPanel("Visualization",
                     sliderInput("minCorr", "Min Correlation for Drawing an Edge", min = 0, max = 1, value = 0.05),
                     actionButton("plotCorrButton", "Plot Correlation Graph"),
                     plotOutput("networkPlot")
            ),
            tabPanel("Metrics",
                     tableOutput("networkStatsTable")
            )
          )
        )
      )
    )
  )
)


# SERVER LOGIC DEFINED HERE #
server <- function(input, output, session) {

  ### Define a function that reads in the sample data ###
  load_sample_data <- reactive({
    req(input$uploadedSampleFile)
    # read in the sample data
    df <- read.table(input$uploadedSampleFile$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # identify the columns to convert to numeric (exclude 'sample_ID')
    numeric_columns <- setdiff(names(df), "sample_ID")
    # convert selected columns to numeric, handling non-numeric values
    df[numeric_columns] <- lapply(df[numeric_columns], function(x) {
      # remove commas and convert to numeric
      as.numeric(gsub(",", "", as.character(x)))
    })
    return(df)
  })
  
  ### Define a function that reads in the counts data (BA9) ###
  load_count_data <- reactive({
    req(input$uploadedSampleFile)
    # read in the sample data
    df <- read.table(input$uploadedCountsFile$datapath, header = TRUE, sep = "\t")
    # rename the gene column
    df <- df %>%
      rename("X" = "gene")
    # define the samples to extract 
    samples <- c(
      "gene", "H_0001", "H_0002", "H_0003",
      "H_0005", "H_0006", "H_0007",
      "H_0008", "H_0009", "H_0010",
      "H_0012", "H_0013", "H_0539",
      "H_0657", "H_0658", "H_0681",
      "H_0695", "H_0700", "H_0726",
      "H_0740", "H_0750", "C_0012",
      "C_0013", "C_0014", "C_0015",
      "C_0016", "C_0017", "C_0018",
      "C_0020", "C_0021", "C_0022",
      "C_0023", "C_0024", "C_0025",
      "C_0026", "C_0029", "C_0031",
      "C_0032", "C_0033",
      "C_0035", "C_0036", "C_0037",
      "C_0038", "C_0039", "C_0050",
      "C_0053", "C_0060", "C_0061",
      "C_0062", "C_0065", "C_0069",
      "C_0070", "C_0071", "C_0075",
      "C_0076", "C_0077", "C_0081",
      "C_0082", "C_0083", "C_0087",
      "C_0002", "C_0003", "C_0004",
      "C_0005", "C_0006", "C_0008",
      "C_0009", "C_0010", "C_0011"
    )
    # selected the required samples
    selected_df <- df %>% select(one_of(samples))
    return(selected_df)
  })
  
  ### define a function that loads in the DE data ###
  load_de_data <- reactive({
    data <- read.table('/Users/sophiemarcotte/Desktop/deseq2_results.txt', header = TRUE, sep = "\t")
    return(data)
  })
  
  ### load corr data and filter ###
  corr_filtered_data <- reactive({
      data <- load_count_data()
      # check if the input is not NULL
      if (!is.null(input$geneInput)) {
        # split the input into a vector of gene names
        selected_genes <- strsplit(input$geneInput, "\n", fixed = TRUE)[[1]]
        # remove any leading or trailing whitespace
        selected_genes <- trimws(selected_genes)
        # filter data based on selected genes
        subset_data <- data[data$gene %in% selected_genes, ]
        return(subset_data)
      }else {
        return(NULL)
      }
  })
  
  ### Define a function that returns the summary table for each column ###
  summary_table <-
    function(data) {
      # extract group information from the sample ID (control versus Huntingtons)
      data <- data %>%
        mutate(Condition = ifelse(grepl("^H", sample_ID), "Huntington's", "Control"))
      # create a summary data frame for each column
      summary_data <- data %>%
        # group by column and condition
        gather(key = "Column", value = "Value", -sample_ID, -Condition) %>%
        group_by(Column, Condition) %>%
        summarise(
          DataType = class(Value),
          'Mean (sd)' = switch(
            class(Value),
            # for numeric column, compute sd and mean - format with parenthesis
            "numeric" = ifelse(all(is.na(Value)), NA,
                               paste0(format(mean(as.numeric(Value), na.rm = TRUE), digits = 2),
                                      " (", format(sd(as.numeric(Value), na.rm = TRUE), digits = 2), ")")),
            NA #if no other data type (none in this case)
          )
        )
      # return the data frame
      return(summary_data)
    }
  
  ### define a function that filters the sample data based on input value (H vs C) ###
  filter_samples <- 
    function(data, keyword) {
    keyword <- toupper(keyword)
    filtered_data <- data %>%
      filter(grepl(paste0("^", substr(keyword, 1, 1), "_"), sample_ID, ignore.case = TRUE))
    return(filtered_data)
    }
  
  ### define function to create histograms for each column (excluding the first one) ###
  create_histograms <- function(data) {
    # exclude the first column
    numeric_columns <- names(data)[-1]
    # list to store histograms
    histograms_list <- list()
    # create a histogram for each numeric column
    for (col in numeric_columns) {
      plot <- ggplot(data, aes(x = !!sym(col))) +
        geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
        labs(title = col) +
        theme_minimal()
      # add the plot to the list
      histograms_list[[col]] <- plot
    }
    return(histograms_list)
  }
  
  ### define function to filter genes based on variance and non-zero samples ###
  filter_genes <- function(data, var_percentile, min_nonzero_samples) {
    # remove rows with any missing values
    data <- na.omit(data)
    # calculate variance for each gene
    gene_var <- apply(data[-1], 1, var, na.rm = TRUE)
    # calculate the percentile threshold based on var_percentile
    var_threshold <- quantile(gene_var, var_percentile / 100, na.rm = TRUE)
    # calculate the number of non-zero samples for each gene
    non_zero_samples <- apply(data[-1] != 0, 1, sum, na.rm = TRUE)
    # filter genes based on criteria
    selected_genes <- which(gene_var >= var_threshold & non_zero_samples >= min_nonzero_samples)
    # create a summary table
    summary_table <- data.frame(
      "Number of Samples" = ncol(data) - 1,  # subtract 1 for the gene column
      "Total Number of Genes" = nrow(data),
      "Number of Genes Passing Filter" = length(selected_genes),
      "Percentage of Genes Passing Filter" = (length(selected_genes) / nrow(data)) * 100,
      "Number of Genes Not Passing Filter" = nrow(data) - length(selected_genes),
      "Percentage of Genes Not Passing Filter" = ((nrow(data) - length(selected_genes)) / nrow(data)) * 100
    )
    # return the summary table
    return(summary_table)
  }
  
  ### define a function that creates a diagnostic plot for zero filter ###
  create_zeros_plot <- function(filtered_data, min_nonzero_samples) {
    # calculate median count
    median_count <- apply(filtered_data[-1], 1, median, na.rm = TRUE)
    # calculate the number of non-zero samples for each gene
    non_zero_samples <- apply(filtered_data[-1] != 0, 1, sum, na.rm = TRUE)
    # filter genes based on the minimum number of non-zero samples
    selected_genes <- which(non_zero_samples >= min_nonzero_samples)
    # create a data frame for the plot
    plot_data <- data.frame(
      MedianCount = median_count[selected_genes],
      Zeros = apply(filtered_data[selected_genes, -1] == 0, 1, sum, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    # create the scatter plot
    plot_zeros <- ggplot(plot_data, aes(x = Zeros, y = log(MedianCount))) +
      geom_point() +
      labs(title = "Median Count vs Number of Zeros",
           x = "Number of Zeros",
           y = "Log(Median Count)")
    # return the plot
    return(plot_zeros)
  }
  
  ### define a function that creates a diagnostic plot for variance filter ###
  create_variance_plot <- function(filtered_data, min_variance_percentile) {
    # calculate median count
    median_count <- apply(filtered_data[-1], 1, median, na.rm = TRUE)
    # calculate the threshold based on the minimum percentile of variance
    variance_threshold <- quantile(apply(filtered_data[-1], 1, var, na.rm = TRUE), min_variance_percentile / 100, na.rm = TRUE)
    # filter genes based on the minimum percentile of variance
    selected_genes <- which(apply(filtered_data[-1], 1, var, na.rm = TRUE) >= variance_threshold)
    # create a data frame for the plot
    plot_data <- data.frame(
      MedianCount = median_count[selected_genes],
      Variance = apply(filtered_data[selected_genes, -1], 1, var, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    # create the scatter plot
    plot_variance <- ggplot(plot_data, aes(x = log(Variance), y = log(MedianCount))) +
      geom_point() +
      labs(title = "Median Count vs Variance",
           x = "Log(Variance)",
           y = "Log(Median Count)")
    # return the plot
    return(plot_variance)
  }
  
  ### define a function that produces a heatmap based on slider var and nonzero values ###
  create_clustered_heatmap <- function(filtered_data, min_variance_percentile, min_nonzero_samples) {
    # filter genes based on the minimum percentile of variance
    variance_threshold <- quantile(apply(filtered_data[-1], 1, var, na.rm = TRUE), min_variance_percentile / 100, na.rm = TRUE)
    selected_genes <- which(apply(filtered_data[-1], 1, var, na.rm = TRUE) >= variance_threshold)
    # filter genes based on the minimum number of non-zero samples
    non_zero_samples <- apply(filtered_data[-1] != 0, 1, sum, na.rm = TRUE)
    selected_genes <- intersect(selected_genes, which(non_zero_samples >= min_nonzero_samples))
    # extract the selected data
    selected_data <- filtered_data[selected_genes, ]
    # subset the numeric matrix (excluding the first column, assuming it contains sample IDs)
    numeric_matrix <- as.matrix(selected_data[, -1])
    # log-transform the counts
    log_transformed_data <- log1p(numeric_matrix)  # log1p to handle zeros
    # create a clustered heatmap
    pheatmap(
      log_transformed_data,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "row",
      col = colorRampPalette(c("white", "blue"))(50),
      main = "Clustered Heatmap",
      show_rownames = FALSE,  # remove row labels
      fontsize_row = 8,       # adjust font size of the row labels
      fontsize_col = 8        # adjust font size of the column labels
    )
  }
  
  ### define a function to produce the pca plot ###
  plot_pca <- function(data, pc_components) {
    # set genes as the row names
    rownames(data) <- data[, 1]
    data <- data[, -1]
    # transpose the data
    transposed_data <- t(data)
    # perform PCA
    pc <- prcomp(transposed_data)
    # define components based on user input
    if (pc_components == "PC1 vs PC2") {
      pc_components <- c(1, 2)
    } else if (pc_components == "PC1 vs PC3") {
      pc_components <- c(1, 3)
    } else if (pc_components == "PC2 vs PC3") {
      pc_components <- c(2, 3)
    } else {
      stop("Invalid PCA option")
    }
    # create PCA results data frame
    pca_results <- data.frame(pc$x[, pc_components])
    pca_results$sample_ID <- rownames(pc$x)
    # extract sample type from sample_ID
    pca_results$sample_type <- ifelse(grepl("^H_", pca_results$sample_ID), "Huntingtons", "Control")
    # calculate the percentage of variance explained by each selected component
    pc_var_percent <- pc$sdev[pc_components]^2 / sum(pc$sdev^2) * 100
    # create scatter plot
    pca_plot <- ggplot(pca_results, aes(x = !!sym(paste0("PC", pc_components[1])),
                                      y = !!sym(paste0("PC", pc_components[2])),
                                      color = sample_type)) + # color by treat vs control
      geom_point() +
      labs(
        x = paste("PC", pc_components[1], " (", round(pc_var_percent[1], 2), "% variance)"),
        y = paste("PC", pc_components[2], " (", round(pc_var_percent[2], 2), "% variance)")
      )
    
    return(pca_plot)
  }
  
  ### define a plot to filter the DE results based on a threshold
  filter_de_data <- function(data, threshold) {
    filter(data, padj <= threshold)
  }
  
  ### create a plot that plots the pvals of the DE results ###
  plot_pvals <- function(dataf) {
    p <- ggplot(dataf, aes(x=padj)) + 
      theme_classic() +
      geom_histogram()
    return(p)
  }
  
  ### create a plot that plots the log2foldchange of DE results ###
  plot_log2fc <- function(dataf) {
    dataf <- dataf %>%
      mutate(volc_plot_status = case_when(
        padj < 0.1 & log2FoldChange > 0 ~ "UP",
        padj < 0.1 & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      ))
    plot_sig <- dataf[dataf$volc_plot_status != 'NS',]
    p <- ggplot(plot_sig, aes(x=log2FoldChange)) + 
      theme_classic() +
      geom_histogram()
    return(p)
  }
  
  ### create a volcano plot of the DE results ###
  plot_volcano <- function(dataf) {
    dataf <- dataf %>%
      mutate(log10neg = -log10(padj))
    dataf <- dataf %>%
      mutate(Status = case_when(
        log2FoldChange > 0.5 ~ "Up",
        log2FoldChange < -0.5 ~ "Down",
        TRUE ~ "Not Significant"
      ))
    p <- ggplot(dataf, aes(x=log2FoldChange, y=log10neg, color=Status)) +
      geom_point()
    return(p)
  }
  
  ### define a function to create a heatmap for corr analysis ###
  plot_corrHeatmap <- function(dataf){
    rownames(dataf) <- dataf[, 1]
    # Remove the first column (as it's now the row names)
    dataf <- dataf[, -1]
    dataf <- as.matrix(dataf)
    # log-transform the counts
    log_transformed_data <- log1p(dataf)  # log1p to handle zeros
    # create a clustered heatmap
    pheatmap(
      log_transformed_data,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "row",
      col = colorRampPalette(c("white", "blue"))(50),
      main = "Clustered Heatmap",
      show_rownames = TRUE,
      fontsize_row = 8,       # adjust font size of the row labels
      fontsize_col = 8)        # adjust font size of the column labels
  }
  
  ### define a function to compute the correlation matrix and vizualize with a graph ###
  correlation_matrix_viz <- function(selected_genes){
    # set first col as rownames
    rownames(selected_genes) <- selected_genes[, 1]
    # Remove the first column (as it's now the row names)
    selected_genes <- selected_genes[, -1]
    # calculate the correlation matrix
    cor_matrix <- cor(t(selected_genes))
    # filter the correlation matrix based on the threshold
    filtered_correlation_matrix <- ifelse(abs(cor_matrix) >= input$minCorr, cor_matrix, NA)
    # create the igraph object
    graph <- graph_from_adjacency_matrix(filtered_correlation_matrix, mode = "undirected", weighted = TRUE)
    # plot the network
    plot(graph, 
         vertex.color = "lightblue",    # vertex color
         vertex.label.cex = 0.8,        # vertex label size
         vertex.label.color = "black",  # vertex label color
         edge.arrow.size = 0.5,         # edge arrow size
         layout = layout_with_fr)       # graph layout
  }
  
  ### define a function to calculate metrics calculated for each gene in the input ###
  corr_matrix_metrics <- function(selected_genes){
    # set first col as rownames
    rownames(selected_genes) <- selected_genes[, 1]
    # Remove the first column (as it's now the row names)
    selected_genes <- selected_genes[, -1]
    # calculate the correlation matrix
    cor_matrix <- cor(t(selected_genes))
    # create the igraph object
    graph <- graph_from_adjacency_matrix(cor_matrix, mode = "undirected", weighted = TRUE)
    # calculate degree, closeness centrality, and betweenness centrality
    degree_values <- degree(graph)
    closeness_values <- closeness(graph)
    betweenness_values <- betweenness(graph)
    # create a data frame with the results
    results_table <- data.frame(
      Gene = V(graph)$name,
      Degree = degree_values,
      Closeness = closeness_values,
      Betweenness = betweenness_values
    )
    return(results_table)
  }
  
  # THESE ARE THE RENDER OUTPUT FUNCTIONS #
  output$sampleTable <- renderTable({
    dataf <- load_sample_data()
    summary_table(dataf)
  })
  
  output$sampleFilteredTable<- renderTable({
    dataf <- load_sample_data()
    filter_samples(dataf, input$filterChoice)
  })
  
  output$sampleHistogram <- renderPlot({
    dataf <- load_sample_data()
    histograms <- create_histograms(dataf)
    if (input$samplePlotChoice == "PMI"){
      return(histograms[1])
    }
    else if (input$samplePlotChoice == "age_of_death"){
      return(histograms[2])
    }
    else if (input$samplePlotChoice == "RIN"){
      return(histograms[3])
    }
    else if (input$samplePlotChoice == "age_of_onset"){
      return(histograms[5])
    }
    else if (input$samplePlotChoice == "duration"){
      return(histograms[6])
    }
    else if (input$samplePlotChoice == "CAG"){
      return(histograms[7])
    }
    else if (input$samplePlotChoice == "vonsattel_grade"){
      return(histograms[8])
    }
    else if (input$samplePlotChoice == "HV_striatal_score"){
      return(histograms[9])
    }
    else{
      return(histograms[10])
    }
  })
  
  output$countsFilteredTable <- renderTable({
    dataf <- load_count_data()
    filter_genes(dataf, input$varianceSlider, input$nonzeroSlider)
  })
  
  output$varDiagPlot <- renderPlot({
    dataf <- load_count_data()
    create_variance_plot(dataf, input$varianceSlider)
  })
  
  output$zeroDiagPlot <- renderPlot({
    dataf <- load_count_data()
    create_zeros_plot(dataf, input$nonzeroSlider)
  })
  
  output$heatmapDiagPlot <- renderPlot({
    dataf <- load_count_data()
    create_clustered_heatmap(dataf, input$varianceSlider, input$nonzeroSlider)
  })
  
  output$pca_plot <- renderPlot({
    datac <- load_count_data()
    plot_pca(datac, input$pc_components)
  })
  
  output$detable <- renderTable({
    dataf <- load_de_data()
    filter_de_data(dataf, input$threshold)
  })
  
  output$plotPval <- renderPlot({
    dataf <- load_de_data()
    plot_pvals(dataf)
  })
  
  output$log2FCplot <- renderPlot({
    dataf <- load_de_data()
    plot_log2fc(dataf)
  })
  
  output$volcanoPlot <- renderPlot({
    dataf <- load_de_data()
    plot_volcano(dataf)
  })
  
  output$corrHeatmapPlot <- renderPlot({
    if (input$plotHeatmapButton > 0) {
      dataf <- corr_filtered_data()
      plot_corrHeatmap(dataf)
    }
  })
  
  output$networkPlot <- renderPlot({
    if (input$plotCorrButton > 0 ){
      dataf <- corr_filtered_data()
      correlation_matrix_viz(dataf)
    }
  })
  
  output$networkStatsTable <- renderTable({
    dataf <- corr_filtered_data()
    corr_matrix_metrics(dataf)
  })
}
# Run the application
shinyApp(ui = ui, server = server)
