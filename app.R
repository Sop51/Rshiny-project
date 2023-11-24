## Author: Sophia Marcotte
## marcotts@bu.edu
## BU BF591
## Final Project

# import the required packages
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) 
library(dplyr)
library(tidyr)
library(DESeq2)

# increase the shiny upload size option for the counts data
options(shiny.maxRequestSize=30*1024^2) 

# UI FOR THE APPLICATION DEFINED HERE #
ui <- fluidPage(
  tabsetPanel(
    sidebarLayout(
      sidebarPanel(
        fileInput(inputId = "uploadedSampleFile", label = "Upload a file containing sample information:"),
        fileInput(inputId = "uploadedCountsFile", label = "Upload a file containing counts information:")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Samples",
                   tabsetPanel(
                     tabPanel("Summary", 
                              tableOutput("sampleTable")),
                     tabPanel("Table", 
                              radioButtons("filterChoice", "Choose Sample Data to View",
                                           choices = c("Control", "Huntington's")),
                              tableOutput("sampleFilteredTable")),
                     tabPanel("Plots", 
                              radioButtons("samplePlotChoice", "Choose Sample Information to View",
                                           choices = c("PMI","age_of_death", "RIN", "age_of_onset",
                                                       "duration", "CAG", "vonsattel_grade", "HV_striatal_score", "HV_cortical_score")),
                              plotOutput("sampleHistogram"))
                   )
          ),
          tabPanel("Counts",
                   tabsetPanel(
                     tabPanel("Filter Information",
                              sliderInput("varianceSlider", "Variance Percentile Threshold", min = 0, max = 100, value = 50),
                              sliderInput("nonzeroSlider", "Min Non-Zero Samples", min = 0, max = 5000, value = 10),
                              tableOutput("countsFilteredTable")),
                     tabPanel("Diagnostic Plot"),
                     tabPanel("Heatmap"),
                     tabPanel("PCA")
                   )),
          tabPanel("DE", "content")
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
  
  ### Define a function that reads in the counts data + normalizes ###
  load_count_data <- reactive({
    req(input$uploadedSampleFile)
    # read in the sample data
    df <- read.table(input$uploadedCountsFile$datapath, header = TRUE, sep = "\t")
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
    selected_df <- df %>% select(one_of(samples))
    return(selected_df)
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
    # return the filtered data
    return(data[selected_genes, , drop = FALSE])
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
  
}
# Run the application
shinyApp(ui = ui, server = server)
