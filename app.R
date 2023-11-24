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
library(gridExtra)
library(png)

# define the ui for the application
ui <- fluidPage(
  tabsetPanel(
    sidebarLayout(
      sidebarPanel(
        fileInput(inputId = "uploadedSampleFile", label = "Upload a CSV file containing sample information:")
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
          tabPanel("Counts", "content"),
          tabPanel("DE", "content")
        )
      )
    )
  )
)


# Define server logic required
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
  
  # define a function that filters the sample data based on input value (H vs C)
  filter_samples <- 
    function(data, keyword) {
    keyword <- toupper(keyword)
    filtered_data <- data %>%
      filter(grepl(paste0("^", substr(keyword, 1, 1), "_"), sample_ID, ignore.case = TRUE))
    return(filtered_data)
    }
  
  # Function to create histograms for each column (excluding the first one)
  create_histograms <- function(data) {
    # Exclude the first column
    numeric_columns <- names(data)[-1]
    
    # List to store histograms
    histograms_list <- list()
    
    # Create a histogram for each numeric column
    for (col in numeric_columns) {
      plot <- ggplot(data, aes(x = !!sym(col))) +
        geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
        labs(title = col) +
        theme_minimal()
      
      # Add the plot to the list
      histograms_list[[col]] <- plot
    }
    
    return(histograms_list)
  }
  
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
}
# Run the application
shinyApp(ui = ui, server = server)
