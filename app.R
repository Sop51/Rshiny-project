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
                     tabPanel("Plots", "Content for Subtab 1.3")
                   )
          ),
          tabPanel("Counts", "content"),
          tabPanel("DE", "content")
        )
      )
    )
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # Define a function that reads in the sample data
  load_sample_data <- reactive({
    req(input$uploadedSampleFile)
    
    # Read in the sample data
    df <- read.table(input$uploadedSampleFile$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    # Identify the columns to convert to numeric (exclude 'sample_ID')
    numeric_columns <- setdiff(names(df), "sample_ID")
    
    # Convert selected columns to numeric, handling non-numeric values
    df[numeric_columns] <- lapply(df[numeric_columns], function(x) {
      # Remove commas and convert to numeric
      as.numeric(gsub(",", "", as.character(x)))
    })
    
    return(df)
  })
  
  #create a function that returns the summary table
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
  
  # create a function that filters the data based on input value (H vs C)
  filter_samples <- 
    function(data, keyword) {
    keyword <- toupper(keyword)
    filtered_data <- data %>%
      filter(grepl(paste0("^", substr(keyword, 1, 1), "_"), sample_ID, ignore.case = TRUE))
    return(filtered_data)
  }
  
  output$sampleTable <- renderTable({
    dataf <- load_sample_data()
    summary_table(dataf)
  })
  
  output$sampleFilteredTable<- renderTable({
    dataf <- load_sample_data()
    filter_samples(dataf, input$filterChoice)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
