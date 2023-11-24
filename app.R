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
                     tabPanel("Table", "Content for Subtab 1.2"),
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
  
  # define a function that reads in the sample data
  load_sample_data <- reactive({
    req(input$uploadedSampleFile)
    df <- read.table(input$uploadedSampleFile$datapath, header = TRUE, sep = "\t")
    # Identify the columns to convert to numeric (exclude 'Sample')
    numeric_columns <- setdiff(names(df), "sample_ID")
    # Convert selected columns to numeric
    df[numeric_columns] <- apply(df[numeric_columns], 2, function(x) as.numeric(as.character(x)))
    return(df)
  })
  
  #create a function that returns the summary table
  summary_table <-
    function(data) {
      # Extract group information from the sample ID
      data <- data %>%
        mutate(Condition = ifelse(grepl("^H", sample_ID), "H", "C"))
      
      # Create a summary data frame for each column
      summary_data <- data %>%
        gather(key = "Column", value = "Value", -sample_ID, -Condition) %>%
        group_by(Column, Condition) %>%
        summarise(
          DataType = class(Value),
          'Mean (sd)' = switch(
            class(Value),
            "numeric" = ifelse(all(is.na(Value)), NA,
                               paste0(format(mean(as.numeric(Value), na.rm = TRUE), digits = 2),
                                      " (", format(sd(as.numeric(Value), na.rm = TRUE), digits = 2), ")")),
            NA
          )
        )
      return(summary_data)
    }
  
  # create a function that returns a volcano plot
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
    }
  

  draw_table <- function(dataf, slider) {
  }
  
  output$volcano <- renderPlot({
  })
  
  output$sampleTable <- renderTable({
    dataf <- load_sample_data()
    summary_table(dataf)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
