# Interactive visualization of Arabidopsis RNASeq data
# Note: before running, working directory must be set to parent directory of app.R

# Load necessary libraries, functions, and data ----
packageList <- c("shiny", "shinyjs", "shinyWidgets", "dplyr", "DT", "ggplot2")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages) > 0){
  install.packages(newPackages)
}

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(dplyr)
library(DT)
library(ggplot2)

setwd("C:/Users/Bryce/Documents/interactive-data/arabidopsis-RNASeq/")
source("helpers.R")

load("data/RNAseq.rda")


# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Arabidopsis RNAseq"),
  fluidRow(
    column(4, wellPanel(
      helpText("Interactive tool for quick and easy subsetting of Arabidopsis RNAseq data"),
      searchInput(
        inputId="gene_ID_search",
        label="Enter a gene id to search for",
        placeholder="ex: AT3G01345",
        btnSearch=icon("search"),
        width="100%"
      )
    )),
    column(8,
      plotOutput("genePlot")
    ),
    
    column(12,
      h4("Select genes to plot"),
      DTOutput("selectedGenes"),
      DTOutput("expressionData")
    )
  )
)

server <- function(input, output, session){
  
  output$expressionData = renderDT(RNAseq, server=TRUE, filter="top", selection="multiple", 
    options=list(
      sDom="<'top'>lrt<'bottom'>ip",
      columnDefs=list(list(searchable=FALSE, targets=c(3, 4, 5, 6))))
  )
  
  output$selectedGenes = renderDT(RNAseq[input$expressionData_rows_selected, ], filter="none", selection="none")
  
  output$genePlot <- renderPlot({
    if(length(input$expressionData_rows_selected)>0){
      subsetData <- RNAseq[input$expressionData_rows_selected, ]
      ggplot(data=subsetData, mapping=aes(x=gene_id, y=ratio_FPKM)) +
        geom_col()
    }else{
      return(NA)
    }
    
    
  })
  
}



shinyApp(ui=ui, server=server)
