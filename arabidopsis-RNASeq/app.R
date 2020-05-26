# Interactive visualization of Arabidopsis RNASeq data
# Note: before running, working directory must be set to parent directory of app.R

# Load necessary libraries, functions, and data ----
packageList <- c("shiny", "shinyjs", "shinyWidgets", "dplyr", "DT", "ggplot2")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
#if(length(newPackages) > 0){
#  install.packages(newPackages)
#}

library(shiny)
library(shinyjs)
library(dplyr)
library(DT)
library(ggplot2)

source("helpers.R")

load("data/RNAseq.rda")

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  fluidRow(
    div(style="margin-top:0px", 
      column(3, 
        h2("RNAseq data viewer"),
        wellPanel(
          p("Tool for interactive viewing and subsetting of Arabidopsis RNAseq data. 
            Click on genes in the lower table to add them to the plot.
            Filter and search columns by using the table header.")
        )
      )
    ),
    div(style="margin-top:15px",
      column(5,
        plotOutput("genePlot", height="300px")
      ),
      column(4,
        DTOutput("selectedGenes")
      )
    ),
    column(12,
      DTOutput("expressionData")
    )
  )
)

server <- function(input, output, session){
  
  output$expressionData <- renderDT(RNAseq, 
    server=TRUE, filter="top", selection="multiple", class="compact", 
    options=list(
      sDom="<'top'>lrt<'bottom'>ip",
      autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="10%", targets=c(4, 5, 6, 7))),
      pageLength=10,
      lengthMenu=c(10, 15, 20, 25))
  )
  
  output$selectedGenes <- renderDT(RNAseq[input$expressionData_rows_selected, c(1, 2, 3)],
    class="compact", selection="none",
    options=list(
      dom="rti",
      autoWidth=TRUE,
      scrollY="250px"
    )
  )
  
  output$genePlot <- renderPlot({
    if(length(input$expressionData_rows_selected)>0){
      subsetData <- RNAseq[input$expressionData_rows_selected, ]
      WT_data <- subsetData[ , 1:5]
      WT_data$genotype <- "WT"
      colnames(WT_data)[4:5] <- c("mean_FPKM", "se_FPKM")
      Ref5_data <- subsetData[ , c(1:3, 6:7)]
      Ref5_data$genotype <- "Ref5"
      colnames(Ref5_data)[4:5] <- c("mean_FPKM", "se_FPKM")
      subsetData <- rbind(WT_data, Ref5_data)
      
      ggplot(data=subsetData, mapping=aes(x=gene_id, y=mean_FPKM, fill=genotype)) +
        geom_col(position=position_dodge()) +
        geom_errorbar(mapping=aes(ymin=mean_FPKM-se_FPKM, ymax=mean_FPKM+se_FPKM), color="black",
                      width=0.2, position=position_dodge(0.9)) +
        theme(legend.position="right")
    }else{
      return(NA)
    }
  })
  
}

shinyApp(ui=ui, server=server)
