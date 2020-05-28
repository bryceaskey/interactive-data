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
        tags$h2("RNAseq data viewer"),
        wellPanel(
          tags$p("Tool for interactive viewing and subsetting of Arabidopsis RNAseq data."),
          tags$br(),
          tags$p(tags$strong("To use:")),
          tags$ol(
            tags$li("Search and filter for genes in the lower table with the header row of the
                    table."),
            tags$li("Select genes to plot from the lower table by clicking on them."),
            tags$li("Genes being plotted will appear in the table to the right of the plot 
                    in the order that they were selected."),
            tags$li("To remove genes from the plot, click on them in the right table")
               
          )
        )
      )
    ),
    tags$div(style="margin-top:60px; margin-bottom:470px",
      column(6,
        plotOutput("genePlot", height="390px")
      ),
      column(3,
        DTOutput("selectedGenes")
      )
    ),
    column(12,
      DTOutput("expressionData")
    )
  )
)

if (!exists("default_search_columns")) default_search_columns <- NULL

server <- function(input, output, session){
  RD <- reactiveValues(allData=RNAseq, selectedData=NULL)
  
  output$expressionData <- renderDT(RD$allData[, c(1:3, 4, 6, 8, 10, 12, 14)], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    colnames=c("Gene_ID", "Gene_name", "Gene_family", "WT_FPKM", "Ref5_FPKM", "Ref2_FPKM", "Med5_FPKM", "Ref5Med5_FPKM", "Ref2Med5_FPKM"),
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="7.5%", targets=c(4, 5, 6, 7, 8, 9)), list(width="15%", targets=c(1, 2))),
      pageLength=15,
      searchCols=default_search_columns,
      stateSave=FALSE)
  )
  
  proxy <- dataTableProxy("expressionData")
  
  observeEvent(input$expressionData_cell_clicked, {
    rowClicked <- input$expressionData_cell_clicked$row
    if(is.null(rowClicked)){
      return()
    }else{
      RD$selectedData <- rbind(RD$selectedData, RD$allData[rowClicked, ])
      RD$allData <- RD$allData[-rowClicked, ]
      
      isolate({
        default_search_columns <- c("", input$expressionData_search_columns)
        proxy %>% updateSearch(keywords=list(columns=default_search_columns))
      })
    }
  })

  output$selectedGenes <- renderDT(RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    colnames=c("Gene_ID", "Gene_name", "Gene_family"),
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="320px",
      searchCols=default_search_columns,
      stateSave=FALSE
    )
  )
  
  observeEvent(input$selectedGenes_cell_clicked, {
    rowClicked <- input$selectedGenes_cell_clicked$row
    if(is.null(rowClicked)){
      return()
    }else{
      RD$allData <- rbind(RD$allData, RD$selectedData[rowClicked, ])
      RD$selectedData <- RD$selectedData[-rowClicked, ]
      
      isolate({
        default_search_columns <- c("", input$expressionData_search_columns)
        proxy %>% updateSearch(keywords=list(columns=default_search_columns))
      })
    }
  })

  output$genePlot <- renderPlot({
    if(!is.null(RD$selectedData)){
      WT_data <- RD$selectedData[ , 1:5]
      WT_data$Genotype <- "WT"
      colnames(WT_data)[4:5] <- c("FPKM", "stdError")
      
      Ref5_data <- RD$selectedData[ , c(1:3, 6:7)]
      Ref5_data$Genotype <- "Ref5"
      colnames(Ref5_data)[4:5] <- c("FPKM", "stdError")
      
      Ref2_data <- RD$selectedData[ , c(1:3, 8:9)]
      Ref2_data$Genotype <- "Ref2"
      colnames(Ref2_data)[4:5] <- c("FPKM", "stdError")
      
      Med5_data <- RD$selectedData[ , c(1:3, 10:11)]
      Med5_data$Genotype <- "Med5"
      colnames(Med5_data)[4:5] <- c("FPKM", "stdError")
      
      Ref5Med5_data <- RD$selectedData[ , c(1:3, 12:13)]
      Ref5Med5_data$Genotype <- "Ref5Med5"
      colnames(Ref5Med5_data)[4:5] <- c("FPKM", "stdError")
      
      Ref2Med5_data <- RD$selectedData[ , c(1:3, 14:15)]
      Ref2Med5_data$Genotype <- "Ref2Med5"
      colnames(Ref2Med5_data)[4:5] <- c("FPKM", "stdError")
      
      plotData <- rbind(WT_data, Ref5_data, Ref2_data, Med5_data, Ref5Med5_data, Ref2Med5_data)
      colnames(plotData)[1:3] <- c("Gene_ID", "Gene_name", "Gene_family")
      plotData$Genotype <- factor(as.character(plotData$Genotype), levels=c("WT", "Ref5", "Ref2", "Med5", "Ref5Med5", "Ref2Med5"))
      plotData$Gene_ID <- factor(as.character(plotData$Gene_ID), levels=RD$selectedData$gene_id)
      
      ggplot(data=plotData, mapping=aes(x=Gene_ID, y=FPKM, fill=Genotype)) +
        geom_col(position=position_dodge()) +
        geom_errorbar(mapping=aes(ymin=FPKM-stdError, ymax=FPKM+stdError), color="black",
                      width=0.2, position=position_dodge(0.9)) +
        theme_bw() +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14),
              legend.position="right")
    }else{
      return(NA)
    }
  })
  
}

shinyApp(ui=ui, server=server)
