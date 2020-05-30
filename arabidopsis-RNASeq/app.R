# Interactive visualization of Arabidopsis RNASeq data
# Note: before running, working directory must be set to parent directory of app.R

# Load necessary libraries, functions, and data ----
packageList <- c("shiny", "shinyjs", "dplyr", "DT", "ggplot2")
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
colnames(RNAseq) <- c("Gene_ID", "Gene_name", "Gene_family", "WT_FPKM", "WT_FPKM_SE", 
  "Ref5_FPKM", "Ref5_FPKM_SE", "Ref2_FPKM", "Ref2_FPKM_SE", "Med5_FPKM", "Med5_FPKM_SE",
  "Ref5Med5_FPKM", "Ref5Med5_FPKM_SE", "Ref2Med5_FPKM", "Ref2Med5_FPKM_SE")

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Arabidopsis RNAseq data viewer"),
  fluidRow(
    column(3, 
      wellPanel(
        tags$p("Tool for interactive subsetting of Arabidopsis RNAseq data."),
        tags$p(tags$strong("To use:")),
        tags$ol(
          tags$li("Search and filter for genes in the lower table with the header row of the
                  table."),
          tags$li("Select genes to plot from the lower table by clicking on them."),
          tags$li("Genes being plotted will appear in the table to the right of the plot 
                  in the order that they were selected."),
          tags$li("To remove genes from the plot, click on them in the right table")
        )
      ),
      checkboxGroupInput("selectedGenotypes",
        label=tags$p(tags$strong("Select genotypes to plot:")),
        inline=FALSE,
        choices=list("WT"="WT", "Ref5"="Ref5", "Ref2"="Ref2", "Med5"="Med5", "Ref5Med5"="Ref5Med5", "Ref2Med5"="Ref2Med5"),
        selected=c("WT", "Ref5", "Ref2", "Med5", "Ref5Med5", "Ref2Med5"))
    ),
    column(6,
      plotOutput("genePlot", height="500px")
    ),
    column(3,
      div(style="float:right; margin-bottom:20px",
        downloadButton(style="display:inline; margin-right: 10px", "downloadGraph_pdf", "Download graph"),
        downloadButton(style="display:inline", "downloadData_csv", "Download data")
      ),
      DTOutput("selectedGenes")
    )
  ),
  fluidRow(style="padding-top:20px",
    column(12,
      DTOutput("expressionData")
    )
  )
)

if (!exists("default_search_columns")) default_search_columns <- NULL

# TODO: To add genotype checkbox functionality, just show/hide columns in lower dataframe, and
# add if statements to plotting call. Much easier than changing allData/selectedData

server <- function(input, output, session){

  RD <- reactiveValues(allData=RNAseq, selectedData=NULL)
  
  output$expressionData <- renderDT(RD$allData[, c(1:3, getGenotypeCols(input$selectedGenotypes))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="15%", targets=c(1, 2, 3))),
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
        proxy %>% updateSearch(keywords=list(columns=default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
    }
  })

  output$downloadGraph_pdf <- downloadHandler(
    filename=function(){paste("Arabidopsis_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(RD$selectedData, input$selectedGenotypes)
      plotData$Genotype <- factor(as.character(plotData$Genotype), levels=c("WT", "Ref5", "Ref2", "Med5", "Ref5Med5", "Ref2Med5"))
      plotData$Gene_ID <- factor(as.character(plotData$Gene_ID), levels=RD$selectedData$Gene_ID)
      print(ggplot(data=plotData, mapping=aes(x=Gene_ID, y=FPKM, fill=Genotype)) +
              geom_col(position=position_dodge()) +
              geom_errorbar(mapping=aes(ymin=FPKM-SE, ymax=FPKM+SE), color="black",
                            width=0.2, position=position_dodge(0.9)) +
              theme_bw() +
              theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    legend.position="right")
      )
      dev.off()
    }
  )
  
  output$downloadData_csv <- downloadHandler(
    filename=function(){paste("Arabidopsis_RNAseq", ".csv", sep="")},
    content=function(file){
      write.csv(RD$selectedData[, c(1:3, getGenotypeCols2(input$selectedGenotypes))], file, row.names=FALSE)
    }
  )
  
  output$selectedGenes <- renderDT(RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="390px",
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
        proxy %>% updateSearch(keywords=list(columns=default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
      
      if(nrow(RD$selectedData)==0){
        RD$selectedData <- NULL
      }
    }
  })

  output$genePlot <- renderPlot({
    if(!is.null(RD$selectedData)){
      
      plotData <- convertToTidy(RD$selectedData, input$selectedGenotypes)
      
      plotData$Genotype <- factor(as.character(plotData$Genotype), levels=c("WT", "Ref5", "Ref2", "Med5", "Ref5Med5", "Ref2Med5"))
      plotData$Gene_ID <- factor(as.character(plotData$Gene_ID), levels=RD$selectedData$Gene_ID)
      
      ggplot(data=plotData, mapping=aes(x=Gene_ID, y=FPKM, fill=Genotype)) +
        geom_col(position=position_dodge()) +
        geom_errorbar(mapping=aes(ymin=FPKM-SE, ymax=FPKM+SE), color="black",
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
