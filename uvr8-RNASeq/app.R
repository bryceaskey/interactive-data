# Interactive visualization of Arabidopsis RNASeq data
# Note: before running, working directory must be set to parent directory of app.R

# Load necessary libraries, functions, and data ----
#packageList <- c("shiny", "shinyjs", "dplyr", "DT", "ggplot2")
#newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
#if(length(newPackages) > 0){
#  install.packages(newPackages)
#}

library(shiny)
library(shinyjs)
library(dplyr)
library(DT)
library(ggplot2)

source("helpers.R")

load("data/allData.rda")

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  titlePanel("Arabidopsis RNAseq data viewer - PRJNA546251"),
  fluidRow(
    column(3, 
      wellPanel(
        tags$p("Raw reads downloaded from bioproject PRJNA546251. Reads were prefiltered and trimmed 
               with adapterremoval/2.2.2, and then aligned to TAIR10 reference genome with 
               hisat2/2.2.0. Expression counts were generated with subread/2.0.0, and edgeR/3.30.3 was
               used for differential expression analysis."),
        tags$p(tags$strong("To use:")),
        tags$ol(
          tags$li("Search and filter for genes in the lower table with the header row of the
                  table."),
          tags$li("Select genes to plot from the lower table by clicking on them."),
          tags$li("Genes being plotted will appear in the table to the right of the plot 
                  in the order that they were selected. To remove genes from the plot, click on them
                  in the right table"),
          tags$li("Download a copy of the current graph or data with the with the buttons in the top
                  right corner of the window.")
        )
      ),
      checkboxGroupInput("selectedSamples",
        label=tags$p(tags$strong("Select samples to plot:")),
        inline=FALSE,
        choices=list("Col-0 white light (Col0_WL)"="Col0_WL",
                     "Col-0 white light + UVB 6 hours (Col0_WL_UVB6h)"="Col0_WL_UVB6h",
                     "uvr8 white light (uvr8_WL)"="uvr8_WL",
                     "uvr8 white light + UVB 6 hours (uvr8_WL_UVB6h)"="uvr8_WL_UVB6h"),
        selected=c("Col0_WL", "Col0_WL_UVB6h", "uvr8_WL", "uvr8_WL_UVB6h"))
    ),
    column(6,
      plotOutput("genePlot", height="580px")
    ),
    column(3,
      div(style="float:right; margin-bottom:20px",
        downloadButton(style="display:inline; margin-right: 10px", "PRJNA546251_RNAseq.pdf", "Download graph"),
        downloadButton(style="display:inline; margin-top: 10px", "PRJNA546251_RNAseq.csv", "Download data")
      ),
      DTOutput("selectedGenes")
    )
  ),
  fluidRow(style="padding-top:15px; padding-bottom:5px",
    column(12,
      DTOutput("expressionData")
    )
  )
)

if (!exists("default_search_columns")) default_search_columns <- NULL

server <- function(input, output, session){

  RD <- reactiveValues(allData=allData, selectedData=NULL)
  
  output$expressionData <- renderDT(RD$allData[, c(1:3, getGenotypeCols(input$selectedSamples))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="15%", targets=c(1, 2)), list(width="30%", targets=3)),
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

  output$PRJNA546251_RNAseq.pdf <- downloadHandler(
    filename=function(){paste("PRJNA546251_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(RD$selectedData, input$selectedSamples)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("Col0_WL", "Col0_WL_UVB6h", "uvr8_WL", "uvr8_WL_UVB6h"))
      plotData$locus <- factor(as.character(plotData$locus), levels=RD$selectedData$locus) ###############################################
      print(ggplot(data=plotData, mapping=aes(x=locus, y=RPKM, fill=sample)) +
              geom_col(position=position_dodge()) +
              geom_errorbar(mapping=aes(ymin=RPKM-SE, ymax=RPKM+SE), color="black",
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
  
  output$PRJNA546251_RNAseq.csv <- downloadHandler(
    filename=function(){paste("PRJNA546251_RNAseq", ".csv", sep="")},
    content=function(file){
      write.csv(RD$selectedData[, c(1:3, getGenotypeCols2(input$selectedSamples))], file, row.names=FALSE)
    }
  )
  
  output$selectedGenes <- renderDT(RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="470px",
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
      
      plotData <- convertToTidy(RD$selectedData, input$selectedSamples)
      
      plotData$sample <- factor(as.character(plotData$sample), levels=c("Col0_WL", "Col0_WL_UVB6h", "uvr8_WL", "uvr8_WL_UVB6h"))
      plotData$locus <- factor(as.character(plotData$locus), levels=RD$selectedData$locus)
      
      ggplot(data=plotData, mapping=aes(x=locus, y=RPKM, fill=sample)) +
        geom_col(position=position_dodge()) +
        geom_errorbar(mapping=aes(ymin=RPKM-SE, ymax=RPKM+SE), color="black",
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
