library(shiny)
library(shinyjs)

setwd("C:/Users/Bryce/Documents/interactive-data/RNAseq/")

library(dplyr)
library(DT)
library(ggplot2)

source("helpers.R")

PRJNA549285_data <- readRDS("data/PRJNA549285.rds")
PRJNA546251_data <- readRDS("data/PRJNA546251.rds")
PRJNA230565_data <- readRDS("data/PRJNA230565.rds")
PRJNA388948_data <- readRDS("data/PRJNA388948.rds")
Cyp79A2_data <- readRDS("data/Cyp79A2.rds")


ui <- navbarPage("RNAseq data viewer",
# PRJNA549285 ui ----
  tabPanel("PRJNA549285", # hy5
    useShinyjs(),
    fluidRow(
      column(3, 
        wellPanel(
          tags$p("Raw reads downloaded from bioproject PRJNA549285. Reads were prefiltered and trimmed 
          with trimmomatic/0.39, and then aligned to TAIR10 reference genome with 
          tophat/2.1.2. Expression counts were generated with subread/2.0.0, and edgeR/3.30.3 was
          used for differential expression analysis."),
          tags$p(tags$strong("To use:")),
          tags$ol(
            tags$li("Search and filter for genes in the lower table with the header row of the
              table."),
            tags$li("Select genes to plot from the lower table by clicking on them."),
            tags$li("Genes being plotted will appear in the table to the right of the plot 
              in the order that they were selected. To remove genes from the plot, click on them in the right table"),
            tags$li("Download a copy of the current graph or data with the with the buttons in the top
              right corner of the window.")
          )
        ),
        checkboxGroupInput("PRJNA549285_selectedSamples",
          label=tags$p(tags$strong("Select samples to plot:")),
          inline=FALSE,
          choices=list("Col-0 dark 3 days (Col0_d3d)"="Col0_d3d",
                       "hy5 dark 3 days (hy5_d3d)"="hy5_d3d",
                       "Col-0 dark 3 days + light 1.5 hours (Col0_d3d_l1.5h)"="Col0_d3d_l1.5h",
                       "hy5 dark 3 days + light 1.5 hours (hy5_d3d_l1.5h)"="hy5_d3d_l1.5h"),
          selected=c("Col0_d3d", "hy5_d3d", "Col0_d3d_l1.5h", "hy5_d3d_l1.5h"))
        ),
        column(6,
          plotOutput("PRJNA549285_genePlot", height="580px")
        ),
        column(3,
          div(style="float:right; margin-bottom:20px",
            downloadButton(style="display:inline; margin-right: 10px", "PRJNA549285_RNAseq.pdf", "Download graph"),
            downloadButton(style="display:inline; margin-top: 10px", "PRJNA549285_RNAseq.csv", "Download data")
        ),
        DTOutput("PRJNA549285_selectedGenes")
      )
    ),
    fluidRow(style="padding-top:15px; padding-bottom:5px",
      column(12,
        DTOutput("PRJNA549285_expressionData")
      )
    )
  ), 
  
# PRJNA546251 ui ----
  tabPanel("PRJNA546251",
    useShinyjs(),
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
      checkboxGroupInput("PRJNA546251_selectedSamples",
        label=tags$p(tags$strong("Select samples to plot:")),
        inline=FALSE,
        choices=list("Col-0 white light (Col0_WL)"="Col0_WL",
                    "uvr8 white light (uvr8_WL)"="uvr8_WL",
                    "Col-0 white light + UVB 6 hours (Col0_WL_UVB6h)"="Col0_WL_UVB6h",
                    "uvr8 white light + UVB 6 hours (uvr8_WL_UVB6h)"="uvr8_WL_UVB6h"),
        selected=c("Col0_WL", "uvr8_WL", "Col0_WL_UVB6h", "uvr8_WL_UVB6h"))
      ),
      column(6,
        plotOutput("PRJNA546251_genePlot", height="580px")
      ),
      column(3,
        div(style="float:right; margin-bottom:20px",
          downloadButton(style="display:inline; margin-right: 10px", "PRJNA546251_RNAseq.pdf", "Download graph"),
          downloadButton(style="display:inline; margin-top: 10px", "PRJNA546251_RNAseq.csv", "Download data")
        ),
        DTOutput("PRJNA546251_selectedGenes")
      )
    ),
    fluidRow(style="padding-top:15px; padding-bottom:5px",
      column(12,
        DTOutput("PRJNA546251_expressionData")
      )
    )    
  ),
  
# PRJNA230565 ui ---- 
  tabPanel("PRJNA230565", # NAA treated
    useShinyjs(),
    fluidRow(
      column(3, 
        wellPanel(
          tags$p("Raw reads downloaded from bioproject PRJNA230565 Reads were prefiltered and trimmed 
            with trimmomatic/0.39, and then aligned to TAIR10 reference genome with 
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
        selectInput("PRJNA230565_unit", "Select expression unit:", 
          choices=list("TMM"="tmm",
                       "FPKM"="fpkm")
        ),
        checkboxGroupInput("PRJNA230565_selectedSamples",
          label=tags$p(tags$strong("Select samples to plot:")),
          inline=FALSE,
          choices=list("control"="control",
                       "NAA_treated"="NAA"),
          selected=c("control", "NAA")
        )
      ),
      column(6,
        plotOutput("PRJNA230565_genePlot", height="580px")
      ),
      column(3,
        div(style="float:right; margin-bottom:20px",
          downloadButton(style="display:inline; margin-right: 10px", "PRJNA230565_RNAseq.pdf", "Download graph"),
          downloadButton(style="display:inline; margin-top: 10px", "PRJNA230565_RNAseq.csv", "Download data")
        ),
        DTOutput("PRJNA230565_selectedGenes")
      )
    ),
    fluidRow(style="padding-top:15px; padding-bottom:5px",
      column(12,
        DTOutput("PRJNA230565_expressionData")
      )
    )    
  ), 
# PRJNA388948 ui ----
  tabPanel("PRJNA388948",
  useShinyjs(),
  fluidRow(
    column(3, 
      wellPanel(
        tags$p("Raw reads downloaded from bioproject PRJNA388948 Reads were prefiltered and trimmed 
          with trimmomatic/0.39, and then aligned to TAIR10 reference genome with 
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
      selectInput("PRJNA388948_unit", "Select expression unit:", 
        choices=list("TMM"="tmm",
                     "FPKM"="fpkm")
      ),
      checkboxGroupInput("PRJNA388948_selectedSamples",
        label=tags$p(tags$strong("Select samples to plot:")),
        inline=FALSE,
        choices=list("WT"="WT",
                     "ref5"="ref5",
                     "ref2"="ref2"),
        selected=c("WT", "ref5", "ref2")
      )
    ),
    column(6,
      plotOutput("PRJNA388948_genePlot", height="580px")
    ),
    column(3,
      div(style="float:right; margin-bottom:20px",
        downloadButton(style="display:inline; margin-right: 10px", "PRJNA388948_RNAseq.pdf", "Download graph"),
        downloadButton(style="display:inline; margin-top: 10px", "PRJNA388948_RNAseq.csv", "Download data")
      ),
      DTOutput("PRJNA388948_selectedGenes")
      )
    ),
    fluidRow(style="padding-top:15px; padding-bottom:5px",
      column(12,
         DTOutput("PRJNA388948_expressionData")
      )
    )    
  ), 

# Cyp79A2 ui ----
  tabPanel("Cyp79A2",
    useShinyjs(),
    fluidRow(
      column(3, 
        wellPanel(
          tags$p("Raw reads were prefiltered and trimmed 
            with trimmomatic/0.39, and then aligned to TAIR10 reference genome with 
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
        selectInput("Cyp79A2_unit", "Select expression unit:", 
          choices=list("TMM"="tmm",
                       "FPKM"="fpkm")
          ),
        checkboxGroupInput("Cyp79A2_selectedSamples",
           label=tags$p(tags$strong("Select samples to plot:")),
             inline=FALSE,
             choices=list("WT"="WT",
                          "WT_PAOx"="PAOx",
                          "Cyp79A2"="Cyp79A2"),
             selected=c("WT", "PAOx", "Cyp79A2")
        )
      ),
      column(6,
        plotOutput("Cyp79A2_genePlot", height="580px")
      ),
      column(3,
        div(style="float:right; margin-bottom:20px",
          downloadButton(style="display:inline; margin-right: 10px", "Cyp79A2_RNAseq.pdf", "Download graph"),
          downloadButton(style="display:inline; margin-top: 10px", "Cyp79A2_RNAseq.csv", "Download data")
        ),
        DTOutput("Cyp79A2_selectedGenes")
      )
    ),
    fluidRow(style="padding-top:15px; padding-bottom:5px",
      column(12,
         DTOutput("Cyp79A2_expressionData")
      )
    )           
  )
)
# ----
if (!exists("PRJNA549285_default_search_columns")) PRJNA549285_default_search_columns <- NULL
if (!exists("PRJNA546251_default_search_columns")) PRJNA546251_default_search_columns <- NULL
if (!exists("PRJNA230565_default_search_columns")) PRJNA230565_default_search_columns <- NULL
if (!exists("PRJNA388948_default_search_columns")) PRJNA388948_default_search_columns <- NULL
if (!exists("Cyp79A2_default_search_columns")) Cyp79A2_default_search_columns <- NULL

server <- function(input, output, session){
# PRJNA549285 server ----
  PRJNA549285_RD <- reactiveValues(allData=PRJNA549285_data, selectedData=NULL)
  
  output$PRJNA549285_expressionData <- renderDT(
    PRJNA549285_RD$allData[, c(1:4, getGenotypeCols(input$PRJNA549285_selectedSamples, PRJNA="PRJNA549285"))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="10%", targets=c(1, 2)), list(width="25%", targets=3), list(width="15%", targets=4)),
      pageLength=15,
      searchCols=PRJNA549285_default_search_columns,
      stateSave=FALSE)
  )
  
  PRJNA549285_proxy <- dataTableProxy("PRJNA549285_expressionData")
  
  observeEvent(input$PRJNA549285_expressionData_cell_clicked, {
    PRJNA549285_rowClicked <- input$PRJNA549285_expressionData_cell_clicked$row
    if(is.null(PRJNA549285_rowClicked)){
      return()
    }else{
      PRJNA549285_RD$selectedData <- rbind(PRJNA549285_RD$selectedData, PRJNA549285_RD$allData[PRJNA549285_rowClicked, ])
      PRJNA549285_RD$allData <- PRJNA549285_RD$allData[-PRJNA549285_rowClicked, ]
      
      isolate({
        PRJNA549285_default_search_columns <- c("", input$PRJNA549285_expressionData_search_columns)
        PRJNA549285_proxy %>% updateSearch(keywords=list(columns=PRJNA549285_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
    }
  })
  
  output$PRJNA549285_RNAseq.pdf <- downloadHandler(
    filename=function(){paste("PRJNA549285_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(PRJNA549285_RD$selectedData, input$PRJNA549285_selectedSamples, PRJNA="PRJNA549285")
      plotData$sample <- factor(as.character(plotData$sample), levels=c("Col0_d3d", "hy5_d3d", "Col0_d3d_l1.5h", "hy5_d3d_l1.5h"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA549285_RD$selectedData$short_name)
      print(ggplot(data=plotData, mapping=aes(x=short_name, y=RPKM, fill=sample)) +
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
  
  output$PRJNA549285_RNAseq.csv <- downloadHandler(
    filename=function(){paste("PRJNA549285_RNAseq", ".csv", sep="")},
    content=function(file){
      write.csv(PRJNA549285_RD$selectedData[, c(1:3, getGenotypeCols2(input$PRJNA549285_selectedSamples, PRJNA="PRJNA549285"))], file, row.names=FALSE)
    }
  )
  
  output$PRJNA549285_selectedGenes <- renderDT(PRJNA549285_RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="470px",
      searchCols=PRJNA549285_default_search_columns,
      stateSave=FALSE
    )
  )
  
  observeEvent(input$PRJNA549285_selectedGenes_cell_clicked, {
    PRJNA549285_rowClicked <- input$PRJNA549285_selectedGenes_cell_clicked$row
    if(is.null(PRJNA549285_rowClicked)){
      return()
    }else{
      PRJNA549285_RD$allData <- rbind(PRJNA549285_RD$allData, PRJNA549285_RD$selectedData[PRJNA549285_rowClicked, ])
      PRJNA549285_RD$selectedData <- PRJNA549285_RD$selectedData[-PRJNA549285_rowClicked, ]
      
      isolate({
        PRJNA549285_default_search_columns <- c("", input$PRJNA549285_expressionData_search_columns)
        PRJNA549285_proxy %>% updateSearch(keywords=list(columns=PRJNA549285_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
      
      if(nrow(PRJNA549285_RD$selectedData)==0){
        PRJNA549285_RD$selectedData <- NULL
      }
    }
  })
  
  output$PRJNA549285_genePlot <- renderPlot({
    if(!is.null(PRJNA549285_RD$selectedData)){
      
      plotData <- convertToTidy(PRJNA549285_RD$selectedData, input$PRJNA549285_selectedSamples, PRJNA="PRJNA549285")
      
      plotData$sample <- factor(as.character(plotData$sample), levels=c("Col0_d3d", "hy5_d3d", "Col0_d3d_l1.5h", "hy5_d3d_l1.5h"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA549285_RD$selectedData$short_name)
      
      ggplot(data=plotData, mapping=aes(x=short_name, y=RPKM, fill=sample)) +
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
  
# PRJNA546251 server ----
  PRJNA546251_RD <- reactiveValues(allData=PRJNA546251_data, selectedData=NULL)
  
  output$PRJNA546251_expressionData <- renderDT(PRJNA546251_RD$allData[, c(1:4, getGenotypeCols(input$PRJNA546251_selectedSamples, PRJNA="PRJNA546251"))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="10%", targets=c(1, 2)), list(width="25%", targets=3), list(width="15%", targets=4)),
      pageLength=15,
      searchCols=PRJNA546251_default_search_columns,
      stateSave=FALSE)
  )
  
  PRJNA546251_proxy <- dataTableProxy("PRJNA546251_expressionData")
  
  observeEvent(input$PRJNA546251_expressionData_cell_clicked, {
    PRJNA546251_rowClicked <- input$PRJNA546251_expressionData_cell_clicked$row
    if(is.null(PRJNA546251_rowClicked)){
      return()
    }else{
      PRJNA546251_RD$selectedData <- rbind(PRJNA546251_RD$selectedData, PRJNA546251_RD$allData[PRJNA546251_rowClicked, ])
      PRJNA546251_RD$allData <- PRJNA546251_RD$allData[-PRJNA546251_rowClicked, ]
      
      isolate({
        PRJNA546251_default_search_columns <- c("", input$PRJNA546251_expressionData_search_columns)
        PRJNA546251_proxy %>% updateSearch(keywords=list(columns=PRJNA546251_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
    }
  })
  
  output$PRJNA546251_RNAseq.pdf <- downloadHandler(
    filename=function(){paste("PRJNA546251_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(PRJNA546251_RD$selectedData, input$PRJNA546251_selectedSamples, PRJNA="PRJNA546251")
      plotData$sample <- factor(as.character(plotData$sample), levels=c("Col0_WL", "uvr8_WL", "Col0_WL_UVB6h", "uvr8_WL_UVB6h"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA546251_RD$selectedData$short_name)
      print(ggplot(data=plotData, mapping=aes(x=short_name, y=RPKM, fill=sample)) +
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
      write.csv(PRJNA546251_RD$selectedData[, c(1:3, getGenotypeCols2(input$PRJNA546251_selectedSamples, PRJNA="PRJNA546251"))], file, row.names=FALSE)
    }
  )
  
  output$PRJNA546251_selectedGenes <- renderDT(PRJNA546251_RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="470px",
      searchCols=PRJNA546251_default_search_columns,
      stateSave=FALSE
    )
  )
  
  observeEvent(input$PRJNA546251_selectedGenes_cell_clicked, {
    PRJNA546251_rowClicked <- input$PRJNA546251_selectedGenes_cell_clicked$row
    if(is.null(PRJNA546251_rowClicked)){
      return()
    }else{
      PRJNA546251_RD$allData <- rbind(PRJNA546251_RD$allData, PRJNA546251_RD$selectedData[PRJNA546251_rowClicked, ])
      PRJNA546251_RD$selectedData <- PRJNA546251_RD$selectedData[-PRJNA546251_rowClicked, ]
      
      isolate({
        PRJNA546251_default_search_columns <- c("", input$PRJNA546251_expressionData_search_columns)
        PRJNA546251_proxy %>% updateSearch(keywords=list(columns=PRJNA546251_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
      
      if(nrow(PRJNA546251_RD$selectedData)==0){
        PRJNA546251_RD$selectedData <- NULL
      }
    }
  })
  
  output$PRJNA546251_genePlot <- renderPlot({
    if(!is.null(PRJNA546251_RD$selectedData)){
      
      plotData <- convertToTidy(PRJNA546251_RD$selectedData, input$PRJNA546251_selectedSamples, PRJNA="PRJNA546251")
      
      plotData$sample <- factor(as.character(plotData$sample), levels=c("Col0_WL", "uvr8_WL", "Col0_WL_UVB6h", "uvr8_WL_UVB6h"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA546251_RD$selectedData$short_name)
      
      ggplot(data=plotData, mapping=aes(x=short_name, y=RPKM, fill=sample)) +
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
  
# PRJNA230565 server ----
  PRJNA230565_RD <- reactiveValues(allData=PRJNA230565_data, selectedData=NULL)
  
  output$PRJNA230565_expressionData <- renderDT(PRJNA230565_RD$allData[, c(1:4, getGenotypeCols(input$PRJNA230565_selectedSamples, PRJNA="PRJNA230565", unit=input$PRJNA230565_unit))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="10%", targets=c(1, 2)), list(width="25%", targets=3), list(width="15%", targets=4)),
      pageLength=15,
      searchCols=PRJNA230565_default_search_columns,
      stateSave=FALSE)
  )
  
  PRJNA230565_proxy <- dataTableProxy("PRJNA230565_expressionData")
  
  observeEvent(input$PRJNA230565_expressionData_cell_clicked, {
    PRJNA230565_rowClicked <- input$PRJNA230565_expressionData_cell_clicked$row
    if(is.null(PRJNA230565_rowClicked)){
      return()
    }else{
      PRJNA230565_RD$selectedData <- rbind(PRJNA230565_RD$selectedData, PRJNA230565_RD$allData[PRJNA230565_rowClicked, ])
      PRJNA230565_RD$allData <- PRJNA230565_RD$allData[-PRJNA230565_rowClicked, ]
      
      isolate({
        PRJNA230565_default_search_columns <- c("", input$PRJNA230565_expressionData_search_columns)
        PRJNA230565_proxy %>% updateSearch(keywords=list(columns=PRJNA230565_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
    }
  })
  
  output$PRJNA230565_RNAseq.pdf <- downloadHandler(
    filename=function(){paste("PRJNA230565_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(PRJNA230565_RD$selectedData, input$PRJNA230565_selectedSamples, PRJNA="PRJNA230565", unit=input$PRJNA230565_unit)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("control", "NAA"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA230565_RD$selectedData$short_name)
      print(ggplot(data=plotData, mapping=aes(x=short_name, y=exp, fill=sample)) +
              geom_col(position=position_dodge()) +
              theme_bw() +
              labs(y=input$PRJNA230565_unit) +
              theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    legend.position="right")
      )
      dev.off()
    }
  )
  
  output$PRJNA230565_RNAseq.csv <- downloadHandler(
    filename=function(){paste("PRJNA230565_RNAseq", ".csv", sep="")},
    content=function(file){
      write.csv(PRJNA230565_RD$selectedData[, c(1:3, getGenotypeCols2(input$PRJNA230565_selectedSamples, PRJNA="PRJNA230565", unit=input$PRJNA230565_unit))], file, row.names=FALSE)
    }
  )
  
  output$PRJNA230565_selectedGenes <- renderDT(PRJNA230565_RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="470px",
      searchCols=PRJNA230565_default_search_columns,
      stateSave=FALSE
    )
  )
  
  observeEvent(input$PRJNA230565_selectedGenes_cell_clicked, {
    PRJNA230565_rowClicked <- input$PRJNA230565_selectedGenes_cell_clicked$row
    if(is.null(PRJNA230565_rowClicked)){
      return()
    }else{
      PRJNA230565_RD$allData <- rbind(PRJNA230565_RD$allData, PRJNA230565_RD$selectedData[PRJNA230565_rowClicked, ])
      PRJNA230565_RD$selectedData <- PRJNA230565_RD$selectedData[-PRJNA230565_rowClicked, ]
      
      isolate({
        PRJNA230565_default_search_columns <- c("", input$PRJNA230565_expressionData_search_columns)
        PRJNA230565_proxy %>% updateSearch(keywords=list(columns=PRJNA230565_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
      
      if(nrow(PRJNA230565_RD$selectedData)==0){
        PRJNA230565_RD$selectedData <- NULL
      }
    }
  })
  
  output$PRJNA230565_genePlot <- renderPlot({
    if(!is.null(PRJNA230565_RD$selectedData)){
      plotData <- convertToTidy(PRJNA230565_RD$selectedData, input$PRJNA230565_selectedSamples, PRJNA="PRJNA230565", unit=input$PRJNA230565_unit)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("control", "NAA"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA230565_RD$selectedData$short_name)
      ggplot(data=plotData, mapping=aes(x=short_name, y=exp, fill=sample)) +
        geom_col(position=position_dodge()) +
        theme_bw() +
        labs(y=input$PRJNA230565_unit) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14),
              legend.position="right")
    }else{
      return(NA)
    }
  })
  
# PRJNA388948 server ----
  PRJNA388948_RD <- reactiveValues(allData=PRJNA388948_data, selectedData=NULL)
  
  output$PRJNA388948_expressionData <- renderDT(PRJNA388948_RD$allData[, c(1:4, getGenotypeCols(input$PRJNA388948_selectedSamples, PRJNA="PRJNA388948", unit=input$PRJNA388948_unit))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="10%", targets=c(1, 2)), list(width="25%", targets=3), list(width="15%", targets=4)),
      pageLength=15,
      searchCols=PRJNA388948_default_search_columns,
      stateSave=FALSE
    )
  )
  
  PRJNA388948_proxy <- dataTableProxy("PRJNA388948_expressionData")
  
  observeEvent(input$PRJNA388948_expressionData_cell_clicked, {
    PRJNA388948_rowClicked <- input$PRJNA388948_expressionData_cell_clicked$row
    if(is.null(PRJNA388948_rowClicked)){
      return()
    }else{
      PRJNA388948_RD$selectedData <- rbind(PRJNA388948_RD$selectedData, PRJNA388948_RD$allData[PRJNA388948_rowClicked, ])
      PRJNA388948_RD$allData <- PRJNA388948_RD$allData[-PRJNA388948_rowClicked, ]
      
      isolate({
        PRJNA388948_default_search_columns <- c("", input$PRJNA388948_expressionData_search_columns)
        PRJNA388948_proxy %>% updateSearch(keywords=list(columns=PRJNA388948_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
    }
  })
  
  output$PRJNA388948_RNAseq.pdf <- downloadHandler(
    filename=function(){paste("PRJNA388948_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(PRJNA388948_RD$selectedData, input$PRJNA388948_selectedSamples, PRJNA="PRJNA388948", unit=input$PRJNA388948_unit)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("WT", "ref5", "ref2"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA388948_RD$selectedData$short_name)
      print(ggplot(data=plotData, mapping=aes(x=short_name, y=exp, fill=sample)) +
              geom_col(position=position_dodge()) +
              theme_bw() +
              labs(y=input$PRJNA388948_unit) +
              theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    legend.position="right")
      )
      dev.off()
    }
  )
  
  output$PRJNA388948_RNAseq.csv <- downloadHandler(
    filename=function(){paste("PRJNA388948_RNAseq", ".csv", sep="")},
    content=function(file){
      write.csv(PRJNA388948_RD$selectedData[, c(1:3, getGenotypeCols2(input$PRJNA388948_selectedSamples, PRJNA="PRJNA388948", unit=input$PRJNA388948_unit))], file, row.names=FALSE)
    }
  )
  
  output$PRJNA388948_selectedGenes <- renderDT(PRJNA388948_RD$selectedData[, c(1, 2, 3)],
    server=TRUE, class="compact hover", selection="none",
    options=list(
      ordering=FALSE,
      paging=FALSE,
      dom="rti",
      autoWidth=TRUE,
      scrollY="470px",
      searchCols=PRJNA388948_default_search_columns,
      stateSave=FALSE
    )
  )
  
  observeEvent(input$PRJNA388948_selectedGenes_cell_clicked, {
    PRJNA388948_rowClicked <- input$PRJNA388948_selectedGenes_cell_clicked$row
    if(is.null(PRJNA388948_rowClicked)){
      return()
    }else{
      PRJNA388948_RD$allData <- rbind(PRJNA388948_RD$allData, PRJNA388948_RD$selectedData[PRJNA388948_rowClicked, ])
      PRJNA388948_RD$selectedData <- PRJNA388948_RD$selectedData[-PRJNA388948_rowClicked, ]
      
      isolate({
        PRJNA388948_default_search_columns <- c("", input$PRJNA388948_expressionData_search_columns)
        PRJNA388948_proxy %>% updateSearch(keywords=list(columns=PRJNA388948_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
      
      if(nrow(PRJNA388948_RD$selectedData)==0){
        PRJNA388948_RD$selectedData <- NULL
      }
    }
  })
  
  output$PRJNA388948_genePlot <- renderPlot({
    if(!is.null(PRJNA388948_RD$selectedData)){
      plotData <- convertToTidy(PRJNA388948_RD$selectedData, input$PRJNA388948_selectedSamples, PRJNA="PRJNA388948", unit=input$PRJNA388948_unit)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("WT", "ref5", "ref2"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=PRJNA388948_RD$selectedData$short_name)
      ggplot(data=plotData, mapping=aes(x=short_name, y=exp, fill=sample)) +
        geom_col(position=position_dodge()) +
        theme_bw() +
        labs(y=input$PRJNA388948_unit) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14),
              legend.position="right")
    }else{
      return(NA)
    }
  })
  
# Cyp79A2 server ----
  Cyp79A2_RD <- reactiveValues(allData=Cyp79A2_data, selectedData=NULL)
  
  output$Cyp79A2_expressionData <- renderDT(Cyp79A2_RD$allData[, c(1:4, getGenotypeCols(input$Cyp79A2_selectedSamples, PRJNA="Cyp79A2", unit=input$Cyp79A2_unit))], 
    server=TRUE, filter="top", selection="none", class="compact hover",
    options=list(
      sDom="<'top'>rt<'bottom'>ip",
      #autoWidth=TRUE,
      scrollX=TRUE,
      columnDefs=list(list(width="10%", targets=c(1, 2)), list(width="25%", targets=3), list(width="15%", targets=4)),
      pageLength=15,
      searchCols=Cyp79A2_default_search_columns,
      stateSave=FALSE
    )
  )
  
  Cyp79A2_proxy <- dataTableProxy("Cyp79A2_expressionData")
  
  observeEvent(input$Cyp79A2_expressionData_cell_clicked, {
    Cyp79A2_rowClicked <- input$Cyp79A2_expressionData_cell_clicked$row
    if(is.null(Cyp79A2_rowClicked)){
      return()
    }else{
      Cyp79A2_RD$selectedData <- rbind(Cyp79A2_RD$selectedData, Cyp79A2_RD$allData[Cyp79A2_rowClicked, ])
      Cyp79A2_RD$allData <- Cyp79A2_RD$allData[-Cyp79A2_rowClicked, ]
      
      isolate({
        Cyp79A2_default_search_columns <- c("", input$Cyp79A2_expressionData_search_columns)
        Cyp79A2_proxy %>% updateSearch(keywords=list(columns=Cyp79A2_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
    }
  })
  
  output$Cyp79A2_RNAseq.pdf <- downloadHandler(
    filename=function(){paste("Cyp79A2_RNAseq", ".pdf", sep="")},
    content=function(file){pdf(file, height=6, width=12)
      plotData <- convertToTidy(Cyp79A2_RD$selectedData, input$Cyp79A2_selectedSamples, PRJNA="Cyp79A2", unit=input$Cyp79A2_unit)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("WT", "PAOx", "Cyp79A2"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=Cyp79A2_RD$selectedData$short_name)
      print(ggplot(data=plotData, mapping=aes(x=short_name, y=exp, fill=sample)) +
              geom_col(position=position_dodge()) +
              theme_bw() +
              labs(y=input$Cyp79A2_unit) +
              theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14),
                    legend.position="right")
      )
      dev.off()
    }
  )
  
  output$Cyp79A2_RNAseq.csv <- downloadHandler(
    filename=function(){paste("Cyp79A2_RNAseq", ".csv", sep="")},
    content=function(file){
      write.csv(Cyp79A2_RD$selectedData[, c(1:3, getGenotypeCols2(input$Cyp79A2_selectedSamples, PRJNA="Cyp79A2", unit=input$Cyp79A2_unit))], file, row.names=FALSE)
    }
  )
  
  output$Cyp79A2_selectedGenes <- renderDT(Cyp79A2_RD$selectedData[, c(1, 2, 3)],
     server=TRUE, class="compact hover", selection="none",
     options=list(
       ordering=FALSE,
       paging=FALSE,
       dom="rti",
       autoWidth=TRUE,
       scrollY="470px",
       searchCols=Cyp79A2_default_search_columns,
       stateSave=FALSE
     )
  )
  
  observeEvent(input$Cyp79A2_selectedGenes_cell_clicked, {
    Cyp79A2_rowClicked <- input$Cyp79A2_selectedGenes_cell_clicked$row
    if(is.null(Cyp79A2_rowClicked)){
      return()
    }else{
      Cyp79A2_RD$allData <- rbind(Cyp79A2_RD$allData, Cyp79A2_RD$selectedData[Cyp79A2_rowClicked, ])
      Cyp79A2_RD$selectedData <- Cyp79A2_RD$selectedData[-Cyp79A2_rowClicked, ]
      
      isolate({
        Cyp79A2_default_search_columns <- c("", input$Cyp79A2_expressionData_search_columns)
        Cyp79A2_proxy %>% updateSearch(keywords=list(columns=Cyp79A2_default_search_columns)) %>% reloadData(resetPaging=FALSE)
      })
      
      if(nrow(Cyp79A2_RD$selectedData)==0){
        Cyp79A2_RD$selectedData <- NULL
      }
    }
  })
  
  output$Cyp79A2_genePlot <- renderPlot({
    if(!is.null(Cyp79A2_RD$selectedData)){
      plotData <- convertToTidy(Cyp79A2_RD$selectedData, input$Cyp79A2_selectedSamples, PRJNA="Cyp79A2", unit=input$Cyp79A2_unit)
      plotData$sample <- factor(as.character(plotData$sample), levels=c("WT", "PAOx", "Cyp79A2"))
      plotData$short_name <- factor(as.character(plotData$short_name), levels=Cyp79A2_RD$selectedData$short_name)
      ggplot(data=plotData, mapping=aes(x=short_name, y=exp, fill=sample)) +
        geom_col(position=position_dodge()) +
        theme_bw() +
        labs(y=input$Cyp79A2_unit) +
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