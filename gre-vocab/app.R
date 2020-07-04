# Set working directory
setwd("C:/Users/Bryce/Documents/interactive-data/gre-vocab/")

# Load necessary libraries, functions, and data ----
packageList <- c("shiny", "shinyjs", "dplyr")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]

if(length(newPackages) > 0){
  install.packages(newPackages)
}

library(shiny)
library(shinyjs)
library(dplyr)

source("helpers.R")

vocab <- read.csv("data/gre-vocab.csv", header=TRUE)
colnames(vocab) <- c("word", "definition")

# Define button styles for correct and incorrect answers
correctBtnStyle <- "background-color:#47AD61; color:#FFFFFF; border-color:#000000; border-style:solid; border-width:1px; border-radius:5%; font-size:18px;"
incorectBtnStyle <- "background-color:#D45252; color:#FFFFFF; border-color:#000000; border-style:solid; border-width:1px; border-radius:5%; font-size:18px;"

# Define UI ----
ui <- fluidPage(
  useShinyjs(),
  titlePanel("GRE Vocab Practice"),
  fluidRow(
    column(3, 
      wellPanel(
       tags$p("Tool for practicing GRE vocabulary words."),
       tags$p("Match the word with its definition.")
      )
    ),
    column(9,
      textOutput("definition"),
      actionButton("option1", "Option1"),
      actionButton("option2", "Option2"),
      actionButton("Option3", "Option3"),
      actionButton("Option4", "Option4"),
      actionButton("Option5", "Option5")
    )
  )
)

server <- function(input, output, session){
  
}

shinyApp(ui=ui, server=server)