#
# This is a Demo Shiny web application for exploring DE genes and their nearest enhancers. You can run the application by clicking
# the 'Run App' button above.
#
#Author: https://github.com/kartiek

library(shiny)
library(dplyr)
library(readr)
library(magrittr)
library(GenomicFeatures)
library(plotly)
library(shinythemes)

# Define UI for application
ui <- fluidPage(theme = shinytheme("cerulean"),
                
                # Application title
                titlePanel("Explore differentially expressed genes"),
                
                # Sidebar with a slider input for distance from gene
                verticalLayout(
                  plotlyOutput("plot1"),
                  hr(),
                  verbatimTextOutput("event"),
                  wellPanel(
                    h4("Filter"), 
                    sliderInput("distnc",
                                "Choose the distance(kb)",
                                min = 0,
                                max = 4500,
                                step = 5,
                                value = c(100,1800),
                                ticks = TRUE,
                                post = 'kb')
                    )
                  
                  # Show a plot of the selected genes
                  # mainPanel(
                    
                    # )
                  )
                )

# Define server logic required to draw a scatter plot

# Read the data
deGenes <- read_tsv('Th17DE.txt') %>% dplyr::filter(abs(logFC_72h) > 1)
cPeaks <- as.data.frame(read_tsv('Th17_72_LE_homer.txt'))
genLocs <- read_tsv('geneLocations.txt')
gL1 <- genLocs %>% dplyr::filter(geneSymbol %in% deGenes$geneSymbol)
gL1$chr <- paste0('chr',gL1$chr)
f1 <- with(gL1,GRanges(chr,IRanges(start = start,end = end)))
f2 <- with(cPeaks,GRanges(chr,IRanges(start = start,end = end)))
f3 <- as.data.frame(distanceToNearest(f1,f2))
expDis <- data_frame(geneSymbol=gL1$geneSymbol,dis=f3$distance) %>% 
  inner_join(deGenes,by='geneSymbol')

# Serve the data
server <- function(input, output, session) {
  # Note the start time
  anDt <- data_frame(START = Sys.time())
  
  # Reactive expression to get the filtered data 
  expDisTemp <- reactive({
    minDis <- (input$distnc[1])*1e3
    maxDis <- (input$distnc[2])*1e3
    expDis %>% filter(dis >= minDis,dis <= maxDis) %>% as.data.frame
  })
  # Reactive expression with plotly plot
  ppt <- reactive({
    expDisTemp() %>%
      plot_ly(x=expDisTemp()$logFC_72h,y=(expDisTemp()$dis/1e3),
              size = (expDisTemp()$dis/1e3),text = paste0('Gene Name: ',expDisTemp()$geneSymbol),
              type = 'scatter',mode='markers') %>%
      plotly::layout(title= 'Distance to nearest lineage-specific enhancer',
                     yaxis=list(zeroline=FALSE,title='distance(kb)'),
                     xaxis=list(zeroline=FALSE,title='logFC'))
  })
  output$plot1 <- renderPlotly({
    if (is.null(expDisTemp)) {
      return(NULL)
    }
    ppt()
  })
  output$event <- renderPrint({
    d <- event_data("plotly_hover")
    if (is.null(d)) "Hover on a point!" else d
  })
  
  # Note the end time and write to local directory
  session$onSessionEnded(function() {
    anDt <- anDt %>% mutate(END = Sys.time())
    write_tsv(anDt, 'xplore_analytics.txt', append = TRUE)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

