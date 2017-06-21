#
# This is a shiny app for plotting heatmaps from gene expression data.
#
# Author: https://github.com/kartiek
#

library(shiny)
library(heatmaply)
library(shinythemes)
library(magrittr)

# Read the data
expData <- readRDS('ThCell.rds')
expGenes <- unique(expData$geneSymbol)

# Define UI for application 
ui <- fluidPage(theme = shinytheme("cerulean"),
                # Application title
                titlePanel("Generate heatmap"),
                tags$head(tags$style(HTML(".multicol{height:auto;
                                          -webkit-column-count: 2;
                                          -moz-column-count: 2;
                                          column-count: 2;
                                          }
                                          div.checkbox {margin-top: 0px;}"))),
                sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      'tPlot',h4('Input type'),
                      choices = c('Type', 'Upload')),
                    
                    uiOutput("ui"),
                    h4('Heatmap options'),
                    checkboxInput('rowDend','Row dendrogram',FALSE),
                    checkboxInput('colDend','Column dendrogram',FALSE),
                    radioButtons('sclP','Scaling',choices = c('Row'='row','Column'='column','None'='none'),inline = TRUE,selected = 'none'),
                    h4('Select time points'), 
                    tags$div(align = "left", 
                             class = "multicol",
                             checkboxGroupInput(
                               'timeGroup', label = NULL,
                               choices = c('0h (Thp)' = '0H', '0.5h' = '05H', '1h' = '1H',
                                           '2h' = '2H', '4h' = '4H', '6h' = '6H',
                                           '12h' = '12H', '24h' = '24H', '48h' = '48H',
                                           '72h' = '72H'),
                               selected = c('0H','05H','1H','2H','4H','6H','12H','24H','48H','72H')
                             )),
                    h4('Select samples'),
                    checkboxGroupInput(
                      'thGroup', label = NULL, #h4('Th1 & Th2'),
                      choices = c('Th1' = 'Th1',
                                  'Th2' = 'Th2',
                                  'Th0-012' = 'Th0-012'),
                      selected = c('Th1','Th2','Th0-012'),inline = TRUE),
                    checkboxGroupInput(
                      'th17Group', label = NULL, #h4('Th17'),
                      choices = c('Th17' = 'Th17',
                                  'Th0-17' = 'Th0-17'),
                      selected = c('Th17','Th0-17'),inline = TRUE),
                    checkboxGroupInput(
                      'tregGroup', label = NULL, #h4('iTreg'),
                      choices = c('iTreg' = 'iTreg',
                                  'Th0-iTr' = 'Th0-iTr'),
                      selected = c('iTreg','Th0-iTr'),inline = TRUE),
                    actionButton("go", "Plot", icon('chevron-right')),
                    width = 3
                  ),
                  mainPanel(
                    plotlyOutput("plot1"), width = 9
                  )
                )
                )

# Serve the data
server <- function(input, output, session) {
  # Note the start time
  anDt <- dplyr::data_frame(START = Sys.time())
  
  output$ui <- renderUI({
    if (is.null(input$tPlot))
      return()
    switch(input$tPlot,
           'Type' = selectizeInput(
             'genes', h4('Type your gene(s) of choice'),
             choices = NULL, multiple = TRUE),
           'Upload' = fileInput('file1', 'Choose INPUT File',
                                accept=c('text/csv', 
                                         'text/comma-separated-values,text/plain', 
                                         '.csv')))
  })
  # Updating selectize input
  observe({
  updateSelectizeInput(session, 'genes', choices = expGenes, server = TRUE)})
  
  #Reactive expression to get the filtered data 
  expDataTemp <- eventReactive(input$go, {
    if(input$tPlot == 'Type'){
      fTest <- dplyr::filter(expData,geneSymbol %in% input$genes & timP %in% input$timeGroup &
                               subset %in% c(input$thGroup,input$th17Group,input$tregGroup)) %>% 
        tidyr::unite(sub_tim,subset,timP,sep='_') %>% 
        dplyr::select(sub_tim,geneSymbol,meanE)
      fTest$sub_tim <- factor(fTest$sub_tim,levels = c('Th0-012_0H','Th0-012_05H','Th0-012_1H','Th0-012_2H','Th0-012_4H','Th0-012_6H','Th0-012_12H','Th0-012_24H','Th0-012_48H','Th0-012_72H',
                                                     'Th1_0H','Th1_05H','Th1_1H','Th1_2H','Th1_4H','Th1_6H','Th1_12H','Th1_24H','Th1_48H','Th1_72H',
                                                     'Th2_0H','Th2_05H','Th2_1H','Th2_2H','Th2_4H','Th2_6H','Th2_12H','Th2_24H','Th2_48H','Th2_72H',
                                                     'Th0-17_0H','Th0-17_05H','Th0-17_1H','Th0-17_2H','Th0-17_4H','Th0-17_6H','Th0-17_12H','Th0-17_24H','Th0-17_48H','Th0-17_72H',
                                                     'Th17_0H','Th17_05H','Th17_1H','Th17_2H','Th17_4H','Th17_6H','Th17_12H','Th17_24H','Th17_48H','Th17_72H',
                                                     'Th0-iTr_0H','Th0-iTr_05H','Th0-iTr_1H','Th0-iTr_2H','Th0-iTr_4H','Th0-iTr_6H','Th0-iTr_12H','Th0-iTr_24H','Th0-iTr_48H','Th0-iTr_72H',
                                                     'iTreg_0H','iTreg_05H','iTreg_1H','iTreg_2H','iTreg_4H','iTreg_6H','iTreg_12H','iTreg_24H','iTreg_48H','iTreg_72H'))
      fTest <- fTest %>% tidyr::spread(sub_tim, meanE) %>% as.data.frame
      fTest <- fTest %>% magrittr::set_rownames(fTest$geneSymbol) %>% dplyr::select(-geneSymbol)
      fTest
    }
    else if(input$tPlot == 'Upload'){
      inFile <- input$file1
      fRed <- readr::read_csv(inFile$datapath, col_names = FALSE)
      fTest <- dplyr::filter(expData,geneSymbol %in% fRed$X1 & timP %in% input$timeGroup &
                               subset %in% c(input$thGroup,input$th17Group,input$tregGroup)) %>% 
        tidyr::unite(sub_tim,subset,timP,sep='_') %>% 
        dplyr::select(sub_tim,geneSymbol,meanE)
      fTest$sub_tim <- factor(fTest$sub_tim,levels = c('Th0-012_0H','Th0-012_05H','Th0-012_1H','Th0-012_2H','Th0-012_4H','Th0-012_6H','Th0-012_12H','Th0-012_24H','Th0-012_48H','Th0-012_72H',
                                                       'Th1_0H','Th1_05H','Th1_1H','Th1_2H','Th1_4H','Th1_6H','Th1_12H','Th1_24H','Th1_48H','Th1_72H',
                                                       'Th2_0H','Th2_05H','Th2_1H','Th2_2H','Th2_4H','Th2_6H','Th2_12H','Th2_24H','Th2_48H','Th2_72H',
                                                       'Th0-17_0H','Th0-17_05H','Th0-17_1H','Th0-17_2H','Th0-17_4H','Th0-17_6H','Th0-17_12H','Th0-17_24H','Th0-17_48H','Th0-17_72H',
                                                       'Th17_0H','Th17_05H','Th17_1H','Th17_2H','Th17_4H','Th17_6H','Th17_12H','Th17_24H','Th17_48H','Th17_72H',
                                                       'Th0-iTr_0H','Th0-iTr_05H','Th0-iTr_1H','Th0-iTr_2H','Th0-iTr_4H','Th0-iTr_6H','Th0-iTr_12H','Th0-iTr_24H','Th0-iTr_48H','Th0-iTr_72H',
                                                       'iTreg_0H','iTreg_05H','iTreg_1H','iTreg_2H','iTreg_4H','iTreg_6H','iTreg_12H','iTreg_24H','iTreg_48H','iTreg_72H'))
      fTest <- fTest %>% tidyr::spread(sub_tim, meanE) %>% as.data.frame 
      fTest <- fTest %>% magrittr::set_rownames(fTest$geneSymbol) %>% dplyr::select(-geneSymbol)
      fTest
    }})
  
  #Reactive expression with heatmaply plot
  ppt <- eventReactive(input$go, {
      p <- expDataTemp() %>% heatmaply(scale = input$sclP,Rowv = input$rowDend, Colv = input$colDend,plot_method = 'plotly')
  })
  
  output$plot1 <- renderPlotly({
    if (is.null(expDataTemp)) {
      return(NULL)
    }
    ppt()
  })
  
  # Note the end time and write to local directory
  session$onSessionEnded(function() {
    anDt <- dplyr::mutate(anDt, END = Sys.time())
    readr::write_tsv(anDt, 'heatGene_analytics.txt', append = TRUE)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)