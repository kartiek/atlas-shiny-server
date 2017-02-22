#
# This is a shiny app for plotting gene expression data.
#
# Author: https://github.com/kartiek
#

library(shiny)
library(tidyverse)
library(plotly)
library(shinythemes)

# Read the data
expData <- readRDS('ThCell.rds')
expGenes <- unique(expData$geneSymbol)

# Define UI for application 
ui <- fluidPage(theme = shinytheme("cerulean"),
   
   # Application title
   titlePanel("Plotting gene expression data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        selectizeInput(
          'genes', 'Type your gene(s) of choice', choices = NULL, multiple = TRUE
        ),
        hr(),
        checkboxGroupInput('timeGroup', label = h3('Select time points'), 
                           choices = c('0h' = '0H', '0.5h' = '05H', '1h' = '1H',
                                          '2h' = '2H', '4h' = '4H', '6h' = '6H',
                                          '12h' = '12H', '24h' = '24H', '48h' = '48H',
                                          '72h' = '72H'),
                           selected = c('0H','05H','1H','2H','4H','6H','12H','24H','48H','72H'),
                           inline = TRUE, width = '220px'), 
        hr(),
        checkboxGroupInput('cellGroup', label = h3('Select samples'), 
                           choices = c('Th1' = 'Th1', 'Th2' = 'Th2', 'Th17' = 'Th17',
                                          'iTreg' = 'iTreg'),
                           selected = c('Th1','Th2','Th17','iTreg'),inline = TRUE)),
      # Show a plot of the generated distribution
      mainPanel(
        plotlyOutput("plot1")
      )
   )
)

# Define server logic required to draw a histogram


# Serve the data
server <- function(input, output, session) {
  # Updating selectize input
  updateSelectizeInput(session, 'genes', choices = expGenes, server = TRUE) 
  
  #Reactive expression to get the filtered data 
  expDataTemp <- reactive({
    fTest <- filter(expData,geneSymbol %in% input$genes & timP %in% input$timeGroup &
                      subset %in% input$cellGroup) %>% 
      unite(sub_tim,subset,timP,sep='_')
    colnames(fTest) <- c('Sample','geneSymbol','maxExp','minExp','meanExp')
    fTest
  })
  #Reactive expression with plotly plot
  ppt <- reactive({
    p <- expDataTemp() %>% ggplot(data=.) + geom_line(aes(Sample,meanExp,group=1,col=geneSymbol)) +
      geom_ribbon(aes(x=Sample,ymin=minExp,ymax=maxExp,group=1,fill=geneSymbol),alpha=0.3) +
      theme_minimal() + labs(y='Expression',x=NULL) +
      theme(axis.text.x=element_text(angle=65, hjust=1),axis.ticks=element_blank())
    ggplotly(p,tooltip = c('Sample','meanExp','geneSymbol','maxExp','minExp'))
  })
  
  output$plot1 <- renderPlotly({
    if (is.null(expDataTemp)) {
      return(NULL)
    }
    ppt()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)