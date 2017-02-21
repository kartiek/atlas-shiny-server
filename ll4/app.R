#
# This is a Demo Shiny web application for performing log-logistic regession. You can run the application by clicking
# the 'Run App' button above.
#
#Author: https://github.com/kartiek


library(shiny)
library(drc)

# Define UI for application that accepts a delimited file and perform 4p log-logistic regression
ui <- fluidPage(
  titlePanel("4p logistic regression"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose INPUT File',
                accept=c('text/csv', 
                         'text/comma-separated-values,text/plain', 
                         '.csv')),
      tags$hr(),
      checkboxInput('header', 'Header', TRUE),
      radioButtons('sep', 'Separator',
                   c(Comma=',',
                     Tab='\t'),
                   ','),
      radioButtons('quote', 'Quote',
                   c(None='',
                     'Double Quote'='"',
                     'Single Quote'="'"),
                   '"'),
      numericInput("obs", "Enter OD value for prediction", 1, min=0, max=100),
      submitButton("Predict!!!")
    ),
    mainPanel(
      h4('Plotting model'),
      plotOutput('plot1'),
      
      h4('Summary'),
      verbatimTextOutput('summary'),
      
      h4('Prediction'),
      tableOutput('view')
      
    )
  )
)

# Define server logic required to read a file and do ll4
server <- function(input, output, session) {
  
  dataModel <- reactive({
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)}
    # dRed <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                     # quote=input$quote)
    dRed <- readr::read_delim(inFile$datapath, col_names = input$header, delim = input$sep,
                       quote = input$quote)
    # colnames(dRed) <- c('OD','CONC')
    # dRed$CONC <- log10(dRed$CONC)
    dRed <- dplyr::mutate(dRed, CONC = log10(CONC))
    dRed <- dplyr::filter(dRed, !is.infinite(CONC))
    dMod <- drm(OD~CONC, fct=LL.4(),data=dRed)
    dMod
  })
  
  output$plot1 <- renderPlot({
    if (is.null(dataModel())) {
      return(NULL)
    }
    plot(dataModel())
  })
  
  output$summary <- renderPrint({
    if (is.null(dataModel())) {
      return(NULL)
    }
    summary(dataModel())
  })
  
  output$view <- renderTable({
    if (is.null(dataModel()))
      return(NULL)
    as.data.frame(ED(dataModel(),input$obs,type = 'absolute',logBase = 10))
    # as.data.frame(ED(dataModel(),input$obs,type = 'absolute'))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)