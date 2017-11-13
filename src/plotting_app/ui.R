#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.


library(shiny)
library(shinydashboard)
if (!require(plotly)){
  install.packages('plotly')
  library(plotly)
}

######
dashboardPage(
  dashboardHeader(title = "Swift/T Genomic Variant Calling workflow provenance trace",
                  titleWidth = 450),
  dashboardSidebar(
    fileInput(inputId = "logfile",
              label = "Timing trace file:",
              buttonLabel = "Upload ... "
    ),
    
    strong('Zoom-in options:'),
    
    checkboxInput('zoomSample', 'Zoom to sample', value = FALSE),
    selectInput(inputId = "sample", 
                label = NULL,
                choices = c(''),
                selected = "NULL"
    ),
    
    checkboxInput('zoomChromo', 'Zoom to chromosome', value = FALSE),
    selectInput(inputId = "chromosome", 
                label = NULL,
                choices = c(''),
                selected = "NULL"
    ),
    
    checkboxInput('zoomApp', 'Zoom to application', value = FALSE),
    selectInput(inputId = "app", 
                label = NULL,
                choices = c(''),
                selected = "NULL"
    ),

    br(),
    img(src = 'logo_uofk.png', height = 50, width = 50),
    img(src = 'logo_uiuc.png', height = 50, width = 50),
    img(src = 'logo_swift.png', height = 50, width = 80),
    br(),
    br(),
    'This App analyzes a timing trace coming from',
    a('The Genomic Swift/T variant calling workflow',
      href = 'https://github.com/jacobrh91/Swift-T-Variant-Calling')
    , 'which was developed in collaboration between these institutes.'
    
  ),
  
  dashboardBody(
    # Boxes need to be put in a row (or column)
    tabsetPanel(type = "tabs",
                tabPanel("Provenance Plot", 
                         # p('Showing provenance plot for Sample'), 
                         # span(textOutput('textSample'), style = "color:blue"),
                         # p('and Chromosome'),
                         # span(textOutput('textChromosome'), style = "color:blue"),
                    
                         box( status = "primary",  solidHeader = TRUE,
                              br(),
                              br(),
                              plotlyOutput("PlotlyProvenancePlot"),
                              downloadButton('saveFig', 'Save figure', style="float:right"),
                              width = 9),
                         
                         h4("Selet the proper timeline resolution (x-axis) in either:"),
                         # box(checkboxInput("daysRes", "Days", FALSE),
                         #     sliderInput(inputId = "dayResValue",
                         #                 label = NULL,
                         #                 min = 1, max = 15, step = 1, value = 1), width = 3),
                         box(solidHeader = TRUE,
                             checkboxInput("hoursRes", "Hours", TRUE),
                             sliderInput(inputId = "hourResValue",
                                         label = NULL,
                                         min = 1, max = 24, step = 1, value = 1), width = 3),
                         box(solidHeader = TRUE,
                             checkboxInput("minResol", "Minutes", FALSE),
                             sliderInput(inputId = "minResValue",
                                         label = NULL,
                                         min = 1, max = 60, value = 30), width = 3)
                         
                ),
                
                tabPanel("Provenance Table", 
                         br(),
                         tableOutput("provenanceTable"))
                
                
    )
  )
)
