#
# This is the user-interface definition of the visualization Shiny web application. 
#

library(shiny)

if (!require(shinydashboard)) {
  install.packages('shinydashboard')
  library(shinydashboard)
}

if (!require(plotly)){
  install.packages('plotly')
  library(plotly)
}

######
dashboardPage(
  dashboardHeader( title = "Swift/T Genomic Variant Calling workflow provenance trace",
                   titleWidth = 450,
                   tags$li(a(
                     img(src = 'logo_swift.png',
                         title = "Swift/T manual", height = "30px"),
                     href = 'http://www.swift-lang.org/Swift-T/',
                     style = "padding-top:10px; padding-bottom:10px;"),
                     class = "dropdown"),
                   tags$li(a(
                     img(src = 'github-mark.png',
                         title = "Workflow code", height = "30px"),
                     href = 'https://github.com/ncsa/Swift-T-Variant-Calling',
                     style = "padding-top:10px; padding-bottom:10px;"),
                     class = "dropdown"),
                   tags$li(a(
                     img(src = 'readthedocs.png',
                         title = "Workflow manual", height = "30px"),
                     href = 'http://swift-t-variant-calling.readthedocs.io/en/latest/',
                     style = "padding-top:10px; padding-bottom:10px;"),
                     class = "dropdown")
                   
                   ),
  
  
  dashboardSidebar(
    # tags$head(tags$style(".wrapper {overflow: visible !important;}")),
    br(),
    br(),
    fileInput(inputId = "logfile",
              label = "Timing trace file:",
              buttonLabel = "Upload ... "
    ),
    br(),
    
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
    )

    # br(),
    # img(src = 'logo_uofk.png', height = 50, width = 50),
    # img(src = 'logo_uiuc.png', height = 50, width = 50),
    # img(src = 'logo_swift.png', height = 50, width = 80),
    # br(),
    # br(),
    # 'This App analyzes a timing trace coming from',
    # a('The Genomic Swift/T variant calling workflow',
    #   href = 'https://github.com/jacobrh91/Swift-T-Variant-Calling')
    # , 'which was developed in collaboration between these institutes.'
    # 
  ),
  
  dashboardBody(
    # Boxes need to be put in a row (or column)
    tabsetPanel(type = "tabs",
                tabPanel("Overall Summary",
                         br(),
                         h3("Number of samples processed by each application and stage:"),
                         tableOutput("Samples_Summary"),
                         
                         verbatimTextOutput("Chromosomes_Summary"),
                         
                         h3("Details of chromosome processed per sample:"),
                         tableOutput("Chromosomes_table")
               ),
                
                tabPanel("Provenance Plot", 
                         tags$head(tags$style(HTML(' 
                            .form-group, .selectize-control {
                                margin-bottom: 0px;
                            }
                            .box-body {
                                padding-bottom: 0px;
                                padding-top: 10px;
                            }'))),
                    
                         box( status = "primary",  solidHeader = TRUE,
                              br(),
                              br(),
                              plotlyOutput("PlotlyProvenancePlot", width = "200%"),
                              # plotOutput("simplePlot", click = clickOpts(id ="plot_click")),
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
                                         min = 1, max = 60, value = 30), width = 3),
                         box(
                           width = 3, solidHeader = T, status = "primary", align = "center",
                           a(
                             img(src = 'logo_uofk.png', height = 50, width = 50, align = "center",
                                 title = "University of Khartoum"),
                             href = 'http://cbsb.uofk.edu/'
                           ),
                           a(
                             img(src = 'logo_uiuc.png', height = 50, width = 50, align = "center",
                                 title = "University of Illinois Urbana-Champaign"),
                             href = 'https://wiki.ncsa.illinois.edu/display/LH/HPC+for+Computational+Genomics'
                           ),
                           
                           br(),
                           br(),
                           'This App analyzes a timing trace coming from',
                           a('The Genomic Swift/T variant calling workflow',
                             href = 'https://github.com/ncsa/Swift-T-Variant-Calling')
                           , 'which was developed in collaboration between centres in these institutes.'
                         )
                         
                ),
                
                tabPanel("Provenance Table", 
                         br(),
                         tableOutput("provenanceTable"))
                
                
    )
  )
)
