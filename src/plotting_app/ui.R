#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.


library(shiny)

# Define UI for application that draws a histogram
shinyUI( fluidPage(
  
  # Application title
  titlePanel("Swift/T Genomic Variant Calling workflow provenance trace"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
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
      ),
     
      br(),
      br(),
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
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Provenance Plot", 
                           # p('Showing provenance plot for Sample'), 
                           # span(textOutput('textSample'), style = "color:blue"),
                           # p('and Chromosome'),
                           # span(textOutput('textChromosome'), style = "color:blue"),
                           br(),
                           plotOutput("provenancePlot"),
                           downloadButton('saveFig', 'Save figure', style="float:right")
                           ),
                           
                  tabPanel("Provenance Table", 
                           br(),
                           tableOutput("provenanceTable"))
                  
                  
      )
      
    )
  )
)

)