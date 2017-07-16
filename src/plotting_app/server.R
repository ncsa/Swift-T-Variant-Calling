# Loading libraries --------------------------------------------------------

library(shiny)
if (!require(tidyverse)){
  install.packages('tidyverse')
  library(tidyverse)
}

if (!require(lubridate)){
  install.packages('lubridate')
  library(lubridate)
}

if (!require(forcats)){
  install.packages('forcats')
  library(forcats)
}

#  Defining server logic --------------------------------------------------

shinyServer(function(input, output, session) {
  
  rawdataInput <- reactive ({
    inFile <- input$logfile
    if (is.null(inFile))
      return(NULL)
    readr::read_tsv(inFile$datapath)
  })
  
  dataInput <- eventReactive (input$logfile, {
  data <- rawdataInput() %>%
    mutate(Time = as_datetime(Time)) %>%
    separate(`App status`, c('App', 'Status'), sep = ' ')

  start <- data %>%
    filter(Status == 'start') %>%
    select(-Status) %>%
    purrr::set_names( c('Sample', 'Chromosome', 'start_App', 'start_time'))

  end <- data %>%
    filter(Status == 'end') %>%
    select(-Status, -Sample, -Chromosome) %>%
    set_names( c('end_App', 'end_time'))

  data <- cbind(start, end) %>%
    select(-end_App) %>%
    mutate(start_App = fct_reorder2(start_App, start_time, end_time))
  })
  
  dataPlot <- reactive({
    data <- dataInput() 
    
    if (input$zoomSample ==T )  data <- data %>% filter(Sample == input$sample)
    if (input$zoomChromo == T)  data <- data %>% filter(Chromosome == input$chromosome)
    if (input$zoomApp ==T )     data <- data %>% filter(start_App == input$app)
    
    plot <- data %>% 
      ggplot(aes(start_App))  +
      geom_linerange(aes(ymin = start_time, ymax = end_time, color = Sample)) +
      coord_flip()
  })
  
  output$textSample <- renderText({
    input$sample
  })
  output$textChromosome <- renderText({
    input$chromosome
  })
  

  output$provenanceTable <- renderTable({
    rawdataInput()
  })
  output$provenancePlot <- renderPlot({
    print(dataPlot())
     # dataInput() %>%
     #  ggplot(aes(start_App))  +
     #  geom_linerange(aes(ymin = start_time, ymax = end_time, color = Sample)) +
     #  coord_flip()
  })
  
  observeEvent(input$zoomSample,
    updateSelectInput(session, "sample",
                      choices = unique(dataInput()$Sample) )
    )
  
  observeEvent(input$zoomChromo,
    updateSelectInput(session, "chromosome",
                      choices = unique(dataInput()$Chromosome) )
  )
  
  observeEvent(input$zoomApp,
    updateSelectInput(session, "app",
                      choices = unique(dataInput()$start_App) )
  )
  

  output$saveFig <- downloadHandler(
    filename = function() {
      paste(input$logfile , Sys.time(), '.png', sep='')},
    content = function(file) {
      ggsave(file,  plot = dataPlot(), device = 'png')
    }
  )
  
})  



