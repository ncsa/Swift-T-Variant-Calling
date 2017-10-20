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

if (!require(plotly)){
  install.packages('plotly')
  library(plotly)
}

if (!require(scales)){
  install.packages('scales')
  library(scales)
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
    mutate(Application = start_App) %>%
    select(-end_App, -start_App) %>%
    mutate(Application = fct_reorder2(Application, start_time, end_time))
  })
  
  dataPlot <- reactive({
    data <- dataInput() 
    
    if (input$zoomSample ==T )  data <- data %>% filter(Sample == input$sample)
    if (input$zoomChromo == T)  data <- data %>% filter(Chromosome == input$chromosome)
    if (input$zoomApp ==T )     data <- data %>% filter(Application == input$app)
    
    plot <- data %>% 
      ggplot(aes(Application))  +
      geom_linerange(aes(ymin = start_time, ymax = end_time, color = Sample, text = Chromosome),
                     position = position_dodge(width = 0.5), size = 2) +
      coord_flip() + scale_y_datetime(breaks = date_breaks("3 hour"),
                                      labels = date_format("%F\n %H:%M"),
                                      expand = expand_scale(mult = 0.1) ) +
      theme(axis.text.x = element_text(angle=45))
    
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

  output$PlotlyProvenancePlot <- renderPlotly({
    ggplotly(dataPlot(), tooltip = c("colour", "text", "x"))
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
                      choices = unique(dataInput()$Application) )
  )
  

  output$saveFig <- downloadHandler(
    filename = function() {
      paste(input$logfile , Sys.time(), '.png', sep='')},
    content = function(file) {
      ggsave(file,  plot = dataPlot(), device = 'png')
    }
  )
  
})  



