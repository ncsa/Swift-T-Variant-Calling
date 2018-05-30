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


if (!require(stringr)){
  install.packages('stringr')
  library(stringr)
}
#  Defining server logic --------------------------------------------------

shinyServer(function(input, output, session) {
  
  rawdataInput <- reactive ({
    inFile <- input$logfile
    if (is.null(inFile))
      return(NULL)
    readr::read_tsv(inFile$datapath) %>% unique()
  })
  
  dataInput <- eventReactive (input$logfile, {
  data <- rawdataInput() %>%
    mutate(Time = as_datetime(Time)) %>%
    separate(`App status`, c('App', 'Status'), sep = ' ')

  start <- data %>%
    filter(Status == 'start') %>%
    select(-Status) %>%
    purrr::set_names( c('Sample', 'Chromosome', 'start_App', 'start_time', 'start_Stage'))

  end <- data %>%
    filter(Status == 'end') %>%
    select(-Status, -Sample, -Chromosome) %>%
    set_names( c('end_App', 'end_time', 'end_Stage'))

  data <- cbind(start, end) %>%
    mutate(Application = start_App, Stage = start_Stage) %>%
    select(-end_App, -start_App, -start_Stage, -end_Stage) %>%
    mutate(Application = fct_reorder2(Application, start_time, end_time)) %>%
    group_by(Stage, Sample, Chromosome) %>% 
    arrange(end_time) %>%
    mutate(start_time =  lag(end_time, default = first(start_time))) %>%
    ungroup() %>%
    mutate(start_time =  if_else(str_detect(Application, "GenotypeGVCFs"), 
                                 max(lag(end_time, default = first(start_time))), start_time))%>%
    mutate(Application = fct_reorder2(Application, start_time, end_time)) 
  })
  
  dataPlot <- reactive({
    data <- dataInput() 
    hour_res <- input$hourResValue
    minute_res <- input$minResValue
    day_res <- input$dayResValue
    
    if (input$minResol == T) resol = paste(minute_res, "min") else resol = paste(hour_res, "hour")
    # if (input$daysRes == T) resol = paste(day_res, "day") else resol = paste(hour_res, "hour")
    
    if (input$zoomSample ==T )  data <- data %>% filter(Sample == input$sample)
    if (input$zoomChromo == T)  data <- data %>% filter(Chromosome == input$chromosome)
    if (input$zoomApp ==T )     data <- data %>% filter(Application == input$app)
    
    plot <- data %>% 
      ggplot(aes(Application))  +
      geom_linerange(aes(ymin = start_time, ymax = end_time, color = Sample, text = Chromosome),
                     position = position_dodge(width = 0.5), size = 2) +
      coord_flip() + scale_y_datetime(breaks = date_breaks(resol),
                                      labels = date_format("%F\n %H:%M")
                                      # expand = expand_scale(mult = 0.1) 
                                      ) +
      labs(x = "Application", y = "Timeline") +
      theme(axis.text.x = element_text(angle=45), 
            axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=16),
            axis.line = element_line(colour = "darkblue", 
                                     size = 1, linetype = "solid"))
    
  })

  output$textSample <- renderText({
    input$sample
  })
  output$textChromosome <- renderText({
    input$chromosome
  })
  

  output$provenanceTable <- renderTable({
    dataInput() %>% 
      mutate(start_time = as.character(start_time)) %>%
      mutate(end_time = as.character(end_time)) %>%
      select(Stage, Application, Sample, Chromosome, start_time, end_time)
  })

  output$PlotlyProvenancePlot <- renderPlotly({
    ggplotly(dataPlot(), tooltip = c("colour", "text", "x"))
  })
  
  output$Run_Summary <- renderTable(
    dataInput()  %>% arrange(desc(start_time)) %>%group_by(Application) %>% 
      summarize(Processed_Samples = n())
  )
  
  # output$simplePlot <- renderPlot(
  #   print(dataPlot())
  # )
  
  # output$info <- renderText({
  #  
  #   if(!is.null(input$plot_click)){
  #     data <- dataInput() 
  #     hover <- input$plot_click
  #     
  #     point <- nearPoints(data, hover, yvar = "start_time", xvar = "Application",
  #                         threshold = 5, maxpoints = 1, addDist = TRUE)
  #     if (nrow(point) == 0) return(NULL)
  #     
  #     
  #     paste0("Sample=", point, 
  #            "\nChromosome = ", point,
  #            "\nApplication = ")
  #   }
  #   
  #  
  # })
  
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
  


  # output$saveFig <- downloadHandler(
  #   filename = function() {
  #     paste(input$logfile , Sys.time(), '.png', sep='')},
  #   content = function(file) {
  #     ggsave(file,  plot = dataPlot(), device = 'png')
  #   }
  # )
  
})  



