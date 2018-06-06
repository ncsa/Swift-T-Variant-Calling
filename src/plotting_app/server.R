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
    data <- rawdataInput()
    header <- c('Sample', 'Chromosome', 'App status', 'Time', 'Stage')
    if (!identical(names(data), header)){
      raw_1 <- as_tibble(matrix(data = names(data), nrow = 1)) %>% 
        setNames(names(data )) 
      raw_1[,4] <- as.integer(raw_1[,4])
      data <- bind_rows(data , raw_1) %>% setNames(header)
    }
    

    data <- data %>%
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
  
  output$tmp <- renderTable(dataInput())

  output$provenanceTable <- renderTable({
    dataInput() %>% 
      mutate(start_time = as.character(start_time)) %>%
      mutate(end_time = as.character(end_time)) %>%
      select(Stage, Application, Sample, Chromosome, start_time, end_time)
  })

  output$PlotlyProvenancePlot <- renderPlotly({
    ggplotly(dataPlot(), tooltip = c("colour", "text", "x"))
  })
  
  output$Samples_Summary <- renderTable(
    dataInput()  %>% 
      distinct(Sample, Stage, Application) %>%
      mutate(Stage = factor(Stage, levels = unique(Stage)), 
             Application = factor(Application, levels = unique(Application))) %>%
      group_by(Stage, Application) %>%  summarise(Processed_Samples = n()) %>%
      ungroup() %>%
      mutate(Stage = as.character(Stage)) %>%
      mutate(Stage = case_when(!duplicated(Stage) ~ Stage,
                               T ~ "")) 
  )
  
  chromosomes_summary <- eventReactive(input$logfile,{
    summary <- dataInput() %>%
      filter(Chromosome != "ALL") %>%
      mutate(Stage = factor(Stage, levels = unique(Stage)),
             Application = factor(Application, levels = unique(Application))) %>%
      group_by(Sample,Application) %>%
      summarise(Processed_chromosomes = n()) %>% ungroup() %>%
      mutate(Sample = case_when(!duplicated(Sample) ~ Sample,
                                T ~ ""))
  })
  
  output$Chromosomes_Summary <- renderText({
    if (length(unique(chromosomes_summary()$Processed_chromosomes)) == 1)
      msg <- "*All* chromosomes in *ALL* samples were processed successfully."
    else
      msg <- "*Some* chromosomes in *Some* samples were not processed- See table below for details"
    if (length(unique(chromosomes_summary()$Processed_chromosomes)) == 0)
      msg <- "*No* chromosomes were processed at all"
    msg
  })
  output$Chromosomes_table <- renderTable(
    chromosomes_summary() 
  )
  
  # output$simplePlot <- renderPlot(
  #   print(dataPlot())
  # )
  

  
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
  
  
})  



