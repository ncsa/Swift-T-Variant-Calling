#!/usr/bin/env Rscript

# This script accepts a `Timing.log` file that is either complete (from a successful pipeline)
# run, or a concatenation of temporary timing log files (in case of partial pipeline run). 

# Loading libraries --------------------------------------------------------
if (!require(tidyverse)){
  install.packages('tidyverse')
  library(tidyverse)
}

if (!require(lubridate)){
  install.packages('lubridate')
  library(lubridate)
}

if (!require(gdata)){
  install.packages('gdata')
  library(gdata)
}


# Clean-up of data --------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
inFile <- args[1]
rawdataInput <- readr::read_tsv(inFile) %>% unique()
header <- c('Sample', 'Chromosome', 'App status', 'Time', 'Stage')

if (!identical(names(rawdataInput), header)){
  raw_1 <- as_tibble(matrix(data = names(rawdataInput), nrow = 1)) %>% 
    setNames(names(rawdataInput )) 
  raw_1[,4] <- as.integer(raw_1[,4])
  rawdataInput <- bind_rows(rawdataInput , raw_1) %>% setNames(header)
}


rawdataInput <- rawdataInput %>%
  mutate(Time = as_datetime(Time)) %>%
  separate(`App status`, c('App', 'Status'), sep = ' ')

start <- rawdataInput %>%
  filter(Status == 'start') %>%
  select(-Status) %>%
  purrr::set_names( c('Sample', 'Chromosome', 'start_App', 'start_time', 'start_Stage'))

end <- rawdataInput %>%
  filter(Status == 'end') %>%
  select(-Status, -Sample, -Chromosome) %>%
  set_names( c('end_App', 'end_time', 'end_Stage'))

dataInput <- cbind(start, end) %>%
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



# Output tables and messages ----------------------------------------------

Samples_Summary <- dataInput  %>% 
  distinct(Sample, Stage, Application) %>%
  mutate(Stage = factor(Stage, levels = unique(Stage)), 
         Application = factor(Application, levels = unique(Application))) %>%
  group_by(Stage, Application) %>%  summarise(Processed_Samples = n()) %>%
  ungroup() %>%
  mutate(Stage = as.character(Stage)) %>%
  mutate(Stage = case_when(!duplicated(Stage) ~ Stage, T ~ "")) 

chromosomes_summary <- dataInput %>%
  filter(Chromosome != "ALL") %>%
  mutate(Stage = factor(Stage, levels = unique(Stage)),
         Application = factor(Application, levels = unique(Application))) %>%
  group_by(Sample,Application) %>%
  summarise(Processed_chromosomes = n()) %>% ungroup() %>%
  mutate(Sample = case_when(!duplicated(Sample) ~ Sample, T ~ ""))

Chromosomes_msg <- ifelse(length(unique(chromosomes_summary$Processed_chromosomes)) == 1,
                          "*All* chromosomes in *ALL* samples were processed successfully.",
                          ifelse(length(unique(chromosomes_summary$Processed_chromosomes)) == 0,
                                 "*No* chromosomes were processed at all",
                                 "*Some* chromosomes in *Some* samples were not processed- See table below for details"))


# Generation of the output file -------------------------------------------

output_file <- "pipeline_run_summary.txt"
write_file("Number of samples processed by each application and stage:\n\n", output_file)
write.fwf(x = as.data.frame(Samples_Summary), output_file, append = T, sep = '|')

write_file("\n\n\t*********************************************************", output_file, append = T)
write_file("\n\nDetails of chromosome processed per sample:\n\t", output_file, append = T)
write_file(Chromosomes_msg, path = output_file, append = T)

write_file("\n\n", output_file, append = T)
write.fwf(as.data.frame(chromosomes_summary), output_file, append = T, sep = '|')

