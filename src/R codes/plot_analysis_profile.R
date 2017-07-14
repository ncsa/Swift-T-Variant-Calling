suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(lubridate, quietly = TRUE)))
suppressWarnings(suppressMessages(library(forcats, quietly = TRUE)))

options(tibble.print_max = 10, tibble.print_min = 4) # if there are > n row, print the first m

file = commandArgs(trailingOnly = TRUE)[1]
# Or use: file = '../../../../Results/run42/delivery/docs/'

rawdata <- read_tsv(file)

data <- rawdata %>%
  mutate(Time = as_datetime(Time)) %>%
  separate(`App status`, c('App', 'Status'), sep = ' ') 

start <- data %>%
  filter(Status == 'start') %>%
  select(-Status) %>%
  purrr::set_names( c('start_App', 'start_time'))

end <- data %>%
  filter(Status == 'end') %>%
  select(-Status) %>%
  set_names( c('end_App', 'end_time'))

as_tibble(cbind(start, end)) %>%
  select(-end_App) %>%
  mutate(start_App = fct_reorder2(start_App, start_time, end_time)) %>% 
  ggplot(aes(start_App))  +
  coord_flip() + geom_linerange(aes(ymin = start_time, ymax = end_time)) + ggsave('myprofile.pdf')
