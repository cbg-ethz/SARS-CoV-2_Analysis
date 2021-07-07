input_file <- snakemake@input$fname
filter_file <- snakemake@input$fname_samples

filter_flag <- TRUE # flag for whether to filter
ent_div <-  FALSE #TRUE # flag for diversity or entropy
output_file <- snakemake@output$fname

# entropy function
Entropy <- function(v) {
  ifelse(v==0, 0, -v*log(v))
}

# diversity function
Diversity <- function(v, n) {
  v*(v*n - 1)/(n - 1)
}

# read data
covid_table <- read.table(input_file, header = TRUE)

library(tidyverse)

# how many samples?
covid_table %>% pull(SAMPLE) %>% unique() %>% length()

if (filter_flag) { # filter samples
  # read samples to filter
  samples_to_keep <- read.csv(filter_file, header = TRUE)
  covid_table %>% filter(SAMPLE %in% samples_to_keep$sample) -> covid_table
  # how many samples now?
  covid_table %>% pull(SAMPLE) %>% unique() %>% length()
}

if (ent_div) { # compute entropy
  covid_table %>%
  mutate(ENT = Entropy(A_freq) + Entropy(C_freq) +
               Entropy(G_freq) + Entropy(T_freq) +
               Entropy(DEL_freq) + Entropy(IN_freq)) %>%
  select(SAMPLE, POS, ENT) -> ent_table
} else { # compute diversity
  covid_table %>%
  mutate(N = X.ADJUSTED_.READ_COUNT,
         DIV = 1 - Diversity(A_freq, N) - Diversity(C_freq, N) -
             Diversity(G_freq, N) - Diversity(T_freq, N) -
             Diversity(DEL_freq, N) - Diversity(IN_freq, N)) %>%
  select(SAMPLE, POS, DIV) -> ent_table
}

write.table(ent_table, file = output_file, row.names = FALSE, quote = FALSE)
