
samples_df <- read.csv(snakemake@input$fname_entropy_samples)

# take top 10

samples_df <- samples_df[order(samples_df$LOG_ENT, decreasing = TRUE), ][1:10, ]

mean_percent <- function(x) {
  round(mean(x)*100, 2)
}

### Process data

library(tidyverse)

# read data
covid_table <- read.table(snakemake@input$fname_vcf, header = TRUE)

n_pos <- 29903

covid_table %>% filter(SAMPLE %in% samples_df$ID) -> covid_table

covid_table %>% group_by(SAMPLE) %>%
  summarise(mut = n()) %>%
  mutate(mut = round(mut/n_pos*100, 2)) -> covid_table

left_join(samples_df, covid_table, by = c("ID" = "SAMPLE")) %>%
  mutate(log_ent = round(LOG_ENT, 2)) %>% select(ID, log_ent, mut) -> samples_table

cbind(samples_table[1:5, ], samples_table[6:10, ]) %>%
  write.table(snakemake@output$fname, sep = ",", quote = FALSE, row.names = F)
