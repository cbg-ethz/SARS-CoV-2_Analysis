input_file1 <- snakemake@input$fname_persample
input_file2 <- snakemake@input$fname_persamplemuts

output_file <- snakemake@output$fname

samples_df <- read.table(input_file1, header = TRUE)
muts_df <- read.table(input_file2, header = TRUE)

# take top 10

samples_df <- samples_df[order(samples_df$log_ENT, decreasing = TRUE), ][1:10, ]

mean_percent <- function(x) {
  round(mean(x)*100, 2)
}

### Process data

library(tidyverse)

n_pos <- 29903

left_join(samples_df, muts_df) %>%
  mutate(mut = round(exp(log_MUT)/n_pos*100, 2),
         log_ent = round(log_ENT, 2)) %>%
  select(SAMPLE, log_ent, mut) %>%
  write.table(output_file, quote = FALSE, row.names = FALSE)
