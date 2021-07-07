filter_flag <- FALSE # flag for whether to filter
input_file1 <- snakemake@input@fname
input_file2 <- snakemake@input@fname_entropy
filter_file <- snakemake@input$samples

output_file <- snakemake@output@fname

base_df <- read.table(input_file2, header = TRUE)

# read data
covid_table <- read.table(input_file1, header = TRUE)

library(tidyverse)

# how many samples?
covid_table %>%
  pull(SAMPLE) %>%
  unique() %>%
  length()

if (filter_flag) { # filter samples
  # read samples to filter
  samples_to_keep <- read.csv(filter_file, header = TRUE)
  covid_table %>% filter(SAMPLE %in% samples_to_keep$sample) -> covid_table
  # how many samples now?
  covid_table %>%
    pull(SAMPLE) %>%
    unique() %>%
    length()
}

# entropy function
Entropy <- function(v) {
  ifelse(v == 0, 0, -v * log(v))
}

mean_percent <- function(x) {
  round(mean(x) * 100, 2)
}

mean_percent_sub <- function(x) { # only include non-zero values
  round(mean(x[which(x > 0)]) * 100, 2)
}

n_samps <- covid_table %>%
  pull(SAMPLE) %>%
  unique() %>%
  length()

covid_table %>% filter(POS %in% c(base_df$position, 23403)) -> covid_table

covid_table %>%
  filter(POS %in% c(23403)) %>%
  group_by(POS) %>%
  summarise(mut = n(), A_cons = sum(A_freq > 0.5)) %>%
  mutate(A_percent = round(A_cons / n_samps * 100, 2), G_percent = 100 - A_percent) %>%
  pull(G_percent)

covid_table %>% filter(POS %in% c(base_df$position)) -> covid_table

covid_table %>%
  mutate(ENT = Entropy(A_freq) + Entropy(C_freq) +
    Entropy(G_freq) + Entropy(T_freq) +
    Entropy(DEL_freq) + Entropy(IN_freq)) %>%
  group_by(POS) %>%
  summarise(
    mut = n(), ent = sum(ENT > 0), ref = names(sort(table(REF_BASE), decreasing = TRUE)[1]),
    A = mean_percent(A_freq), C = mean_percent(C_freq),
    G = mean_percent(G_freq), T = mean_percent(T_freq),
    del = mean_percent(DEL_freq), ins = mean_percent(IN_freq),
    A_frac = mean(A_freq > 0), C_frac = mean(C_freq > 0),
    G_frac = mean(G_freq > 0), T_frac = mean(T_freq > 0),
    del_frac = mean(DEL_freq > 0), in_frac = mean(IN_freq > 0),
    A_sub = mean_percent_sub(A_freq), C_sub = mean_percent_sub(C_freq),
    G_sub = mean_percent_sub(G_freq), T_sub = mean_percent_sub(T_freq),
    del_sub = mean_percent_sub(DEL_freq),
    in_sub = mean_percent_sub(IN_freq)
  ) %>%
  mutate(mut = round(mut / n_samps * 100, 2), ent = round(ent / n_samps * 100, 2)) -> covid_table

left_join(base_df, covid_table, by = c("position" = "POS")) %>%
  select(position, gene, log_ent, mut, ref, A, C, G, T, del, ins) %>%
  write.table(output_file, quote = FALSE, row.names = FALSE)
