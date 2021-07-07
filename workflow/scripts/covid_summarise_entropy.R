# scale_factor <- 1 # for public no scaling of total entropy
scale_factor <- 95719 / 18468 # adjusted sizes - should really be read from data
input_file <- snakemake@input@fname_entropy
output_file1 <- snakemake@output$fname_persample
output_file2 <- snakemake@output$fname_perpos
output_file3 <- snakemake@output$fname_persamplemuts
k <- 4 # which positions in k% of samples to keep for later heatmaps
output_file4 <- snakemake@output$fname_filtered

# read data
ent_table <- read.table(input_file, header = TRUE)

library(tidyverse)

# compute number of samples
ent_table %>%
  pull(SAMPLE) %>%
  unique() %>%
  length() -> n_samples

# summarise entropy per sample
ent_table %>%
  group_by(SAMPLE) %>%
  summarise(log_ENT = log(sum(ENT))) %>%
  filter(log_ENT > -Inf) -> ent_table_SAMPLE

# write to file
write.table(ent_table_SAMPLE, file = output_file1, row.names = FALSE, quote = FALSE)

# summarise entropy per position
ent_table %>%
  group_by(POS) %>%
  summarise(log_ENT = log(sum(ENT) * scale_factor)) %>%
  filter(log_ENT > -Inf) -> ent_table_POS

# write to file
write.table(ent_table_POS, file = output_file2, row.names = FALSE, quote = FALSE)

# count number of mutations per sample
ent_table %>%
  group_by(SAMPLE) %>%
  summarise(log_MUT = log(n())) %>%
  filter(log_MUT > -Inf) -> ent_table_SAMPLE_mut

# write to file
write.table(ent_table_SAMPLE_mut, file = output_file3, row.names = FALSE, quote = FALSE)

# find positions affected in k% of samples
ent_table %>%
  group_by(POS) %>%
  summarise(frac = n() / n_samples) %>%
  filter(frac > k / 100) %>%
  pull(POS) -> filtered_POS

# filter data to those positions
ent_table %>% filter(POS %in% filtered_POS) -> ent_table_filtered

# write to file
write.table(ent_table_filtered, file = output_file4, row.names = FALSE, quote = FALSE)
