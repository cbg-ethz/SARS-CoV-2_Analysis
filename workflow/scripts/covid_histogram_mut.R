
### Process data

library(tidyverse)

# read data
covid_table <- read.table(snakemake@input$fname_vcf, header = TRUE)

outdir <- snakemake@output$outdir
dir.create(outdir, recursive = TRUE)

# add entropy
covid_table %>% mutate(MUT = 1) -> covid_table

### check number of samples

covid_table %>% pull(SAMPLE) %>% unique() %>% length()

### Create wide data for the entropy

covid_table %>% select(SAMPLE, POS, MUT) %>%
  spread(POS, MUT, fill = 0) -> covid_mut_wide

### Store as a matrix

covid_mut_mat <- as.matrix(covid_mut_wide[, -1])
row.names(covid_mut_mat) <- covid_mut_wide[, 1]

### Make histograms

rowSums(covid_mut_mat) %>% data.frame(MUT = .) %>%
  rownames_to_column("ID") -> mut_samples_df

mut_samples_df %>% write.csv(snakemake@output$fname, row.names = FALSE)

mut_samples_df$MUT %>% summary()
n_pos <- 29903
mut_samples_df$MUT %>% summary() %>% {.}/n_pos %>% {.}*100

mut_samples_df %>% mutate(LOG_MUT = log(MUT)) %>%
  ggplot(aes(x = LOG_MUT, y = ..density..)) +
    theme_bw() + theme(text = element_text(size = 16)) +
    geom_histogram(col = "darkorchid4", bins = 30, fill = "darkorchid4", alpha = 0.5) +
    xlab("log number of mutations per sample")
  ggsave(file.path(outdir, "log_mut_histogram.pdf"), width=7, height=3.5)
  ggsave(file.path(outdir, "log_mut_histogram.png"), width=7, height=3.5)

mut_samples_df %>%
  ggplot(aes(x = MUT, y = ..density..)) +
    theme_bw() + theme(text = element_text(size = 16)) +
    geom_histogram(col = "darkorchid4", boundary = 0, bins = 60, fill = "darkorchid4", alpha = 0.5) +
    xlab("number of mutations per sample")
  ggsave(file.path(outdir, "mut_histogram.pdf"), width=7, height=3.5)
  ggsave(file.path(outdir, "mut_histogram.png"), width=7, height=3.5)
