
# entropy function
Entropy <- function(v) {
  ifelse(v==0, 0, -v*log(v))
}

### Process data

library(tidyverse)

# read data
covid_table <- read.table(snakemake@input$fname, header = TRUE)

outdir <- snakemake@output$outdir
dir.create(outdir, recursive = TRUE)

# add entropy
covid_table %>% mutate(ENT = Entropy(A_freq) + Entropy(C_freq) +
                           Entropy(G_freq) + Entropy(T_freq) +
                           Entropy(DEL_freq)) -> covid_table

### check number of samples

covid_table %>% pull(SAMPLE) %>% unique() %>% length()

### Create wide data for the entropy

covid_table %>% select(SAMPLE, POS, ENT) %>%
  spread(POS, ENT, fill = 0) -> covid_ent_wide

### Store as a matrix

covid_ent_mat <- as.matrix(covid_ent_wide[, -1])
row.names(covid_ent_mat) <- covid_ent_wide[, 1]

### Make histograms

rowSums(covid_ent_mat) %>% log %>% data.frame(LOG_ENT = .) %>%
  rownames_to_column("ID") -> log_ent_samples_df

log_ent_samples_df %>% write.csv(snakemake@output$fname_entropy_samples, row.names = FALSE)

log_ent_samples_df %>% ggplot(aes(x = LOG_ENT, y = ..density..)) +
    theme_bw() + theme(text = element_text(size = 16)) +
    geom_histogram(col = "darkorange", bins = 30, fill = "darkorange", alpha = 0.5) +
    xlab("log total entropy per sample")
  ggsave(file.path(outdir, "log_histogram_samples.pdf"), width=7, height=3.5)
  ggsave(file.path(outdir, "log_histogram_samples.png"), width=7, height=3.5)

log_ent_samples_df %>% mutate(ENT = exp(LOG_ENT)) %>%
  ggplot(aes(x = ENT, y = ..density..)) +
    theme_bw() + theme(text = element_text(size = 16)) +
    geom_histogram(col = "darkorange", boundary = 0, bins = 60, fill = "darkorange", alpha = 0.5) +
    xlab("total entropy per sample")
  ggsave(file.path(outdir, "histogram_samples.pdf"), width=7, height=3.5)
  ggsave(file.path(outdir, "histogram_samples.png"), width=7, height=3.5)

colSums(covid_ent_mat) %>% log %>% data.frame(LOG_ENT = .) %>%
  rownames_to_column("POS") -> log_ent_positions_df

log_ent_positions_df %>% write.csv(snakemake@output$fname_entropy_positions, row.names = FALSE)

log_ent_positions_df %>% ggplot(aes(x = LOG_ENT, y = ..density..)) +
    theme_bw() + theme(text = element_text(size = 16)) +
    geom_histogram(col = "firebrick3", bins = 30, fill = "firebrick3", alpha = 0.5) +
    xlab("log total entropy per position")
  ggsave(file.path(outdir, "log_histogram_positions.pdf"), width=7, height=3.5)
  ggsave(file.path(outdir, "log_histogram_positions.png"), width=7, height=3.5)

log_ent_positions_df %>% mutate(ENT = exp(LOG_ENT)) %>%
  ggplot(aes(x = ENT, y = ..density..)) +
    theme_bw() + theme(text = element_text(size = 16)) +
    geom_histogram(col = "firebrick3", boundary = 0, bins = 60, fill = "firebrick3", alpha = 0.5) +
    xlab("total entropy per position")
  ggsave(file.path(outdir, "histogram_positions.pdf"), width=7, height=3.5)
  ggsave(file.path(outdir, "histogram_positions.png"), width=7, height=3.5)
