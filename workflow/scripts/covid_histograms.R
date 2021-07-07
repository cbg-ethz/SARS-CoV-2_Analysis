input_file1 <- snakemake@input$fname_persample
input_file3 <- snakemake@input$fname_perpos
input_file5 <- snakemake@input$fname_persamplemuts

outdir <- snakemake@output$outdir

library(tidyverse)

### Sample wise

log_ent_samples_public_df <- read.table(input_file1, header = TRUE)

log_ent_samples_public_df$cohort <- "Public"

log_ent_samples_df <- rbind(log_ent_samples_public_df)

log_ent_samples_df %>% ggplot(aes(
  x = log_ENT, y = ..density..,
  colour = cohort, fill = cohort
)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_histogram(position = "identity", bins = 30, show.legend = FALSE, alpha = 0.5) +
  xlab("log total entropy per sample") +
  scale_colour_manual(values = c("darkorange", "firebrick3")) +
  scale_fill_manual(values = c("darkorange", "firebrick3")) -> p1

x_max <- 1250
y_max <- 0.04

log_ent_samples_df %>%
  mutate(ENT = exp(log_ENT)) %>%
  ggplot(aes(x = ENT, y = ..density.., colour = cohort, fill = cohort)) +
  xlim(c(0, x_max)) +
  ylim(c(0, y_max)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_histogram(position = "identity", boundary = 0, bins = 60, alpha = 0.5) +
  xlab("total entropy per sample") +
  scale_colour_manual(values = c("darkorange", "firebrick3")) +
  scale_fill_manual(values = c("darkorange", "firebrick3")) -> p2

p2 + annotation_custom(ggplotGrob(p1), xmin = x_max / 6, xmax = x_max, ymin = y_max / 6, ymax = y_max) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "./histogram_samples.pdf"), width = 6, height = 4)
ggsave(file.path(outdir, "./histogram_samples.png"), width = 6, height = 4)


### Position wise

log_ent_positions_public_df <- read.table(input_file3, header = TRUE)

log_ent_positions_public_df$cohort <- "Public"

log_ent_positions_df <- rbind(log_ent_positions_public_df)

log_ent_positions_df %>% ggplot(aes(
  x = log_ENT, y = ..density..,
  colour = cohort, fill = cohort
)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_histogram(position = "identity", bins = 30, show.legend = FALSE, alpha = 0.5) +
  xlab("log total entropy per position") +
  scale_colour_manual(values = c("dodgerblue", "darkorchid4")) +
  scale_fill_manual(values = c("dodgerblue", "darkorchid4")) -> p1

x_max <- 32500
y_max <- 0.00185

log_ent_positions_df %>%
  mutate(ENT = exp(log_ENT)) %>%
  ggplot(aes(x = ENT, y = ..density.., colour = cohort, fill = cohort)) +
  xlim(c(0, x_max)) +
  ylim(c(0, y_max)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_histogram(position = "identity", boundary = 0, bins = 60, alpha = 0.5) +
  xlab("total entropy per position") +
  scale_colour_manual(values = c("dodgerblue", "darkorchid4")) +
  scale_fill_manual(values = c("dodgerblue", "darkorchid4")) -> p2

p2 + annotation_custom(ggplotGrob(p1), xmin = x_max / 6, xmax = x_max, ymin = y_max / 6, ymax = y_max) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "./histogram_positions.pdf"), width = 6, height = 4)
ggsave(file.path(outdir, "./histogram_positions.png"), width = 6, height = 4)


### Mutations per sample

mut_samples_public_df <- read.table(input_file5, header = TRUE)

mut_samples_public_df$cohort <- "Public"

mut_samples_df <- rbind(mut_samples_public_df)

mut_samples_df %>%
  ggplot(aes(x = log_MUT, y = ..density.., colour = cohort, fill = cohort)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_histogram(position = "identity", bins = 30, show.legend = FALSE, alpha = 0.5) +
  xlab("log number of mutations per sample") +
  scale_colour_manual(values = c("darkgoldenrod", "forestgreen")) +
  scale_fill_manual(values = c("darkgoldenrod", "forestgreen")) -> p1

x_max <- 18000
y_max <- 0.0027

mut_samples_df %>%
  mutate(MUT = exp(log_MUT)) %>%
  ggplot(aes(x = MUT, y = ..density.., colour = cohort, fill = cohort)) +
  xlim(c(0, x_max)) +
  ylim(c(0, y_max)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  geom_histogram(position = "identity", boundary = 0, bins = 60, alpha = 0.5) +
  xlab("number of mutations per sample") +
  scale_colour_manual(values = c("darkgoldenrod", "forestgreen")) +
  scale_fill_manual(values = c("darkgoldenrod", "forestgreen")) -> p2

p2 + annotation_custom(ggplotGrob(p1), xmin = x_max / 6, xmax = x_max, ymin = y_max / 6, ymax = y_max) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "./histogram_samples_mut.pdf"), width = 6, height = 4)
ggsave(file.path(outdir, "./histogram_samples_mut.png"), width = 6, height = 4)
