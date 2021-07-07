input_file <- snakemake@input$fname_perpos
output_file1 <- snakemake@output$fname_top_regions
output_file2 <- snakemake@output$fname_top_bases
plot_name1 <- snakemake@output$fname_plot_pdf
plot_name2 <- snakemake@output$fname_plot_png

cut_off_quantile <- 0.95 # where do we think a region stops

library(tidyverse)

# read data
covid_POS <- read.table(input_file, header = TRUE)

gene_pos_df <- read.csv(snakemake@input$fname_genes)
# reorder genes
gene_pos_df$gene_name <- factor(gene_pos_df$gene, levels = gene_pos_df$gene)

find_gene <- function(x, gene_pos_df) {
  start_row <- max(c(0, which(gene_pos_df$start <= as.numeric(x))))
  end_row <- min(c(which(gene_pos_df$end >= as.numeric(x)), nrow(gene_pos_df) + 1))
  if (start_row == end_row) {
    # as.character(gene_pos_df$gene_name[start_row])
    start_row
  } else {
    NA
  }
}

covid_POS$gene_no <- sapply(covid_POS$POS, find_gene, gene_pos_df)

# get average entropy per gene
covid_POS %>%
  filter(!is.na(gene_no)) %>%
  mutate(ENT = exp(log_ENT)) %>%
  group_by(gene_no) %>%
  summarise(tot_ent = sum(ENT)) %>%
  mutate(gene_name = gene_pos_df$gene_name[gene_no]) %>%
  full_join(gene_pos_df) %>%
  mutate(
    length = end - start, av_ent = tot_ent / length,
    log_av_ent = log(av_ent)
  ) -> gene_ent_df

# covid_POS %>% pull(log_ENT) %>% quantile(0.95)

### Find top diverse bases and top diverse regions
cut_off <- covid_POS %>%
  pull(log_ENT) %>%
  quantile(cut_off_quantile) # threshold to be considered consecutive
covid_POS %>%
  filter(log_ENT > cut_off) %>%
  pull(POS) -> high_ent_POS

length(high_ent_POS)

consec <- which(high_ent_POS[-1] - high_ent_POS[-length(high_ent_POS)] == 1)
deletions <- sort(union(high_ent_POS[-1][consec], high_ent_POS[-length(high_ent_POS)][consec]))
firsts <- c(deletions[1], deletions[-1][which(deletions[-1] - deletions[-length(deletions)] > 1)])
lasts <- c(deletions[-length(deletions)][which(deletions[-1] - deletions[-length(deletions)] > 1)], deletions[length(deletions)])
del_df <- data.frame(start = firsts, end = lasts, log_ent = NA)

for (ii in 1:nrow(del_df)) {
  del_df[ii, "log_ent"] <- covid_POS %>%
    filter(POS >= del_df[ii, "start"], POS <= del_df[ii, "end"]) %>%
    pull(log_ENT) %>%
    exp() %>%
    mean() %>%
    log() %>%
    round(2)
}

del_df <- del_df[order(del_df$log_ent, decreasing = TRUE), ]

which(del_df$start <= 23403 & del_df$end >= 23403)

colnames(del_df)[3] <- "log_av_ent"
del_df <- del_df[1:20, ]

del_df$gene_no <- sapply(del_df$start, find_gene, gene_pos_df)

del_df %>%
  mutate(gene_name = gene_pos_df$gene_name[gene_no], length = end - start + 1) %>%
  select(-gene_no) -> del_df

# only keep those in a gene
del_df %>% filter(!is.na(gene_name)) -> del_df
del_df <- del_df[1:10, ]

del_df[, c(1, 2, 5, 4, 3)] %>% write.table(output_file1, quote = FALSE, row.names = F)

covid_POS %>%
  filter(log_ENT > cut_off) %>%
  filter(!POS %in% deletions) %>%
  mutate(log_ent = round(log_ENT, 2)) %>%
  select(-log_ENT) -> base_df
base_df <- base_df[order(base_df$log_ent, decreasing = TRUE), ]

which(base_df$POS == 23403)

colnames(base_df)[1] <- "position"
base_df <- base_df[1:20, ]

base_df %>%
  mutate(gene_name = gene_pos_df$gene_name[gene_no]) %>%
  select(-gene_no) -> base_df

# only keep those in a gene
base_df %>% filter(!is.na(gene_name)) -> base_df
base_df <- base_df[1:10, ]

base_df[, c(1, 3, 2)] %>% write.table(output_file2, quote = FALSE, row.names = F)

ggplot() +
  theme_bw() +
  theme(legend.key.size = unit(0.4, "cm")) +
  geom_rect(data = gene_ent_df, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = log_av_ent, colour = gene_name, fill = gene_name)) +
  scale_fill_brewer(palette = "Spectral") +
  scale_colour_brewer(palette = "Spectral") +
  xlab("genomic position") +
  ylab("diversity (log average entropy)") +
  geom_point(data = del_df, mapping = aes(x = start, y = log_av_ent, colour = gene_name)) +
  geom_point(data = del_df, mapping = aes(x = end, y = log_av_ent, colour = gene_name)) +
  geom_segment(data = del_df, mapping = aes(x = start, xend = start, y = 0, yend = log_av_ent, colour = gene_name), alpha = 0.25, linetype = 3) +
  geom_segment(data = del_df, mapping = aes(x = end, xend = end, y = 0, yend = log_av_ent, colour = gene_name), alpha = 0.25, linetype = 3) +
  geom_segment(data = del_df, mapping = aes(x = start, xend = end, y = log_av_ent, yend = log_av_ent, colour = gene_name), linetype = 3) +
  geom_point(data = base_df, mapping = aes(x = position, y = log_ent, colour = gene_name), shape = 6) +
  geom_segment(data = base_df, mapping = aes(x = position, xend = position, y = 0, yend = log_ent, colour = gene_name), alpha = 0.25, linetype = 3)
ggsave(plot_name1, width = 8, height = 2.5)
ggsave(plot_name2, width = 8, height = 2.5)
