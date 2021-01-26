
ent_div <- FALSE # flag for diversity or entropy

### Select positions that appear more than k times
# and proportion of x labels (position) to keep with keep_frac

k <- 29 # 4%, c. 1k positions
keep_frac <- 5

# entropy function
Entropy <- function(v) {
  ifelse(v==0, 0, -v*log(v))
}

# diversity function
Diversity <- function(v, n) {
  v*(v*n - 1)/(n - 1)
}

### Process data

library(tidyverse)

file_append <- paste0("_", k)
# read data
covid_table <- read.table(snakemake@input$fname, header = TRUE)

covid_samples <- read.csv(snakemake@input$fname_covariates) %>%
  filter(per_base_read_count_median > snakemake@config$thresholds$median_minium) %>% filter(per_base_read_count_lower_quartile > snakemake@config$thresholds$quartile_range[[1]]) %>% filter(per_base_read_count_upper_quartile < snakemake@config$thresholds$quartile_range[[2]]) %>%
  pull(X)

split_accession <- function(x) {
  str_split(x, "__")[[1]][1]
}

covid_table$SAMPLE <- sapply(covid_table$SAMPLE, split_accession)

# filter out low coverage

covid_table %>% filter(SAMPLE %in% covid_samples) -> covid_table

# filter out plus

find_plus <- function(x) {
  any(str_split(x, "_")[[1]] == "plus")
}

covid_table$PLUS <- sapply(covid_table$SAMPLE, find_plus)

covid_table %>% filter(PLUS == FALSE) -> covid_table

if (ent_div) {
# add entropy
covid_table %>% mutate(ENT = Entropy(A_freq) + Entropy(C_freq) +
                           Entropy(G_freq) + Entropy(T_freq) +
                           Entropy(DEL_freq)) -> covid_table
} else {
# add diversity
covid_table %>% mutate(N = X.ADJUSTED_.READ_COUNT,
                       ENT = 1 - Diversity(A_freq, N) - Diversity(C_freq, N) -
                         Diversity(G_freq, N) - Diversity(T_freq, N) -
                         Diversity(DEL_freq, N)) -> covid_table
}

covid_table %>% group_by(POS) %>%
  summarise(n = n()) %>% filter(n > k) %>%
  pull(POS) -> selected_POS

covid_table %>%
  filter(POS %in% selected_POS) -> covid_table_sel

### check number of samples

covid_table_sel %>% pull(SAMPLE) %>% unique() %>% length()

### Create wide data for the entropy

covid_table_sel %>% select(SAMPLE, POS, ENT) %>%
  spread(POS, ENT, fill = 0) -> covid_ent_wide

### Store as a matrix

covid_ent_mat <- as.matrix(covid_ent_wide[, -1])
row.names(covid_ent_mat) <- covid_ent_wide[, 1]

### Read in gene positions

gene_pos_df <- read.csv(snakemake@input$fname_genes)

find_gene <- function(x, gene_pos_df) {
  start_row <- max(c(0, which(gene_pos_df$start <= as.numeric(x))))
  end_row <- min(c(which(gene_pos_df$end >= as.numeric(x)), nrow(gene_pos_df) + 1))
  if (start_row == end_row) {
    #as.character(gene_pos_df$gene_name[start_row])
    start_row
  } else {
    NA
  }
}

### Make heatmap

library(ComplexHeatmap)
library(circlize)

if (ent_div) {
  col_fun = colorRamp2(c(0, 0.16, 1.6), c("white", "#2094c8", "#2040C8"), space = "RGB")
  measure_name <- "entropy"
  at_scale <- 0.4
} else {
  col_fun = colorRamp2(c(0, 0.08, 0.8), c("white", "#9420c8", "#4020C8"), space = "RGB")
  measure_name <- "diversity"
  at_scale <- 0.2
}
col_fun2 = colorRamp2(1:11, RColorBrewer::brewer.pal(11, "Spectral"), space = "RGB")

local_names <- colnames(covid_ent_mat)

local_genes <- sapply(local_names, find_gene, gene_pos_df)

column_ha <- HeatmapAnnotation(gene = local_genes,
                               annotation_name_side = "left",
                               annotation_name_gp = gpar(fontsize = 10),
                               col = list(gene = col_fun2),
                               annotation_legend_param = list(
                                 gene = list(color_bar = "discrete",
                                             at = 11:1,
                                             labels = rev(as.character(gene_pos_df$gene_name)),
                                             title = "\ngene",
                                             labels_gp = gpar(fontsize = 8),
                                             title_gp = gpar(fontsize = 10))
                               ))


if (keep_frac > 1) {
  local_names[which(1:length(local_names)%% keep_frac != 0)] <- ""
}

colnames(covid_ent_mat) <- local_names

hm_ent <- Heatmap(covid_ent_mat, name = measure_name, col = col_fun,
        row_title = "sample", column_title_side = "bottom",
        show_row_names = FALSE,
        #show_column_names = FALSE,
        #column_labels = local_names,
        column_names_gp = gpar(fontsize = 8),
        column_title = "position", cluster_rows = TRUE, cluster_columns = FALSE,
        bottom_annotation = column_ha,
        heatmap_legend_param = list(title = measure_name,
                                    at = c(0:4)*at_scale,
                                    labels_gp = gpar(fontsize = 10),
                                    title_gp = gpar(fontsize = 10)),
        use_raster = TRUE,
        raster_device = "png"
  )

#ht_opt("TITLE_PADDING" = unit(1.5, "mm"))

if (ent_div) {
png(file.path(snakemake@output$outdir, paste0("./figures/heatmap_entropy_basel", file_append, ".png")), width = 8, height = 4, units = "in", res = 300)
  draw(hm_ent, padding = unit(c(0, 0, 1, 1), "mm"), merge_legend = TRUE)
dev.off()
} else {
  png(file.path(snakemake@output$outdir, paste0("./figures/heatmap_diversity_basel", file_append, ".png")), width = 8, height = 4, units = "in", res = 300)
  draw(hm_ent, padding = unit(c(0, 0, 1, 1), "mm"), merge_legend = TRUE)
  dev.off()
}
#pdf(file.path(snakemake@output$outdir, paste0("./figures/heatmap_entropy", file_append, ".pdf")), width = 8, height = 4)
#  draw(hm_ent, merge_legend = TRUE)
#dev.off()
