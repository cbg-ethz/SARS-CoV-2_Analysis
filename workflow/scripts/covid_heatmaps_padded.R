ent_div <- FALSE # TRUE # flag for diversity or entropy
k <- 4 # which positions in k% of samples were kept
n_samps <- 8000 # how many samples to plot in the heatmaps
# 1000 takes about 50 seconds
# 2000 takes about 100 seconds
# 4000 takes about 235 seconds
# 6000 takes about 540 seconds
# 8000 takes about 1000 seconds
input_file1 <- paste0(snakemake@input$fname)
output_file1 <- paste0(snakemake@output$fname)

set.seed(101)


# Start the clock!
ptm <- proc.time()

library(tidyverse)

covid_ent_public <- read.table(input_file1, header = TRUE)
if (ent_div) {
  covid_ent_public %>%
    select(SAMPLE, POS, ENT) %>%
    spread(POS, ENT, fill = 0) -> covid_ent_mat_public
} else {
  covid_ent_public %>%
    select(SAMPLE, POS, DIV) %>%
    spread(POS, DIV, fill = 0) -> covid_ent_mat_public
}
row.names(covid_ent_mat_public) <- covid_ent_mat_public[, 1]
covid_ent_mat_public <- covid_ent_mat_public[, -1]
public_pos <- colnames(covid_ent_mat_public)

## pad the matrices

if (n_samps < nrow(covid_ent_mat_public)) { # public one
  covid_ent_mat_p <- as.matrix(covid_ent_mat_public[sample(nrow(covid_ent_mat_public), n_samps), ])
} else {
  covid_ent_mat_p <- as.matrix(covid_ent_mat_public)
}

### Read in gene positions

gene_pos_df <- read.csv(snakemake@input$fname_genes)

find_gene <- function(x, gene_pos_df) {
  start_row <- max(c(0, which(gene_pos_df$start <= as.numeric(x))))
  end_row <- min(c(which(gene_pos_df$end >= as.numeric(x)), nrow(gene_pos_df) + 1))
  if (start_row == end_row) {
    # as.character(gene_pos_df$gene[start_row])
    start_row
  } else {
    NA
  }
}

### Make heatmap

library(ComplexHeatmap)
library(circlize)

if (ent_div) {
  col_fun <- colorRamp2(c(0, 0.18, 1.8), c("white", "#2094c8", "#2040C8"), space = "RGB")
  measure_name <- "entropy"
  at_scale <- 0.45
} else {
  col_fun <- colorRamp2(c(0, 0.083, 0.83), c("white", "#9420c8", "#4020C8"), space = "RGB")
  measure_name <- "nucleotide\ndiversity"
  at_scale <- 0.2
}

col_fun2 <- colorRamp2(1:11, RColorBrewer::brewer.pal(11, "Spectral"), space = "RGB")

local_names <- colnames(covid_ent_mat_p)

local_genes <- sapply(local_names, find_gene, gene_pos_df)

column_ha <- HeatmapAnnotation(
  gene = local_genes,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  col = list(gene = col_fun2),
  annotation_legend_param = list(
    gene = list(
      color_bar = "discrete",
      at = 11:1,
      labels = rev(as.character(gene_pos_df$gene)),
      title = "\ngene",
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 10)
    )
  )
)


keep_frac <- 20

if (keep_frac > 1) {
  local_names[which(1:length(local_names) %% keep_frac != 0)] <- ""
}

hm_ent_p <- Heatmap(covid_ent_mat_p,
  name = measure_name, col = col_fun,
  row_title = "sample", column_title_side = "bottom",
  show_row_names = FALSE,
  # show_column_names = FALSE,
  column_labels = local_names,
  column_names_gp = gpar(fontsize = 8),
  column_title = "position", cluster_rows = TRUE, cluster_columns = FALSE,
  bottom_annotation = column_ha,
  heatmap_legend_param = list(
    title = measure_name,
    at = c(0:4) * at_scale,
    labels_gp = gpar(fontsize = 10),
    title_gp = gpar(fontsize = 10)
  ),
  use_raster = TRUE,
  raster_device = "png"
)

ht_opt("TITLE_PADDING" = unit(1.5, "mm"))

png(output_file1, width = 8, height = 3.6, units = "in", res = 600)
draw(hm_ent_p, padding = unit(c(0, 0, 1, 1), "mm"), merge_legend = TRUE)
dev.off()

# Stop the clock
runtime <- (proc.time() - ptm)[1]
print(runtime)
