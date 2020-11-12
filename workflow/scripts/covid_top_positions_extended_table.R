
base_df <- read.csv(snakemake@input$fname_top_bases)

# entropy function
Entropy <- function(v) {
  ifelse(v==0, 0, -v*log(v))
}

mean_percent <- function(x) {
  round(mean(x)*100, 2)
}

### Process data

library(tidyverse)

# read data
covid_table <- read.table(snakemake@input$fname_vcf, header = TRUE)

n_samps <- covid_table %>% pull(SAMPLE) %>% unique() %>% length

covid_table %>% filter(POS %in% c(base_df$position, 23403)) -> covid_table

covid_table %>% mutate(ENT = Entropy(A_freq) + Entropy(C_freq) +
                     Entropy(G_freq) + Entropy(T_freq) +
                     Entropy(DEL_freq)) %>% group_by(POS) %>%
  summarise(mut = n(), ent = sum(ENT > 0), ref = unique(REF_BASE),
            A = mean_percent(A_freq), C = mean_percent(C_freq),
            G = mean_percent(G_freq), T = mean_percent(T_freq),
            del = mean_percent(DEL_freq)) %>%
  mutate(mut = round(mut/n_samps*100, 2), ent = round(ent/n_samps*100, 2)) -> covid_table

covid_table %>% filter(POS == 23403)

left_join(base_df, covid_table, by = c("position" = "POS")) %>%
  write.table(snakemake@output$fname, sep = ",", quote = FALSE, row.names = F)
