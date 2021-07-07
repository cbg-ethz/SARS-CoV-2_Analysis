input_file1 <- snakemake@input$fname_entropy
input_file2 <- snakemake@input$fname
output_file1 <- snakemake@output$fname
output_file2 <- snakemake@output$fname_hist
output_file3 <- snakemake@output$fname_regression

library(tidyverse)

# read data
covid_ENT <- read.table(input_file1, header = TRUE)

covid_covariates <- read.csv(input_file2)
covid_covariates %>%
  mutate(SAMPLE = X) %>%
  inner_join(covid_ENT) -> covid_covs
# select data with sex information
covid_covs %>%
  pull(host_sex) %>%
  table()
covid_covs %>% filter(host_sex %in% c("male", "female")) -> covid_cov_sel

# select adult data with age information
covid_cov_sel %>%
  pull(host_age) %>%
  table()
covid_cov_sel %>%
  filter(host_age %in% as.character(18:100)) %>%
  mutate(host_age = as.numeric(as.character(host_age))) -> covid_cov_sel2

# check collection dates
covid_cov_sel2 %>%
  pull(collection_date) %>%
  table()
# pick the exact dates
covid_cov_sel2 %>% filter(!collection_date %in% c("2020", "missing")) -> covid_cov_sel3


covid_cov_sel3 %>%
  pull(country) %>%
  table() # empty
covid_cov_sel3 %>%
  pull(geographic.location..region.and.locality.) %>%
  table() # empty
covid_cov_sel3 %>%
  pull(geographic.location..country.and.or.sea.) %>%
  table() # empty
covid_cov_sel3 %>%
  pull(geo_loc_name) %>%
  table()

covid_cov_sel3$country <- sapply(covid_cov_sel3$geo_loc_name, function(x) strsplit(x, ":")[[1]][1])
covid_cov_sel3 %>%
  pull(country) %>%
  table() # now full

# pick countries with more than 100 samples, so Australia and USA
covid_cov_sel3 %>% filter(country %in% c("Australia", "USA")) -> covid_cov_sel4

# select studies with reasonable numbers
covid_cov_sel4 %>%
  select(study_accession) %>%
  table()
covid_cov_sel4 %>% filter(study_accession %in% c("SRP253798")) -> covid_cov_sel5

# select sequencing type with reasonable numbers
covid_cov_sel5 %>%
  pull(library_selection) %>%
  table()
# all PCR

# check instrument types
covid_cov_sel5 %>%
  pull(instrument) %>%
  table()
# so if we include instrument, we have county and the other ones too

covid_cov_sel5 %>% mutate(col_date = as.Date(collection_date) - as.Date("2020-01-01")) -> covid_cov_sel6
covid_cov_sel6 %>%
  pull(col_date) %>%
  table()


covid_cov_sel6 %>%
  mutate(
    log_count = log(per_base_read_count_median),
    log_IQR = log(per_base_read_count_upper_quartile - per_base_read_count_lower_quartile) - log_count
  ) %>%
  select(
    SAMPLE, log_count, log_IQR, log_ENT,
    host_sex, host_age, instrument, col_date
  ) -> covid_df

# write to file
write.table(covid_df, file = output_file1, row.names = FALSE, quote = FALSE)

# look at the time distribution

cut_off_date <- 135 # mid may

covid_df %>%
  mutate(early_late = col_date > cut_off_date) %>%
  ggplot(aes(x = col_date, y = ..count.., color = early_late, fill = early_late)) +
  theme_bw() +
  theme(text = element_text(size = 16)) + # xlim(c(15, 105)) +
  geom_histogram(boundary = 0, bins = 20, alpha = 0.5) +
  scale_fill_manual(values = c("slategray3", "slateblue3"), labels = c("early", "later")) +
  scale_color_manual(values = c("slategray3", "slateblue3"), labels = c("early", "later")) +
  xlab("time (days)") +
  ylab("cases") +
  labs(fill = "", color = "")
ggsave(output_file2, width = 7, height = 3.5)

covid_lm_base <- lm(log_ENT ~ log_count + log_IQR +
  instrument, covid_df)

### adjusted entropy

lm_coeffs_base <- (summary(covid_lm_base)$coefficients)[, "Estimate"]

covid_df %>% mutate(log_ENT_adj = log_ENT - log_count * lm_coeffs_base["log_count"] -
  log_IQR * lm_coeffs_base["log_IQR"] -
  (instrument == "NextSeq 500") * lm_coeffs_base["instrumentNextSeq 500"] -
  (instrument == "NextSeq 550") * lm_coeffs_base["instrumentNextSeq 550"]) -> covid_df_adj_base

covid_df_adj_base %>%
  mutate(early_late = col_date > cut_off_date) %>%
  ggplot(aes(x = col_date, y = log_ENT_adj, color = early_late)) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  xlim(c(20, 300)) +
  geom_point(alpha = 0.3, shape = 16) +
  xlab("time (days)") +
  ylab("log total entropy (adjusted)") +
  # add regression line
  geom_smooth(method = lm) +
  labs(color = "") +
  scale_color_manual(values = c("slategray3", "slateblue3"), labels = c("early", "later"))
ggsave(output_file3, width = 7, height = 3.5)


# regress on remaining covariates
# (note sequencing, country etc only have one level, so are not included)

# split into early and later time points

covid_df_early <- filter(covid_df, col_date < cut_off_date)

covid_lm_early <- lm(log_ENT ~ host_sex + host_age + col_date + log_count + log_IQR +
  instrument, covid_df_early)

summary(covid_lm_early)

covid_df_later <- filter(covid_df, col_date >= cut_off_date)

covid_lm_later <- lm(log_ENT ~ host_sex + host_age + col_date + log_count + log_IQR +
  instrument, covid_df_later)

summary(covid_lm_later)

covid_df %>%
  pull(log_count) %>%
  sd()

covid_df %>%
  pull(log_IQR) %>%
  sd()

covid_df %>%
  pull(host_sex) %>%
  table()
covid_df %>%
  pull(host_sex) %>%
  table() %>%
  {
    .
  } / nrow(covid_df)

covid_df %>%
  pull(host_age) %>%
  summary()

covid_df %>%
  pull(instrument) %>%
  table()
covid_df %>%
  pull(instrument) %>%
  table() %>%
  {
    .
  } / nrow(covid_df)
