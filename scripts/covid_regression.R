
library(tidyverse)

# read data
covid_ENT <- read.csv(snakemake@input$entropy_samples_fname)
# rename columns
colnames(covid_ENT) <- c("accession", "log_ent")

covid_covariates <- read.csv(snakemake@input$covariates_fname)

covid_covariates %>% unite(country, c(geo_loc_name_country, geographic_location_.country_and.or_sea.), sep = "") %>%
  mutate(log_count = log(per_base_read_count_median),
         log_IQR = log(per_base_read_count_upper_quartile - per_base_read_count_lower_quartile) - log_count) %>%
  select(accession, country, host_sex, host_age, sample_type,
         Assay.Type, LibrarySelection, log_count, log_IQR) -> covid_cov_sel

# select data with sex information
covid_cov_sel %>% pull(host_sex) %>% table
covid_cov_sel %>% filter(host_sex %in% c("male", "female")) -> covid_cov_sel

# select data with age information
covid_cov_sel %>% pull(host_age) %>% table
covid_cov_sel %>% filter(host_age != "missing", host_age != "") %>%
  mutate(host_age = as.numeric(as.character(recode(host_age, "4 years" = "4")))) -> covid_cov_sel

# select countries with reasonable numbers
covid_cov_sel %>% pull(country) %>% table
covid_cov_sel %>% filter(country == "Australia") -> covid_cov_sel

# select sequencing type with reasonable numbers
covid_cov_sel %>% pull(sample_type) %>% table

# select sequencing type with reasonable numbers
covid_cov_sel %>% pull(Assay.Type) %>% table

# select sequencing type with reasonable numbers
covid_cov_sel %>% pull(LibrarySelection) %>% table

# add entropy information to dataframe
covid_cov_sel %>% left_join(covid_ENT) -> covid_df

# regress on remaining covariates
# (note sequencing, country etc only have one level, so are not included)
# covid_lm <- lm(log_ent ~ host_sex + host_age + log_count, covid_df)
# add relative IQR logged

covid_lm <- lm(log_ent ~ host_sex + host_age + log_count + log_IQR, covid_df)

summary(covid_lm)

#covid_lm %>% ggplot(aes(x = .fitted, y = .resid)) +
#  geom_point(aes(colour = log_count), size = 2)

sink(snakemake@output$lm_summary_fname)
summary(covid_lm)
sink()
