# Reassortant Heatmap
library(tidyverse)
library(slider)
# Import existing metadata

reassortant_metadata <- read_csv('./data/metadata/h5_metadata.csv') %>%
  mutate(date = dmy(date) %>%
           as.Date())

reassortant_metadata_formatted <- FormatMetadata(reassortant_metadata) 


lines <- reassortant_metadata_formatted %>%
  mutate(collection.date = ymd(collection.date)) %>%
  mutate(collection.region.name = str_split_i(collection.region.name, pattern = ' ', i = 2)) %>%
  mutate(collection.region.name = factor(collection.region.name, levels = c('europe', 'asia', 'america', 'africa'))) %>%
    filter(!is.na(collection.date)) %>%
    filter(!is.na(collection.region.name)) %>%
  mutate(collection.month = month <- format(collection.date, "%m")) %>%
    select(collection.region.name, collection.date, cluster.profile, collection.month) %>%
    arrange(collection.region.name, collection.date) %>%
    summarise(ratio = n_distinct(cluster.profile)/n(), n_reassortants = n_distinct(cluster.profile), n_sequences = n(), .by = c('collection.date', 'collection.region.name')) %>%
    group_by(collection.region.name) %>%
    mutate(indexed_monthly = slide_index_dbl(
      .x = ratio,                       # calculate on new_cases
      .i = collection.date,       # indexed with date_onset 
      .f = ~ exp(mean(log(.x))),     # function is sum() with missing values removed
      .after = weeks(12))               # window is the DAY and 6 prior DAYS
    ) %>%
    ungroup() 
  
  
ggplot(lines) + 
  geom_col(aes(x = collection.date, y = n_sequences/40)) + 
  geom_line(aes(x = collection.date, y = indexed_monthly, colour= collection.region.name), linewidth = 1)+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*40, name="Number of Sequences"),
                       'Reassortant Ratio',
                       limits = c(0,1.1),
                       breaks = seq(0,1,by = 0.2),
                       expand = c(0,0)) + 
  facet_wrap(.~ collection.region.name)+ 
  scale_color_brewer(palette = 'Set1') + 
  my_theme 



grouped_lines <- reassortant_metadata_formatted %>%
  mutate(collection.date = ymd(collection.date)) %>%
  mutate(collection.region.name = str_split_i(collection.region.name, pattern = ' ', i = 2)) %>%
  mutate(collection.region.name = factor(collection.region.name, levels = c('europe', 'asia', 'america', 'africa'))) %>%
  filter(!is.na(collection.date)) %>%
  filter(!is.na(collection.region.name)) %>%
  mutate(collection.monthyear = format(collection.date, "%Y-%m")) %>%
  select(collection.region.name, collection.date, cluster.profile, collection.monthyear) %>%
  summarise(n_reassortants = n_distinct(cluster.profile), n_sequences = n(), .by = c(collection.monthyear, collection.region.name))

library(brms)
bayes_mod <- brm(n_reassortants ~ n_sequences + collection.region.name,
              family = negbinomial(),
              chains =2,
              iter = 10000,
              data = grouped_lines)

cbind.data.frame(collection.region.name = c('europe', 'asia', 'america', 'africa'), n_sequences = 100) %>%
  as_tibble()%>%
  add_epred_draws(bayes_mod, ndraws = 100) %>%
  ggplot(.,  aes(x = .epred, y = collection.region.name)) +
  stat_halfeye() +
  theme_minimal()


cbind.data.frame(collection.region.name = c('europe', 'asia', 'america', 'africa'), n_sequences = 100) %>%
  as_tibble()%>%
  add_epred_draws(bayes_mod, ndraws = 100) %>%
  mean_hdi()
  

bayes_mod %>%
  emmeans(~ collection.region.name,
          at = list(collection.region.name =  c('europe', 'asia', 'america', 'africa')),
          re_formula = NA) %>%
  contrast(method = 'revpairwise', ref= 'asia', type = 'response') %>%
  gather_emmeans_draws() %>%
  mutate(`.value` = exp(`.value`)) %>% mean_hdi()

  ggplot(., aes(x = `.value`, y = contrast, fill = contrast)) + 
  stat_halfeye(point_interval = mean_hdi,.width = c(.90, .95))
 
reassortant_metadata_formatted %>%
  mutate(collection.date = ymd(collection.date)) %>%
  ggplot(aes(x = collection.date, y = collection.region.name, fill = as.factor(cluster.genome)))+
  geom_density()+
  facet_wrap(.~collection.region.name)