####################################################################################################
####################################################################################################
## Script name: Number and Class of Reassortants Model
##
## Purpose of script: to determine whether reassortants emerging from some regions are more
## likely and more severe than those from others. Variables of interest interest include: region of
## origin, persistence time (persistence time increases uncertainty), number of sequences per region
##
## This script implements a workflow for a multivariate negative binomial model. 
##
## Date created: 2024-xx-xx
##
##
########################################## SYSTEM OPTIONS ##########################################
#options(scipen = 6, digits = 7) 
memory.limit(30000000) 


########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(magrittr)
library(ggmcmc)


# User functions


############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')


############################################## MAIN ################################################
# Data pre-processing
h5_hpai_counts <- read_csv('./2025Jan17_faoHPAI.csv') %>%
  filter(grepl('Avian', Disease) & grepl('H5N', Serotype)) %>%
  mutate(Region = case_when(grepl('Europe', 
                                  Region ) ~ 'europe',
                            grepl('Africa',
                                  Region ) ~ 'africa',
                            grepl('Asia', 
                                  Region ) ~ 'asia',
                            grepl('America', 
                                  Region ) ~ 'central & northern america',
                            .default = NA_character_ )) %>%
  mutate(Observation.date = dmy(Observation.date) %>% format.Date('%Y-%m')) %>%
  filter(Observation.date > '2016-07-31') %>%
  drop_na(Observation.date,Region) %>%
  summarise(n_cases = n(), .by = c(Region, Observation.date)) %>%
  arrange(Observation.date) %>%
  rename(collection_monthyear = Observation.date,
         collection_regionname = Region)



numberofsequence_data <- read_csv('./2024-09-09_meta.csv') %>%
  filter(grepl('[hH]5', virus_subtype)) %>%
  rename(collection_monthyear = collection_datemonth) %>%
  filter(collection_date > '2016-07-31') %>%
  
  select(c(collection_monthyear,
           collection_regionname,
           cluster_profile)) %>%
  drop_na(collection_monthyear) %>%
  mutate(collection_regionname = case_when(grepl('europe', 
                                                 collection_regionname) ~ 'europe',
                                           grepl('africa',
                                                 collection_regionname) ~ 'africa',
                                           grepl('asia', 
                                                 collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', 
                                                 collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname)) %>%

  #mutate(collection_year = 
  #paste0(collection_datemonth, '-01') %>%
  # ymd() %>%
  #decimal_date() %>%
  #round(digits = 1)) %>%
  summarise(n_sequences = n(), .by = c(collection_monthyear, collection_regionname))




# Stratified by region (for supplementary)
counts %>%
  filter(collection_regionname %in% c('europe', 'asia','central & northern america','africa')) %>%
  pivot_longer(starts_with('n_'), names_to = 'var', values_to = 'value') %>%
  mutate(collection_monthyear = ym(collection_monthyear)) %>%
  filter(collection_monthyear > '2017-01-01') %>%
  ggplot(aes(x=collection_monthyear, y= value, colour = var)) + 
  geom_bar(aes(x = ym(collection_monthyear), y = n_reassortants*30), stat = 'identity', data = count_data %>% filter(collection_monthyear > '2017-01-01'), inherit.aes = FALSE) + 
  geom_point(size =1) +
  geom_line(aes(group = var)) +
  facet_grid(rows = vars(collection_regionname)) + 
  scale_x_date()  +
  scale_y_continuous(sec.axis = sec_axis( trans=~./30, name="Reassortment")) +
  global_theme

# Global (for supplementary)
counts %>%
  filter(collection_regionname %in% c('europe', 'asia','central & northern america','africa')) %>%
  pivot_longer(starts_with('n_'), names_to = 'var', values_to = 'value') %>%
  summarise(value = sum(value), .by =  c(collection_monthyear, var)) %>%
  mutate(collection_monthyear = ym(collection_monthyear)) %>%
  filter(collection_monthyear > '2017-01-01') %>%
  ggplot(aes(x=collection_monthyear, y= value, colour = var)) + 
  geom_bar(aes(x = ym(collection_monthyear), y = n_reassortants*30), stat = 'identity', data = count_data %>% filter(collection_monthyear > '2017-01-01'), inherit.aes = FALSE) + 
  geom_point(size =1) +
  geom_line(aes(group = var)) +
  scale_x_date()  +
  scale_y_continuous(sec.axis = sec_axis( trans=~./30, name="Reassortment"))+
  global_theme


reassortant_counts <- combined_data %>%
  filter(segment == 'ha') %>%
  mutate(collection_monthyear = date_decimal(TMRCA) %>% 
           format(., "%Y-%m")) %>%
  mutate(collection_regionname = case_when(grepl('europe', 
                                                 collection_regionname) ~ 'europe',
                                           grepl('africa',
                                                 collection_regionname) ~ 'africa',
                                           grepl('asia', 
                                                 collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', 
                                                 collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname)) %>%
  filter(!is.na(collection_regionname)) %>%
  select(collection_regionname,
         collection_monthyear, 
         cluster_profile,
         group2) %>% # include all month-years over collection period to generate zer countr
  summarise(n_reassortants = n_distinct(cluster_profile), 
            .by = c(collection_monthyear,
                    collection_regionname,
                    group2
                    ))  #%>%
  #mutate(collection_year = collection_monthyear %>%
  #  paste0('-01') %>%
  # ymd() %>%
  #decimal_date() %>%
  #round(., digits = 1))  %>%
  #summarise(n_reassortants = sum(n_reassortants), 
          #  .by = c(collection_monthyear,collection_regionname)) %>%
  
  #arrange(collection_year) %>%
  
count_data <- h5_hpai_counts %>%
  left_join(numberofsequence_data) %>%
  left_join(reassortant_counts) %>%
  replace_na(list(n_case = 0, 
                  n_sequences = 0,
                  n_reassortants = 0,
                  group2 = 'none')) %>%
 
  separate_wider_delim(collection_monthyear, '-', names = c('collection_year', 'collection_month'),   cols_remove = FALSE) %>%
  mutate(across(c('collection_month', 'collection_year'), .fns = ~ as.double(.x))) %>%
  mutate(time = ym(collection_monthyear) %>% decimal_date() %>% subtract(2015)) %>%
  mutate(collection_year = collection_year - 2015)%>%
  filter(n_reassortants > 0)
 # mutate(n_cases = scale(n_cases))  

# Model Formula
#bf1 <- 
#bf2 <- bf(group2 ~ 1 + collection_regionname + (1|collection_year), 
         # family = categorical())

#count_formula <- mvbf(bf1, bf2)


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 8000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472http://127.0.0.1:37603/graphics/plot_zoom_png?width=2281&height=1324
THREADS <- 2


# Fit model to data
# count model only
bf1 <- bf(n_reassortants ~ 1 + collection_regionname +s(collection_year, bs = 'cr') + s(collection_month, bs = 'cr') + (1 |collection_year),
          #zi ~ collection_regionname + s(collection_year, bs = 'cr') + s(collection_month, bs = 'cr'),
          family = zero_inflated_negbinomial())


countmodel_fit_mv <- brm(
  count_formula,
  data = bf1 ,           # Dataframe containing y1, y2, x1, x2http://127.0.0.1:37603/graphics/plot_zoom_png?width=2534&height=934
  # prior = countmodel_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.97)
)



# to consider: spline/autoregression for month
# number of cases
# collection_regionname as predictor
# Model Formula
bf1 <- bf(n_reassortants ~ 1 + collection_regionname + (collection_month+collection_regionname |collection_year),
          #zi ~ collection_regionname + s(collection_year, bs = 'cr') + s(collection_month, bs = 'cr'),
          family = zero_inflated_negbinomial())

bf2 <- bf(group2 ~ 1 + collection_regionname , 
          family = categorical)

count_formula <- mvbf(bf1, bf2)

countmodel_fit_mv <- brm(
  count_formula,
  data = count_data ,           # Dataframe containing y1, y2, x1, x2http://127.0.0.1:37603/graphics/plot_zoom_png?width=2534&height=934
 # prior = countmodel_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.97)
)

#newdata <- expand_grid(collection_regionname = unique(count_data %>% filter(n_reassortants > 0) %>% pull(collection_regionname)),
                    #   collection_monthyear = unique(count_data %>% filter(n_reassortants > 0) %>% pull(collection_monthyear)))  %>%
  #separate_wider_delim(collection_monthyear, '-', names = c('collection_year', 'collection_month'),   cols_remove = FALSE)  %>%
 # mutate(across(c('collection_month', 'collection_year'), .fns = ~ as.double(.x))) %>%
  #mutate(collection_year = collection_year - 2015) %>%
  #mutate(collection_monthyear = ym(collection_monthyear)) 
  

pred_model_continent_only <- countmodel_fit_mv %>%
  epred_draws(newdata = count_data, re_formula = NULL, ndraws = 100) 

pred_model_continent_only %>%
ggplot(aes(x = ym(collection_monthyear), y = .epred)) +
  stat_lineribbon(alpha = 0.5) +
  geom_point(data = count_data , aes(y = n_reassortants, x = ym(collection_monthyear)),
             color = "grey50", size = 3, alpha = 0.5) +
 
  scale_fill_brewer(palette = "Reds") +
  facet_grid(rows = vars(collection_regionname))


# Post-fitting diagnostics (including inspection of ESS, Rhat and posterior predictive)
tidy_countmodel_fit <- tidy(countmodel_fit)
#posteriorpredictive <-pp_check(countmodel_fit, ndraws = 500)


# Misc evaluations
#print(prior_summary(countmodel_fit), show_df = FALSE)
performance(countmodel_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(countmodel_fit) # default output plot of brms showing posterior distributions of
prior_summary(countmodel_fit) #obtain dataframe of priors used in model.
