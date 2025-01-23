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
library(cmdstanr)


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
  filter(Observation.date > '2016-05-31' & Observation.date < '2024-05-01') %>%
  drop_na(Observation.date,Region) %>%
  summarise(n_cases = n(), .by = c(Region, Observation.date)) %>%
  arrange(Observation.date) %>%
  rename(collection_monthyear = Observation.date,
         collection_regionname = Region)



numberofsequence_data <- read_csv('./2024-09-09_meta.csv') %>%
  filter(grepl('[hH]5', virus_subtype)) %>%
  rename(collection_monthyear = collection_datemonth) %>%
  filter(collection_date > '2016-05-31') %>%
  
  
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
  arrange(collection_regionname, ym(collection_monthyear)) %>%
  

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
                    ))  %>%
  
  # within each continent...
  arrange(collection_regionname, collection_monthyear) %>%
  group_by(collection_regionname) %>%

  
  # Get the class of n-1 reassortant
  pivot_wider(names_from = group2, values_from = n_reassortants) %>%
  #rowwise() %>%
  
  
  # Note: explicit choice that THE MOST SIGNIFICANT REASSORTANT CLASS is to be included
  mutate(joint_reassortant_class = case_when(
    minor > 0 & is.na(major) & is.na(dominant) ~ 'minor',
    is.na(minor) & major > 0 & is.na(dominant) ~ 'major',
    is.na(minor)  & is.na(major) & dominant > 0 ~ 'dominant',
    
    minor > 0 & major > 0 & is.na(dominant) ~ 'major',
    is.na(minor) & major > 0 & dominant > 0 ~ 'dominant',
    minor > 0 & is.na(major) & dominant > 0 ~ 'dominant',
    
    
    minor > 0 & major > 0 & dominant > 0 ~ 'dominant',
    .default = NA)) %>%
  
  
  fill(joint_reassortant_class, .direction = 'down') %>%
  mutate(previous_reassortant_class = lag(joint_reassortant_class)) %>%
  
  pivot_longer(cols = c(minor, major,  dominant), 
               names_to = 'reassortant_class',
               values_to = 'n_reassortants',
               values_drop_na = TRUE)


count_data <- expand_grid(collection_monthyear = seq.Date(from = ym('2016-06'),
                                                          to =  ym('2024-06'),
                                                          by = "month"),
                          collection_regionname = unique(reassortant_counts$collection_regionname)) %>%
  mutate(collection_monthyear = format.Date(collection_monthyear, "%Y-%m")) %>%
  
  left_join(reassortant_counts) %>%
  left_join(numberofsequence_data) %>%
  left_join(h5_hpai_counts) %>%

 
  separate_wider_delim(collection_monthyear, '-', names = c('collection_year', 'collection_month'),   cols_remove = FALSE) %>%
  mutate(across(c('collection_month', 'collection_year'), .fns = ~ as.double(.x))) %>%
  mutate(time = ym(collection_monthyear) %>% decimal_date()) %>%
  
  # within each continent...
  group_by(collection_regionname) %>%

  
  # And infer the time since the last dominant reassortant
  mutate(last_dominant = if_else(grepl('dominant', reassortant_class), time, NA)) %>%
  fill(last_dominant) %>%
  mutate(time_since_last_dominant = as.numeric(time - last_dominant)) %>%
  
  # Classify tmrca month according to breeding season
  mutate(collection_season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                                       collection_month %in% c(3,4,5)  ~ 'migrating_spring', # Rename to spring migration
                                       collection_month %in% c(6,7,8)  ~ 'breeding', 
                                       collection_month %in% c(9,10,11)  ~ 'migrating_autumn' # Rename to autumn migration
         )) %>%

  # select variables
  ungroup(collection_regionname) %>%
  select(collection_regionname,
         collection_year,
         collection_month,
         collection_season, 
         collection_monthyear,
         n_cases,
         n_sequences,
         n_reassortants,
         reassortant_class,
         joint_reassortant_class,
         previous_reassortant_class,
         time,
         time_since_last_dominant,
         ) %>%

  
  # scaling
  mutate(across(c(time, collection_year), .fns = ~ subtract(.x, 2015))) %>%
  
  # set default NA values
  replace_na(list(n_cases = 0, 
                  n_sequences = 0,
                  n_reassortants = 0,
                  time_since_last_dominant = 0, # this may be quite a strong assumption (accidentally)
                  reassortant_class = 'none',
                  last_reassortant_class = 'none',
                  joint_reassortant_class = 'none',
                  previous_reassortant_class = 'none'
                  )) 

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
bf1 <- 



# to consider: spline/autoregression for month
# number of cases
# collection_regionname as predictor
# Model Formula
bf1 <- bf(n_reassortants ~ collection_regionname + n_cases  +  collection_season + previous_reassortant_class + (1 + collection_regionname + collection_month|collection_year),
           zi ~ collection_regionname + n_cases  + (1 + collection_regionname + collection_month|collection_year),
           family = zero_inflated_negbinomial())


countmodel_bf1_priors <- get_prior(bf1, count_data)

countmodel_bf1_priors[c(1:9,23:27),1] <- "normal(0,5)"
countmodel_bf1_priors[c(1:7,20:24),1] <- "normal(0,5)"
#countmodel_bf1_priors[c(1:7,22:26),1] <- "normal(0,5)"
countmodel_bf1_priors
  
countmodel_bf1_fit <- brm(
  bf1,
  data = count_data ,           
  prior = countmodel_bf1_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97,
                 max)
)

### Temporary model testing and visulaisation ###
pred_model_continent_only <- countmodel_bf1_fit %>%
  epred_draws(newdata = count_data, re_formula = NULL, ndraws = 100) 

pred_model_continent_only %>%
  ggplot(aes(x = ym(collection_monthyear), y = .epred)) +
  stat_lineribbon(alpha = 0.5) +
  geom_point(data = count_data , aes(y = n_reassortants, x = ym(collection_monthyear)),
             color = "grey50", size = 3, alpha = 0.5) +
  
  scale_fill_brewer(palette = "Reds") +
  facet_grid(rows = vars(collection_regionname))



#no_zero_counts <- count_data %>% filter(n_reassortants > 0) %>% filter(previous_reassortant_class != 'none') 

bf2 <- bf(reassortant_class ~  collection_regionname + previous_reassortant_class + (1 + collection_regionname + collection_month|collection_year), 
          family = categorical(link  = 'logit'))

countmodel_bf2_priors <- get_prior(bf2, count_data)

countmodel_bf2_priors[c(3:9,18:24, 32:38),1] <- "normal(0,5)"

countmodel_bf2_fit<- brm(
  bf2,
  data =  count_data,           
  prior = countmodel_bf2_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = 8000,
  warmup = 800,
  seed = SEED,
  control = list(adapt_delta = 0.97)
)



#newdata <- expand_grid(collection_regionname = unique(count_data %>% filter(n_reassortants > 0) %>% pull(collection_regionname)),
                    #   collection_monthyear = unique(count_data %>% filter(n_reassortants > 0) %>% pull(collection_monthyear)))  %>%
  #separate_wider_delim(collection_monthyear, '-', names = c('collection_year', 'collection_month'),   cols_remove = FALSE)  %>%
 # mutate(across(c('collection_month', 'collection_year'), .fns = ~ as.double(.x))) %>%
  #mutate(collection_year = collection_year - 2015) %>%
  #mutate(collection_monthyear = ym(collection_monthyear)) 

mv <- mvbf(bf1, bf2)

countmodel_mv_priors <- get_prior(mv, count_data)

#countmodel_mv_priors[c(1, 4:15,25:29, 38:44, 53:59, 69:74),1] <- "normal(0,5)"
countmodel_mv_priors
countmodel_mv_priors[c(1, 4:15, 25:29, 38:44, 53:59, 69:74),1] <- "normal(0,5)"


countmodel_bf2_fit<- brm(
  mv,
  data =  count_data,           
  prior = countmodel_mv_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = 16000,
  warmup = 1600,
  seed = SEED,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97))










pred_model_continent_only <- countmodel_fit_mv %>%
  epred_draws(newdata = count_data, re_formula = NULL, ndraws = 100) 

pred_model_continent_only %>%
ggplot(aes(x = ym(collection_monthyear), y = .epred)) +
  stat_lineribbon(alpha = 0.5) +
  geom_point(data = count_data , aes(y = n_reassortants, x = ym(collection_monthyear)),
             color = "grey50", size = 3, alpha = 0.5) +
 
  scale_fill_brewer(palette = "Reds") +
  facet_grid(rows = vars(collection_regionname))

#count_data %>%
  #add_residual_draws(countmodel_fit_mv) %>%
  #median_qi %>%
  #ggplot(aes(sample = .residual)) +
  #geom_qq() + 
  #geom_qq_line()

# Post-fitting diagnostics (including inspection of ESS, Rhat and posterior predictive)
tidy_countmodel_fit <- tidy(countmodel_fit)
#posteriorpredictive <-


# Misc evaluations
#print(prior_summary(countmodel_fit), show_df = FALSE)
performance(countmodel_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(countmodel_fit) # default output plot of brms showing posterior distributions of
prior_summary(countmodel_fit) #obtain dataframe of priors used in model.
