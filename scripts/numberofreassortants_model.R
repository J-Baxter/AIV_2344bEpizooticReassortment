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

FormatContinent <- function(dataframe){
  dataframe %<>%
    mutate(collection_regionname = case_when(grepl('europe', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'europe',
                                             
                                             grepl('africa',
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'africa',
                                             
                                             grepl('asia', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'asia',
                                             
                                             grepl('(central|northern) america', 
                                                   collection_regionname, 
                                                   ignore.case = TRUE) ~ 'central & northern america',
                                             
                                             .default = NA_character_ ))
  
  return(dataframe)
}

############################################## DATA ################################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')


############################################## MAIN ################################################
# Data pre-processing

# Total number of reported H5 Highly Pathogenic Avian Influenza
# (Source: https://empres-i.apps.fao.org/diseases)

h5_hpai_counts <- read_csv('./2025Jan17_faoHPAI.csv') %>%
  filter(grepl('Avian', Disease) & grepl('H5N', Serotype)) %>%
  rename(collection_monthyear = Observation.date,
         collection_regionname = Region) %>%
  FormatContinent() %>%
  mutate(collection_monthyear = dmy(collection_monthyear) %>% format.Date('%Y-%m')) %>%
  filter(collection_monthyear > '2016-05-31' & collection_monthyear < '2024-05-01') %>%
  drop_na(collection_monthyear, collection_regionname) %>%
  summarise(n_cases = n(), .by = c(collection_regionname, collection_monthyear)) %>%
  arrange(collection_monthyear)


# Number of sequences 
numberofsequence_data <- read_csv('./2024-09-09_meta.csv') %>%
  filter(grepl('[hH]5', virus_subtype)) %>%
  rename(collection_monthyear = collection_datemonth) %>%
  filter(collection_date > '2016-05-31') %>%
  select(c(collection_monthyear,
           collection_regionname,
           cluster_profile)) %>%
  drop_na(collection_monthyear) %>%
  FormatContinent() %>%
  arrange(collection_regionname, ym(collection_monthyear)) %>%
  summarise(n_sequences = n(), 
            .by = c(collection_monthyear, collection_regionname))

# Number and TMRCA of reassortants as inferred from our phylodynamic analysis
reassortant_counts <- combined_data %>%
  
  # restrict to HA for now
  filter(segment == 'ha') %>%
  mutate(collection_monthyear = date_decimal(TMRCA) %>% 
           format(., "%Y-%m")) %>%

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
  
  
  # Get the class of t-1 reassortant
  pivot_wider(names_from = group2, values_from = n_reassortants) %>%

  mutate(joint_reassortant_class = case_when(
    minor > 0 & is.na(major) & is.na(dominant) ~ 'minor',
    is.na(minor) & major > 0 & is.na(dominant) ~ 'major',
    is.na(minor)  & is.na(major) & dominant > 0 ~ 'dominant',
    minor > 0 & major > 0 & is.na(dominant) ~ 'minor_major',
    is.na(minor) & major > 0 & dominant > 0 ~ 'major_dominant',
    minor > 0 & is.na(major) & dominant > 0 ~ 'minor_dominant',
    minor > 0 & major > 0 & dominant > 0 ~ 'minor_major_dominant',
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



# Model Formula
bf1 <- bf(n_reassortants ~ collection_regionname + n_cases  +  collection_season + previous_reassortant_class + (1 + collection_regionname + collection_month|collection_year),
          zi ~ collection_regionname + n_cases  + (1 + collection_regionname + collection_month|collection_year),
          family = zero_inflated_negbinomial())



bf2 <- bf(reassortant_class ~  collection_regionname + previous_reassortant_class + (1 + collection_regionname + collection_month|collection_year), 
          family = categorical(link  = 'logit'))
count_formula <- mvbf(bf1, bf2)


# Define Priors
mv <- mvbf(bf1, bf2)

countmodel_mv_priors <- get_prior(mv, count_data)

#countmodel_mv_priors[c(1, 4:15,25:29, 38:44, 53:59, 69:74),1] <- "normal(0,5)"
countmodel_mv_priors
countmodel_mv_priors[c(1, 4:15, 25:29, 38:44, 53:59, 69:74),1] <- "normal(0,5)"


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 16000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472
THREADS <- 2


# Prior Predictive Checks 
# Note that due to weakly informative priors, the prior predictive checks reveal little about the 
# ability of our priors to describe the data. 
countmodel_priorpredictive <- brm(
  mv,
  data = count_data,           
  prior = countmodel_mv_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = 16000,
  warmup = 1600,
  seed = SEED,
  backend = "cmdstanr",
  sample_prior = "only",
  control = list(adapt_delta = 0.95)
)


# Fit model to data
countmodel_mv_fit <- brm(
  mv,
  data = count_data,           
  prior = countmodel_mv_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = 16000,
  warmup = 1600,
  seed = SEED,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97))


# Post-fitting diagnostics (including inspection of ESS, Rhat and posterior predictive)
tidy_countmodel_fit <- tidy(countmodel_mv_fit)
#posteriorpredictive <-pp_check(countmodel_fit, ndraws = 500)


# Misc evaluations
#print(prior_summary(countmodel_fit), show_df = FALSE)
performance(countmodel_mv_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(countmodel_mv_fit) # default output plot of brms showing posterior distributions of
prior_summary(countmodel_mv_fit) #obtain dataframe of priors used in model.

mcmc_countmodel_fit <- ggs(countmodel_mv_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion

posteriorpredictive_group <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'reassortantclass')
posteriorpredictive_nreassortants <-pp_check(countmodel_mv_fit, ndraws = 500, resp = 'nreassortants')





############################################## WRITE ###############################################
saveRDS(countmodel_mv_fit, file = './saved_models/count_model.rds')
write_csv(count_data, './saved_models/count_data.csv')


############################################## END #################################################
####################################################################################################
####################################################################################################
# Deprecated

# Marginal probability density of the number of unique reassortants/region
# ie, irrespective of class
count_prob <- jointmodel_temp %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),resp = "nreassortants",
              re_formula = NA)

count_prob %>% 
  median_hdi(.epred)

ggplot(count_prob, aes(x = .epred, fill = collection_regionname , y = collection_regionname)) +
  stat_halfeye() +
  theme_minimal() + 
  scale_x_continuous('Number of Unique Reassortants/Year') +
  
  # Scales
  scale_y_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values = region_colour %>% pull(Trait, name = Name))  +
  theme(legend.position = 'none') + 
  coord_cartesian(xlim = c(0,10))

# Marginal probability of reassortant class/region
# ie, irrespective of number
class_prob <- jointmodel_temp %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),
              resp = "group2",
              re_formula = NA)
class_prob %>% 
  median_hdi(.epred)

ggplot(class_prob, aes(x = collection_regionname, colour = .category, y = .epred)) +
  geom_boxplot() +
  theme_minimal() + 
  scale_y_continuous('Posterior Probability of Reassortant Class', labels = scales::percent) +
  
  # Scales
  scale_x_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values = riskgroup_colour %>% pull(Trait, name = group2)) +
  theme(legend.position = 'bottom')


#### Contrasts ####
fit_year %>% 
  emmeans(., ~ collection_regionname ,
          epred = TRUE,
          resp = 'group2') %>% 
  filter(rep.meas == 'nreassortants') %>%
  contrast(method = 'revpairwise')



joint_prob <- count_prob %>%
  left_join(class_prob, 
            by = c(".draw",  '.row', "collection_regionname"), 
            suffix = c("_count", "_class")) %>%
  mutate(joint_epred = .epred_count * .epred_class)