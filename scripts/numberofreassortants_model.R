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

numberofsequence_data <- read_csv('./2024-09-09_meta.csv') %>%
  select(c(collection_datemonth,
           collection_regionname,
           cluster_profile)) %>%
  drop_na(collection_datemonth) %>%
  mutate(collection_regionname = case_when(grepl('europe', 
                                                 collection_regionname) ~ 'europe',
                                           grepl('africa',
                                                 collection_regionname) ~ 'africa',
                                           grepl('asia', 
                                                 collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', 
                                                 collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname)) %>%
  mutate(collection_year = 
           paste0(collection_datemonth, '-01') %>%
           ymd() %>%
           decimal_date() %>%
           round(digits = 1)) %>%
  summarise(n = n(), .by = c(collection_year, collection_regionname))


count_data <- combined_data %>%
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
                    group2))  %>%
  mutate(collection_year = collection_monthyear %>%
           paste0('-01') %>%
           ymd() %>%
           decimal_date() %>%
           round(., digits = 1))  %>%
  summarise(n_reassortants = sum(n_reassortants), 
            .by = c(collection_year,collection_regionname, group2)) %>%
  
  arrange(collection_year) %>%
  left_join(numberofsequence_data)


# Model Formula
bf1 <- bf(n_reassortants|resp_trunc(lb = 1) ~ 1 + collection_regionname + (1|collection_year),
          family = poisson())

bf2 <- bf(group2 ~ 1 + collection_regionname + (1|collection_year), 
          family = categorical())

count_formula <- mvbf(bf1, bf2)


# Define Priors
countmodel_priors <- get_prior(count_formula,
                               data = count_data,
                               family = list(zero_inflated_negbinomial(), categorical())) 


countmodel_priors$prior[c(1:6,11:14,19:22)] <- "normal(0,5)"


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472
THREADS <- 2


# Prior Predictive Checks 
# Note that due to weakly informative priors, the prior predictive checks reveal little about the 
# ability of our priors to describe the data. 
countmodel_priorpredictive <- brm(
  count_formula,  # Multivariate outcome
  data = count_data,           # Dataframe containing y1, y2, x1, x2
  prior = countmodel_priors,
  chains = CHAINS,
  cores = CORES, 
  threads = THREADS,
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  sample_prior = "only",
  control = list(adapt_delta = 0.95)
)


# Fit model to data
countmodel_fit <- brm(
  count_formula,  # Multivariate outcome
  data = count_data,           # Dataframe containing y1, y2, x1, x2
  prior = countmodel_priors,
  chains = CHAINS,
  threads = THREADS, 
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95)
)

# Post-fitting diagnostics (including inspection of ESS, Rhat and posterior predictive)
tidy_countmodel_fit <- tidy(countmodel_fit)
#posteriorpredictive <-pp_check(countmodel_fit, ndraws = 500)


# Misc evaluations
#print(prior_summary(countmodel_fit), show_df = FALSE)
performance(countmodel_fit) # tibble output of model metrics including R2, ELPD, LOOIC, RMSE
plot(countmodel_fit) # default output plot of brms showing posterior distributions of
prior_summary(countmodel_fit) #obtain dataframe of priors used in model.

mcmc_countmodel_fit <- ggs(countmodel_fit) # Warning message In custom.sort(D$Parameter) : NAs introduced by coercion

posteriorpredictive_group <-pp_check(countmodel_fit, ndraws = 500, resp = 'group2')
posteriorpredictive_nreassortants <-pp_check(countmodel_fit, ndraws = 500, resp = 'nreassortants')


# Posterior Predictions and Marginal Effects 
# Marginal probability of N reassortants / region 


# posterior probabilities of n reassortants stratified by country
# while ignoring any  deviations of the intercept or slope
n_reassortant_preds <- countmodel_fit %>% 
  predicted_draws(newdata = no_zeros %>%
                    select(collection_regionname) %>%
                    distinct(),
                  resp = "nreassortants",
                  value = 'n',
                  re_formula = NA) 

test <- n_reassortant_preds %>%
  mutate(n = as.character(n)) %>%
  group_by( collection_regionname) %>%
  count(n) %>%
  reframe(freq = nn/sum(nn), n= n) %>%
  ungroup() 

ggplot(test) +
  geom_bar(aes(x = n, y= freq), stat = 'identity') + 
  facet_grid(cols = vars(collection_regionname)) +
  theme_minimal()


# posterior probability of class reassortants in an average year
# while ignoring any  deviations of the intercept or slope

countmodel_fit %>% 
  predicted_draws(newdata = no_zeros %>%
                    select(collection_regionname) %>%
                    distinct(),
                  resp = "group2",
                  value = 'class',
                  re_formula = NA) %>%
  group_by(collection_regionname) %>%
  count(class) %>%
  reframe(freq = n/sum(n), class = class) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = class, y= freq), stat = 'identity') + 
  facet_grid(cols = vars(collection_regionname)) +
  theme_minimal()

countmodel_fit %>% 
  epred_draws(newdata = no_zeros %>%
                select(collection_regionname) %>%
                distinct(),
              resp = "nreassortants",
              value = 'n',
              re_formula = NA) %>%
  ggplot( aes(x = n, y = collection_regionname, height = after_stat(density))) + 
  geom_density_ridges(scale = 0.95, draw_baseline = FALSE)


class_preds <- countmodel_fit %>% 
  epred_draws(newdata = no_zeros %>%
                select(collection_regionname) %>%
                distinct(),
              resp = "group2",
              value = 'p',
              re_formula = NA) 

median_hdci(class_preds) 

ggplot(class_preds) + 
  geom_density_ridges(aes(x = p,
                          y = collection_regionname, 
                          fill = .category,
                          colour = .category),
                      alpha = 0.5, 
                      rel_min_height = 0.01,
                      scale=0.8)+
  theme_minimal() + 
  scale_x_continuous(limits = c(0,1))


# Required outputs: Pre/post epizootic, marginal effect of region (ie how many more reassortants in a vs b),
# 


# Marginal probability of N reassortants / year ~ 1|Region


# Marginal effect of region


# Marginal class probability  


# Marginal class probabilty ~ 1| Region


# Joint Probability



# Joint Probability ~ 1|Region


# Marginal effect of region




############################################## WRITE ###############################################




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
                distinct(),resp = "group2",
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