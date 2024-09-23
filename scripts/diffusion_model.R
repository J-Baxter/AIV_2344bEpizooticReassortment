####################################################################################################
####################################################################################################
# The aim of this model is to determine associations between variables obtained from our phylogenetic
# analysis and difusion coefficient for each reassortant

# variables of interest include: 
# 1. evolutionary rates
# 2. persistence time
# 3. region of origin
# 4. persistence time in wild birds
# 5. persistence time in domestic birds
# 6. total number of species jumps

# This script implements a workflow for a zero-inflated lognormal model. Initial exploratory analysis 
# is conducted using OLS regression, thereafter progressing to Bayesian analysis using BRMS.

########################################## DEPENDENCIES ############################################
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(emmeans)
library(marginaleffects)
library(magrittr)

########################################### IMPORT DATA ############################################
combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

########################################### FORMAT DATA ############################################

diffusion_data <- combined_data %>%
  
  # select variables of interes
  select(
    segment,
    cluster_profile,
    TMRCA,
    group2,
    weighted_diff_coeff,
    original_diff_coeff,
    evoRate,
    persist.time,
    collection_regionname,
    host_simplifiedhost,
    count_cross_species,
    starts_with('median'),
    starts_with('max')) %>%
  
  # Substitute NA values in diffusion coefficient with 0
  mutate(across(where(is.double), .fns = ~ replace_na(.x, 0))) %>%
  
  rename_with(~gsub('-', '_', .x)) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  filter(!grepl('\\+', host_simplifiedhost)) %>%
  
  # join host richness
  left_join(summary_data %>% select(c(cluster_profile, 
                                      host_richness)),
            by = join_by(cluster_profile)) %>%

  # season (breeding, migrating_north, migrating_south, overwintering)
  mutate(collection_month = date_decimal(TMRCA) %>% format(., "%m") %>% as.integer(),
         season = case_when(collection_month %in% c(12,1,2) ~ 'overwintering', 
                            collection_month %in% c(3,4,5)  ~ 'migrating_north', 
                            collection_month %in% c(6,7,8)  ~ 'breeding', 
                            collection_month %in% c(9,10,11)  ~ 'migrating_south'
                            ))

  
################################### INITIAL EXPLORATORY MODELS #####################################
# Plot logged data to show zero/non-zero segregation

model_data %>%
  ggplot() +
  geom_histogram(aes(x = log1p(weighted_diff_coeff), fill = weighted_diff_coeff>0)) +
  scale_fill_brewer(palette = 'Dark2', 'Is Zero') +
  scale_x_continuous('Weighted Diffusion Coefficient') +
  theme_minimal()


# OLS model of 

######################################## DEFINE FORMULA ############################################
# We assume a hurdle lognormal model, in which the hurdle process is determined the season in which 
# the TMRCA is estimated to be, and the lognormal component is determined by host richness, persistence
# in anseriformes and charadriiformes. Both model components are conditional on segment from which 
# the measurement is taken and region of origin

diffusion_formula <- bf(weighted_diff_coeff ~ 1 + host_richness + median_anseriformes_wild +  median_charadriiformes_wild + 
                          (1|segment + collection_regionname),
                        hu ~ 1 +  season + (1|segment + collection_regionname))
####################################### DEFINE PRIORS ########################################

# Set Priors
diffusionmodel1_priors <- c()



####################################### SET MCMC OPTION ########################################

# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472


###################################### PRIOR PREDICTIVE SIM ########################################

# Prior Predictive Checks 

diffusionmodel1_prior <- brm(
  ,
  data = diffusion_data,
  family = hurdle_lognormal(),
  sample_prior = "yes",
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95)
  )


# PRIOR PREDICTIVE CHECKS
color_scheme_set("green")
plot_priorpredictive <- pp_check(diffusionmodel1_priorpredictive, ndraws = 100) + 
  theme_minimal()

plot_priorpredictive$data %<>%
  mutate(value = log1p(value))

plot_priorpredictive



############################################ FIT MODEL #############################################
diffusionmodel1_fit <- brm(
  bf(weighted_diff_coeff ~ 1 + host_richness + median_anseriformes_wild + median_charadriiformes_wild + (1|segment + collection_regionname),
     hu ~ 1 + season + (1|segment + collection_regionname)),
  data = diffusion_data,
  family = hurdle_lognormal(),
  chains = CHAINS,
  cores = CORES, 
  threads = 2, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.99))

# CONVERGENCE CHECK 


# POSTERIOR PREDICTIVE CHECKS
color_scheme_set('red')
plot_posteriorpredictive <- pp_check(diffusionmodel1_fit, ndraws = 100) + 
  theme_minimal()

plot_posteriorpredictive$data %<>%
  mutate(value = log1p(value))

plot_posteriorpredictive


##################################### PREDICTIONS AND EFFECTS ######################################


# predicted probability of diffusion coefficient == 0 per season
condition_data <- plot(conditional_effects(diffusionmodel1_fit, dpar = "hu"),
                       plot = FALSE,
                       re_formula = NULL)[[4]]$data


# Can we compare across random effect region? <<<<<<<< OUTSTANDING TASK


# Marginal effects of host richness, and host persistence on dispersal velocity (mu part of model)
# does not determine = 0 / !=0 process

# This represents the instantaneous slopes at each of these values of host richness in 'mu' part of
# the model.
# e.g at host richness = 3, an +1 increase in host richness is associated with an additional 19597 of
# diffusion coefficient
hostrichness_marginal <-  diffusionmodel1_fit %>%
  emtrends(~ host_richness,
           var = "host_richness", 
           dpar = "mu",
           at = list(host_richness = seq(1, 5, 1)))

# at response scale
hostrichness_marginal_response <-  diffusionmodel1_fit %>%
  emtrends(~ host_richness,
           var = "host_richness", 
           dpar = "mu",
           type = 'response',
           tran ='log',
           regrid = 'response',
           at = list(host_richness = seq(1, 5, 1)))


# This represents the instantaneous slopes at each of these values of charadriiformeswild in 'mu' part of
# the model.
charadriiformeswild_marginal <-  diffusionmodel1_fit %>%
  emtrends(~ median_charadriiformes_wild,
           var = "median_charadriiformes_wild", 
           dpar = "mu",
           at = list(median_charadriiformes_wild = seq(1, 5, 1)))

# at response scale
charadriiformeswild_marginal_response <-  diffusionmodel1_fit %>%
  emtrends(~ median_charadriiformes_wild,
           var = "median_charadriiformes_wild", 
           dpar = "mu",
           type = 'response',
           tran ='log',
           regrid = 'response',
           at = list(median_charadriiformes_wild = seq(1, 5, 1)))


# This represents the instantaneous slopes at each of these values of anseriformeswild in 'mu' part of
# the model.
anseriformeswild_marginal<-  diffusionmodel1_fit %>%
  emtrends(~ median_anseriformes_wild,
           var = "median_anseriformes_wild", 
           dpar = "mu",
           at = list(median_anseriformes_wild = seq(1, 5, 1)))

# at response scale
anseriformeswild_marginal_response <-  diffusionmodel1_fit %>%
  emtrends(~ median_anseriformes_wild,
           var = "median_anseriformes_wild", 
           dpar = "mu",
           type = 'response',
           tran ='log',
           regrid = 'response',
           at = list(median_anseriformes_wild = seq(1, 5, 1)))


##################################### PAIRWISE COMPARISONS ######################################

# pair wise comparisons of season on probability of diffusion coefficient == 0
diffusionmodel1_fit %>%
  emmeans(~season,
          var = 'season',
          dpar = 'hu')%>%
  contrast(method = "revpairwise",
           ref = 'migrating_south', 
           type = 'response')

############################################### PLOTS ############################################

# Conditional effect of MRCA season on the probability that a reassortant is observed to spread at all
# (ie hurdle component)
plot(conditional_effects(diffusionmodel1_fit, dpar = "hu"), plot = FALSE)[[4]]$data %>%
  ggplot(aes(x = effect1__, y = estimate__)) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = 1-lower__, ymax = 1-upper__, x = effect1__)) +
  scale_y_continuous('Predicted probability that diffusion coefficient > 0') + 
  scale_x_discrete('Season of MRCA',
                   labels = c('breeding' = 'Breeding',
                              'migrating_north' = 'Spring Migration',
                              'migrating_south' = 'Autumn Migration',
                              'overwintering' = 'Overwintering')) +
  theme_minimal(base_size = 18) 


# Conditional effect of median persistence time (anser and charadriformes) and host richness on
# the diffusion coefficient
plot(conditional_effects(diffusionmodel1_fit, dpar = "mu"), plot = FALSE)[-4] %>%
  lapply(., function(x)x$data) %>%
  bind_rows(, .id = 'var') %>%
  #filter(var != 'median_anseriformes_wild')  %>%
  ggplot(.,aes(x = effect1__, y = estimate__, colour = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__ , ymax = upper__,, x = effect1__, fill = var), colour = NA, alpha = 0.2) + 
  scale_y_continuous('Dispersal Coefficient') + 
  theme_minimal() +
  coord_cartesian(xlim = c(0,5)) +
  theme(legend.position = 'bottom')


####################################################################################################
####################################################################################################