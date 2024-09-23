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
  mutate(collection_month = date_decimal(TMRCA) %>% format(., "%m"),
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




####################################### START BRMS PIPELINE ########################################

# Set Priors
diffusionmodel1_priors <- c()


# Set MCMC Options
CHAINS <- 4
CORES <- 4
ITER <- 4000
BURNIN <- ITER/10 # Discard 10% burn in from each chain
SEED <- 4472


# Prior Predictive Checks 
# what does weighted diffusion coefficient measure specifically (ie speed vs geographical distance)
diffusionmodel1_prior <- brm(
  bf(weighted_diff_coeff ~ 1 + host_richness + median_anseriformes_wild + median_charadriiformes_wild + (1|segment + collection_regionname),
     hu ~ 1 + host_simplifiedhost + season + (1|segment + collection_regionname)),
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



color_scheme_set("green")
plot_priorpredictive <- pp_check(diffusionmodel1_priorpredictive, ndraws = 100) + 
  theme_minimal()

plot_priorpredictive$data %<>%
  mutate(value = log1p(value))

plot_priorpredictive




# Fit Model
diffusionmodel1_fit <- brm(
  bf(weighted_diff_coeff ~ 1 + host_richness + median_anseriformes_wild + median_charadriiformes_wild + (1|segment + collection_regionname),
     hu ~ 1 + season + (1|segment + collection_regionname)),
  data = diffusion_data,
  family = hurdle_lognormal(),
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.99))


color_scheme_set('red')
plot_posteriorpredictive <- pp_check(diffusionmodel1_fit, ndraws = 100) + 
  theme_minimal()

plot_posteriorpredictive$data %<>%
  mutate(value = log1p(value))

plot_posteriorpredictive
# Conditional effect of host mrca on hurdle
# predicted probability of diffusion coefficient == 0
plot(conditional_effects(diffusionmodel1_fit, dpar = "hu"), plot = FALSE)[[4]]$data %>%
  ggplot(aes(x = effect1__, y = estimate__)) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = lower__, ymax = upper__, x = effect1__)) +
  scale_y_continuous('Predicted Probability Of Diffusion Coefficient == 0') + 
  theme_minimal() 
  
cond <- plot(conditional_effects(diffusionmodel1_fit, dpar = "mu"), plot = FALSE)[-4] %>%
  lapply(., function(x)x$data) %>%
  bind_rows(, .id = 'var') %>%
  filter(var != 'median_anseriformes_wild')  %>%
  ggplot(.,aes(x = effect1__, y = estimate__, colour = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__ , ymax = upper__,, x = effect1__, fill = var), colour = NA, alpha = 0.2) + 
  scale_y_continuous('Dispersal Coefficient') + 
  theme_minimal() +
  coord_cartesian(xlim = c(0,5)) +
  theme(legend.position = 'bottom')

diffusionmodel1_fit |> emtrends(
  ~host_richness,
  var = 'host_richness',
  regrid = "response",  
  tran = 'log',
  type = "response",
  at = list(host_richness = seq(1, 5, 1)),
  dpar = "mu") |> 
  gather_emmeans_draws() |> 
  ggplot(aes(x = host_richness, y = .value)) +
  stat_lineribbon(size = 1)

hurdle_intercept <- tidy(diffusionmodel1_fit) |> 
  filter(term == "hu_(Intercept)") |> 
  pull(estimate)

hurdle_season <- tidy(diffusionmodel1_fit) |> 
  filter(term == "hu_seasonoverwintering") |> 
  pull(estimate)

plogis(hurdle_intercept + hurdle_season) - plogis(hurdle_intercept)


diffusionmodel1_fit |> 
  emtrends(~ host_richness, var = "median_charadriiformes_wild", dpar = "mu",
           at = list(host_richness = seq(1, 5, 1)))



test <- diffusionmodel1_fit %>% 
  epred_draws(newdata = expand_grid(host_richness = 0,
                                    median_anseriformes_wild = 0,
                                    median_charadriiformes_wild = 0,
                                    season = c('overwintering'),
                                    region = unique(diffusion_data$collection_regionname)[-5]), 
            re_formula = ~ (1|region))  # or re_formula = ~ (1 | region)

test <- diffusionmodel1_fit %>% 
  epred_draws(newdata = expand_grid(host_richness = seq(0,6, by = 1),
                                    median_anseriformes_wild = 0,
                                    median_charadriiformes_wild = 0,
                                    season = c('overwintering'),
                                    region = unique(diffusion_data$collection_regionname)[-5]), 
              re_formula =~ (1|region))

ggplot(test, 
       aes(x = host_richness, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Host Richness", y = "Predicted Diffusion Coefficient",
       fill = "Credible interval") +
  facet_wrap(vars(region))

ggplot(test, 
       aes(x = .epred, 
           fill = region)) +
  stat_halfeye() + 
  scale_x_continuous(trans = 'log')
  coord_cartesian(xlim = c(0, 100000))

diffusionmodel1_fit %>%
  predicted_draws(newdata = tibble(collection_regionname = ,
                                   segment = 'ha'))
  emtrends(~ lifeExp, var = "lifeExp", dpar = "mu",
           at = list(lifeExp = seq(30, 80, 10)))


diffusionmodel1_fit %>% 
  emmeans(~ collection_regionname,
          epred = TRUE,  dpar = "hu") %>% 
  contrast(method = "revpairwise")

# Conditional effect of time in charadrifformes, and anseriformes on dispersal velocity
diffusionmodel1_fit |> 
  emtrends(~  median_anseriformes_wild,
           var = " median_anseriformes_wild",
           dpar = "mu", 
           regrid = "response",
           tran = "log", 
           type = "response",
           at = list( median_anseriformes_wild = seq(0, 10, 2)))


diffusionmodel1_fit|> marginaleffects::comparisons(
  newdata = datagrid( median_anseriformes_wild = seq(0, 10, 2)),
  dpar = "mu", 
  transform_pre = "expdydx"
)

model_gdp_hurdle_life %>%
  emmeans(~ lifeExp, var = "lifeExp", epred = TRUE,
          at = list(lifeExp = seq(30, 80, 1))) |> 
  gather_emmeans_draws() |> 
  ggplot(aes(x = lifeExp, y = .value)) +
  stat_lineribbon(size = 1, color = clrs[4]) +
  scale_fill_manual(values = colorspace::lighten(clrs[4], c(0.95, 0.7, 0.4))) +
  scale_y_continuous(labels = label_dollar()) +
  labs(x = "Life expectancy", y = "Predicted GDP per capita",
       subtitle = "Expected values from mu and hu parts (epred = TRUE)",
       fill = "Credible interval")

t <- plot(conditional_effects(diffusionmodel1_fit, dpar = "mu"), plot = FALSE)[[4]]$data %>% median_anseriformes_wild
  ggplot(aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.3, color = NA) +
  geom_line(size = 1, color = clrs[5]) +
  scale_y_continuous(labels = label_percent()) +
  labs(x = "Life expectancy", y = "Predicted probability\nof seeing $0 GDP per capita",
       subtitle = "Hurdle part of the model (dpar = \"hu\")") +
  theme_nice()



####################################################################################################
####################################################################################################