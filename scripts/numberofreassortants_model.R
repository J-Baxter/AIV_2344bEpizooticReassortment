####################################################################################################
####################################################################################################
# The aim of this model is to determine whether reassortants emerging from some regions are more
# likely and more severe than those from others

# variables of interest include: 
# 1. region of origin
# 3. persistence time (persistence time increases uncertainty)
# 4. number of sequences per region
# 5. diveresity measure?

# This script implements a workflow for a negative binomial model. 

########################################## DEPENDENCIES ############################################
library(brms)
library(broom)
library(broom.mixed)
library(tidybayes)
library(bayesplot)


########################################### IMPORT DATA ############################################

combined_data <- read_csv()


########################################### FORMAT DATA ############################################

count_data <- combined_data %>%
  filter(segment == 'ha') %>%
  mutate(collection_monthyear = date_decimal(TMRCA) %>% format(., "%Y-%m")) %>%
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  filter(!is.na(collection_regionname)) %>%
  select(collection_regionname, collection_monthyear, cluster_profile, group2) %>% # include all month-years over collection period to generate zer countr
  summarise(n_reassortants = n_distinct(cluster_profile), 
            .by = c(collection_monthyear, collection_regionname, group2)) %>%
  mutate(collection_monthyear = ym(collection_monthyear)) %>%
  group_by(collection_regionname, group2) %>%
  complete(collection_monthyear = seq(min(collection_monthyear), max(collection_monthyear), by = "month")) %>% 
  replace_na(list(n_reassortants = 0)) %>%
  ungroup() %>%
  mutate(collection_year = collection_monthyear %>% format(., "%Y") %>% as.integer()) %>%
  summarise(n_reassortants = sum(n_reassortants), .by = c(collection_year, collection_regionname , group2)) %>%
  arrange(collection_year) 
  


ggplot(grouped_lines) +
  geom_bar(aes(x = collection_monthyear, y = n_reassortants, fill = group2), stat = 'identity') + 
  facet_grid(rows = vars(collection_regionname),  labeller = labeller(collection_regionname = label_wrap_gen(10))) +
  scale_y_continuous(breaks = seq(0,14, by = 2), 'Number of Reassortants')+
  scale_x_date('Resassortant MRCA') +
  scale_fill_manual(
    'Reassortant Class',
    values = riskgroup_colour %>% pull(Trait, name = group2)) +
    #labels = riskgroup_colour %>% pull(group2)) +
  theme_minimal()  +
  theme(legend.position = 'bottom')


######################################## DEFINE FORMULA ############################################
nor_formula <- mvbind(n_reassortants, group2) ~ 1 + collection_regionname + (1|collection_year)


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


# Prior Predictive Checks 


###################################### PRIOR PREDICTIVE SIM ########################################

# Multivariate model - no correlation
multivar_model <- brm(
  formula_y1 + formula_y2,  # Joint model formula
  data = grouped_lines,           # Dataframe containing y1, y2, x1, x2
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  sample_prior = "yes",
  control = list(adapt_delta = 0.95)
)


############################################ FIT MODEL #############################################

fit_year <- brm(
  nor_formula,  # Multivariate outcome
  family = list(zero_inflated_negbinomial(), categorical()),  # Count and class families
  data = count_data,           # Dataframe containing y1, y2, x1, x2
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  control = list(adapt_delta = 0.95)
)



# CONVERGENCE CHECK 


# POSTERIOR PREDICTIVE CHECKS


##################################### PREDICTIONS AND EFFECTS ######################################
# to-do: 
# re-plot probability of reassortant class output
# update plot ( col scheme) for number of reassortants/region -done
# joint probability of class & reassortants (similar to thompson et al plot?)


# Marginal probability density of the number of unique reassortants/region
# ie, irrespective of class
count_prob <- fit_year %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),resp = "nreassortants",
              re_formula = NA)

count_prob %>% 
  median_hdi(.epred)

ggplot(count_prob, aes(x = .epred, fill = collection_regionname , y = collection_regionname)) +
  stat_halfeye() +
  theme_minimal(base_size = 18) + 
  scale_x_continuous('Number of Unique Reassortants/Year') +
  
  # Scales
  scale_y_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_fill_manual(
    values = c("africa" = "#CC2929CC", 
               "asia" = "#ABCC29CC", 
               "europe" = "#29CC6ACC", 
               "central & northern america" = "#296ACCCC"))  +
  theme(legend.position = 'none') + 
  coord_cartesian(xlim = c(0,10))




scale_fill_manual()
# Marginal probability of reassortant class/region
# ie, irrespective of number
class_prob <- fit_year %>% 
  epred_draws(newdata = count_data %>%
                select(collection_regionname) %>%
                distinct(),resp = "group2",
              re_formula = NA)
class_prob %>% 
  median_hdi(.epred)


ggplot(class_prob, aes(x = .epred, colour = .category , y = collection_regionname)) +
  stat_halfeye(position = position_dodge()) +
  theme_minimal(base_size = 18) + 
  scale_x_continuous('Posterior Probability of Reassortant Class') +
  
  # Scales
  scale_y_discrete('Region of Origin', labels = function(x) str_wrap(x, width = 20) %>% str_to_title())+
  scale_colour_manual(
    'Reassortant Class',
    values =   c('minor'  =  "#2ca02c",
                 'major' = "#1f77b4",
                 'dominant' = "#FF0000" )) +
  theme(legend.position = 'bottom')


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


# joint predict
pred_count <- posterior_epred(fit_year, resp = "nreassortants") #SxN matrix = S is the number of posterior draws, N is the number of observations
pred_class <- posterior_epred(fit_year, resp = "group2")

class_1_probs <- pred_class[, , 1]  # Extract probabilities for class 2 across all draws and observations
class_2_probs <- pred_class[, , 2]  # Extract probabilities for class 2 across all draws and observations
class_3_probs <- pred_class[, , 3]  # Extract probabilities for class 2 across all draws and observations
count_n_probs <- pred_count         # Use all predicted counts for all observations <- this is wrong?

# Element-wise multiplication for all rows (posterior draws) and columns (observations)
joint_probs_class1 <- count_n_probs * class_1_probs
joint_probs_class2 <- count_n_probs * class_2_probs
joint_probs_class3 <- count_n_probs * class_3_probs


##################################### PAIRWISE COMPARISONS ######################################
fit_year %>% 
  emmeans(., ~ collection_regionname ,
          epred = TRUE,
          resp = 'group2') %>% 
  filter(rep.meas == 'nreassortants') %>%
  contrast(method = 'revpairwise')

fit_year %>%
  emmeans(.,~collection_regionname,
          epred = TRUE,
          resp = 'group2')%>%
  contrast(method = "revpairwise",
           type = 'response')

joint_prob <- count_prob %>%
  left_join(class_prob, 
            by = c(".draw",  '.row', "collection_regionname"), 
            suffix = c("_count", "_class")) %>%
  mutate(joint_epred = .epred_count * .epred_class)





############################################### PLOTS ############################################


numberreassortants_priorpredictive <- brm(
  n_reassortants|trunc(lb = 1) ~ 1 + collection_regionname,
  data = grouped_lines,
  family = negbinomial(), 
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED,
  sample_prior = "yes",
  control = list(adapt_delta = 0.95)
)

color_scheme_set("green")
plot_priorpredictive <- pp_check(fit_year, ndraws = 100) + 
  theme_minimal()

plot_priorpredictive$data %<>%
  mutate(value = log1p(value))

plot_priorpredictive

color_scheme_set("red")
fit_priorpredictive <- pp_check(fit_year, ndraws = 100, resp = 'nreassortants') + 
  theme_minimal()



fit_priorpredictive$data %<>%
  mutate(value = log1p(value))

fit_priorpredictive

# Fit Model
diffusionmodel1_fit <- brm(
  n_reassortants ~ n_sequences + collection.region.name,
  data = grouped_lines,
  family = negbinomial(),
  chains = CHAINS,
  cores = CORES, 
  iter = ITER,
  warmup = BURNIN,
  seed = SEED#,
  #opencl = opencl(c(0, 0)) # Enables GPU accelerated computation, remove if not applicable
)


# Posterior Predictive Checks
diffusionmodel1_posteriorpreds <- posterior_predict(diffusionmodel1_fit)
n <- sample(1:nrow(diffusionmodel1_posteriorpreds), 50)

color_scheme_set("green")
ppc_dens_overlay(y = log1p(model_data$weighted_diff_coeff),
                 yrep = log1p(diffusionmodel1_posteriorpreds[1:10,]))


# Extract Model Terms and Marginal/Conditional Effects
diffusionmodel1_fit_tidy <- tidy(model_hurdle)


####################################################################################################
####################################################################################################

#deprecated 
bayes_mod <- 
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
  stat_halfeye(point_interval = mean_hdi,.width = c(.90, .95)) +
   +
  theme_minimal(base_size = 18) 

reassortant_metadata_formatted %>%
  mutate(collection.date = ymd(collection.date)) %>%
  ggplot(aes(x = collection.date, y = collection.region.name, fill = as.factor(cluster.genome)))+
  geom_density()+
  facet_wrap(.~collection.region.name)