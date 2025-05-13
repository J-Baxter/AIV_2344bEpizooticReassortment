####################################################################################################
####################################################################################################
## Script name: Class Model (using cumulative distribution)
##
## Purpose of script:to model the probability that any given reassortant belongs to class X,
## stratified by continent
##
## Date created: 2025-05-13
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)
library(broom)
library(broom.mixed)
library(brms)
library(cmdstanr)


# User functions


############################################## DATA ################################################



############################################## MAIN ################################################
temp <- brm(cluster_class~1+parent_class+time_since_parent*segments_changed,
            data=reassortant_ancestral_changes %>%
              mutate(cluster_class = ordered(cluster_class, levels = c('minor', 'moderate', 'major'))),
            family =cumulative("probit"),
            backend = 'cmdstanr')

marginal_effects(temp, "parent_class", categorical = TRUE)
marginal_effects(temp, "time_since_parent", categorical = TRUE)
marginal_effects(temp, "segments_changed", categorical = TRUE)


te



############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################