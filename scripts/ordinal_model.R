

library(brms)

temp <- brm(cluster_class~1+parent_class+time_since_parent*segments_changed,
    data=reassortant_ancestral_changes %>%
      mutate(cluster_class = ordered(cluster_class, levels = c('minor', 'moderate', 'major'))),
    family =cumulative("probit"),
    backend = 'cmdstanr')

marginal_effects(temp, "parent_class", categorical = TRUE)
marginal_effects(temp, "time_since_parent", categorical = TRUE)
marginal_effects(temp, "segments_changed", categorical = TRUE)


te
