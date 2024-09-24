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
library(broom)
library(broom.mixed)
library(tidyverse)
library(recipes)


########################################### IMPORT DATA ############################################

combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

# import colour schemes
host_colour <- read_csv('./colour_schemes/hostType_cols.csv')
region_colour <- read_csv('./colour_schemes/regionType_cols.csv')
riskgroup_colour <- read_csv('./colour_schemes/riskgroup_cols.csv') %>%
  mutate(kclust = c(2,1, 3),
         group2 = c('major', 'minor', 'dominant'))
subtype_colour <- read_csv('./colour_schemes/SubType_cols.csv')

########################################### FORMAT DATA ############################################
kclust_nophylo_data <- summary_data %>%
  
  # criteria selected by Lu
  select(c(
    cluster_profile,
    group2,
    host_richness,
    if_mammal,
    max_distance_km,
    Num_Sequence,
    number_Conti,
    Length_Between_First_Last_Sample
  )) %>% 
  
  # All numeric NA as 0
  mutate(across(where(is.double), .fns = ~ replace_na(.x, 0))) %>%
  
  # Start tidymodels recipe
  recipe(~ .) %>% 
  
  # normalise numeric vectors
  step_normalize(all_numeric()) %>% 
  
  # create dummy variables for categorical predictor
  #step_dummy(Continent_of_Earliest_Date, one_hot = TRUE) %>%
  #step_dummy(Host_of_Earliest_Sample_Date, one_hot = TRUE) %>%
  
  # run tidymodels recipe
  prep() %>%
  bake(NULL) %>%
  
  # exclude columns not to be used
  #select(-c(cluster_profile, original_diff_coeff, starts_with('median'))) %>%
  drop_na()
  


kclust_updated_data <- combined_data %>%
  
  # select variables of interes
  select(
    segment,
    cluster_profile,
    group2,
    weighted_diff_coeff,
    original_diff_coeff,
    evoRate,
    persist.time,
    collection_regionname,
    host_simplifiedhost,
    count_cross_species,
    starts_with('median'),
    starts_with('max')
    ) %>%
  
  
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
  
  # Model pre-processing
  group_by(segment, cluster_profile) %>%
  slice_sample(n=1) %>%
  ungroup() %>%
  
  # start tidymodels
  recipe(~ .) %>% 
  
  # normalise numeric vectors
  step_normalize(all_numeric()) %>% 
  
  # create dummy variables for categorical predictors
  step_dummy(collection_regionname, one_hot = TRUE) %>%
  step_dummy(host_simplifiedhost, one_hot = TRUE) %>%
  step_dummy(segment, one_hot = TRUE) %>%
  
  # run tidymodels recipe
  prep() %>%
  bake(NULL) %>%
  
  # exclude columns not to be used
  left_join(kclust_nophylo_data %>% select(-starts_with('group')), by = join_by(cluster_profile)) %>%
  select(-c(cluster_profile, original_diff_coeff, starts_with('median'))) %>%
  drop_na() 


####################################### START KCLUST PIPELINE ########################################
set.seed(4472)
kclust_nophylo <- tibble(k = 1:20) %>%
  mutate(
    kclust = map(k, ~kmeans(kclust_nophylo_data %>% select(-c(group2, cluster_profile)), .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, kclust_nophylo_data)
  )

clusters_nophylo <- 
  kclust_nophylo %>%
  unnest(cols = c(tidied))

assignments_nophylo <- 
  kclust_nophylo %>% 
  unnest(cols = c(augmented))

clusterings_nophylo <- 
  kclust_nophylo %>%
  unnest(cols = c(glanced))

clustersnophylo_plot <- assignments_nophylo %>%
  filter(k == 4) %>%
  select(-c(kclust, tidied, glanced)) %>%
  recipe(~ .) %>%
  step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
  prep() %>%
  bake(NULL)%>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = .cluster), size = 4, alpha = 0.8) + 
  theme_minimal(base_size = 20) + 
  scale_color_manual(
    'K-Means Clusters',
    values = riskgroup_colour %>% pull(Trait, name = kclust), 
    labels = riskgroup_colour %>% pull(group2, name = kclust)) +
  theme(legend.position = 'bottom',
        legend.box = "vertical")            #


ggplot(clusterings_nophylo, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_y_continuous('Total within-cluster sum of squares') +
  scale_x_continuous('K', breaks = seq(1,20, by = 1))




###################################################################################################
# data + labels

# data points used for clustering and cluster_profiles
labelled_nophylo <- assignments_nophylo %>%
  filter(k == 4) %>%
  select(-c(group2,
            kclust,
            tidied,
            glanced,
            k,
            .cluster)) 



# generate permutations for all columns (except 1 - the cluster profile label)
kclust_permutations <- list()

for (i in 2:ncol(labelled_nophylo)){
  kclust_permutations[[i]] <- labelled_nophylo %>%
    permutations(permute = all_of(i), times = 10)
}

names(kclust_permutations) <- colnames(labelled_nophylo)

# bind together permutaed data in a single dataframe, and nest
kclust_permutations %<>% 
  bind_rows(., .id = 'permuted_var') %>%
  unite(id, permuted_var, id)  %>%
  mutate(data = map(splits, ~ analysis(.x))) %>%
  select(-splits) %>%
  unnest(data) 



# execute kmeans on nested data (ie one run of kmeans per permuted dataset)
permuted_nophylo <- kclust_permutations %>%
  #select(-cluster_profile) %>%
  nest(., .by = id) %>%
  mutate(
    kclust = map(data,  ~select(.x , -cluster_profile) %>% kmeans(3)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map2(kclust, data, augment)
  ) 
  

# extract results
permuted_clusters_nophylo <- permuted_nophylo  %>%
  unnest(cols = c(tidied)) 

permuted_assignments_nophylo <- permuted_nophylo  %>% 
  unnest(cols = c(augmented)) 
 
permuted_clusterings_nophylo <- permuted_nophylo  %>%
  unnest(cols = c(glanced))


# was cluster assignment correct according to baseline kmeans?
lookup_clusters <- assignments_nophylo %>%
  filter(k == 3) %>%
  select(c(cluster_profile,
           .cluster)) %>%
  rename(original_cluster = .cluster)


test <- permuted_assignments_nophylo %>%
  select(c(cluster_profile,
           .cluster,
           id)) %>%
  left_join(lookup, by = join_by(cluster_profile)) %>%
  mutate(across(ends_with('cluster'), .fns = ~as.integer(.x))) %>%
  mutate(is_correct = case_when(.cluster == original_cluster ~ 1,
                                .cluster != original_cluster ~ 0)) %>%
  separate_wider_delim(id, delim = regex('_(?!.*_)'), 
                       names = c('var', 'permutation')) %>%
  summarise(prop_correct = mean(is_correct),
            .by = c(var, permutation))





# 
########

kclusts <- tibble(k = 1:20) %>%
  mutate(
    kclust = map(k, ~kmeans(kclust_updated_data %>% select(-group2), .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, kclust_updated_data)
  )

clusters <- 
  kclusts %>%
  unnest(cols = c(tidied))

assignments <- 
  kclusts %>% 
  unnest(cols = c(augmented))

clusterings <- 
  kclusts %>%
  unnest(cols = c(glanced))

assignments %>%
  filter(k %in% seq(from = 3, to = 6)) %>%
  select(-c(kclust, tidied, glanced)) %>%
  group_split(k, .keep = TRUE) %>%
  map(function(df) {
    # 2. Apply PCA to each group separately
    df %>%
      recipe(~ .) %>%
      step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
      prep() %>%
      bake(NULL)
  }) %>%
  bind_rows() %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = .cluster, shape = group2), size = 4, alpha = 0.8) + 
  theme_minimal(base_size = 20) +
  facet_wrap(~ k)

clustersphylo_plot <- assignments %>%
  filter(k ==3) %>%
  select(-c(kclust, tidied, glanced)) %>%
  recipe(~ .) %>%
  step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
  prep() %>%
  bake(NULL)%>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = .cluster, shape = group2), alpha = 0.8) + 
  scale_shape('Pre-Assigned Classification') + 
  scale_color_discrete('K-Means Clusters') +
  theme_minimal() +
  theme(legend.position = 'bottom') +               # Stack legends vertically
  guides(
    color = guide_legend(order = 1),             # Ensure color legend is first
    shape = guide_legend(order = 2)              # Shape legend is second
  )



assignments %>%
  filter(k == 3) %>%
  select(.cluster, group2) %>%
  summarise(n = n(), .by = c(.cluster, group2)) %>%
  pivot_wider(names_from = 'group2', values_from = 'n')

# Elbow plot
ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  scale_y_continuous('Total within-cluster sum of squares') +
  scale_x_continuous('K', breaks = seq(1,20, by = 1))


##### before/after plots

library(cowplot)
legend <- get_legend(clustersphylo_plot)
figure_upper <- plot_grid(clustersnophylo_plot + theme(legend.position = 'none'),
                          clustersphylo_plot + theme(legend.position = 'none'), 
                          align = 'hv',
                          labels = 'AUTO')
figure <- plot_grid(figure_upper, legend , nrow = 2, rel_heights = c(0.9,0.1))

####################################################################################################
####################################################################################################