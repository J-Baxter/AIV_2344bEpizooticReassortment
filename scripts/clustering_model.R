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
library(rsample)
library(mclust)


########################################### IMPORT DATA ############################################

combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

meta <- read_csv('./2024-09-09_meta.csv')

key <- meta %>%
  dplyr::select(cluster_profile,
                cluster_label) %>%
  distinct() %>%
  drop_na()

reassortant_offspring <- read_csv('./reassortant_offspring.csv') %>% 
  dplyr::select(-cluster_class) %>%
  rename(cluster_label = name,
         offspring = importance) %>%
  left_join(key)

# import colour schemes
host_colour <- read_csv('./colour_schemes/hostType_cols.csv')
region_colour <- read_csv('./colour_schemes/regionType_cols.csv')
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
  
  left_join(reassortant_offspring ,
            by = join_by(cluster_profile)) %>%
  
  # select variables of interest
  # must include persistence time, diffusion coefficient, evolutionary rate, species jumps
  select(
    segment,
    cluster_profile,
    weighted_diff_coeff,
    evoRate,
    persist.time,
    count_cross_species,
    count_to_mammal,
    offspring
  ) %>%
  
  
  # Substitute NA values in diffusion coefficient with 0
  mutate(across(where(is.double), .fns = ~ replace_na(.x, 0))) %>%
  rename_with(~gsub('-', '_', .x)) %>%
  
  
  # Model pre-processing
  group_by(segment, cluster_profile) %>%
  slice_sample(n=1) %>%
  ungroup() %>% 
  filter(segment %in% c('ha', 'pb2')) %>%
  
  # pivot wider so one row/reassortant
  pivot_wider(names_from = segment,
              values_from = where(is.double)) %>%
  
  # start tidymodels
  recipe(~ .) %>% 
  
  #drop na
  step_naomit(everything()) %>%
  
  # normalise numeric vectors
  step_zv(everything()) %>%
  step_normalize(all_numeric()) %>% 
  
  # create dummy variables for categorical predictors
  #step_dummy(starts_with('collection_regionname'), one_hot = TRUE) %>%
  #step_dummy(starts_with('host_simplifiedhost'), one_hot = TRUE) %>%
  
  # run tidymodels recipe
  prep() %>%
  bake(NULL) %>%
  
  
  
  # exclude columns not to be used
  left_join(kclust_nophylo_data %>% select(-c(starts_with('group'), Length_Between_First_Last_Sample, number_Conti)), by = join_by(cluster_profile)) %>%
  select(-c(starts_with('host_simplifiedhost'), 
            starts_with(c('collection_regionname'))))


####################################### SUMMARY DATA KCLUST ########################################
set.seed(4472)
kclust_nophylo <- tibble(k = 1:20) %>%
  mutate(
    kclust = purrr::map(k, ~kmeans(kclust_nophylo_data %>% select(-c(group2,cluster_profile)), .x)),
    tidied = purrr::map(kclust, tidy),
    glanced = purrr::map(kclust, glance),
    augmented = purrr::map(kclust, augment, kclust_nophylo_data)
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


####################################### SUMMARY DATA PERMUTATIONS ########################################

# data + labels

# data points used for clustering and cluster_profiles
labelled_nophylo <- assignments_nophylo %>%
  filter(k == 3) %>%
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
  mutate(data = purrr::map(splits, ~ analysis(.x))) %>%
  select(-splits) %>%
  unnest(data) 



# execute kmeans on nested data (ie one run of kmeans per permuted dataset)
permuted_nophylo <- kclust_permutations %>%
  #select(-cluster_profile) %>%
  nest(., .by = id) %>%
  mutate(
    kclust = purrr::map(data,  ~select(.x , -cluster_profile) %>% kmeans(3)),
    tidied = purrr::map(kclust, tidy),
    glanced = purrr::map(kclust, glance),
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
  filter(k == 4) %>%
  select(c(cluster_profile,
           .cluster)) %>%
  rename(original_cluster = .cluster)


####################################### SUMMARY DATA ARI ########################################
# Ari - adjusted rand index
# The Rand Index computes a similarity measure between two clusterings by considering all pairs 
# of samples and counting pairs that are assigned in the same or different clusters in the predicted 
# and true clusterings. The raw RI score is then “adjusted for chance” into the ARI score using the
# following scheme: ARI = (RI - Expected_RI) / (max(RI) - Expected_RI)

ari <- permuted_assignments_nophylo %>%
  select(c(cluster_profile,
           .cluster,
           id)) %>%
  left_join(lookup_clusters, 
            by = join_by(cluster_profile)) %>%
  mutate(across(ends_with('cluster'), .fns = ~as.integer(.x))) %>%
  separate_wider_delim(id, delim = regex('_(?!.*_)'), 
                       names = c('var', 'permutation')) %>%
  group_by(var, permutation) %>%
  mutate(ari = mclust::adjustedRandIndex(original_cluster, .cluster) )%>%
  summarise(ari = mean(ari)) %>%
  ungroup()

ari_summary <- ari %>%
  mutate(importance = 1 - ari) %>% # Higher impact means lower ARI score
  summarise(mean_importance = mean(importance),
            lower_ci = mean(importance) - 1.96 * sd(importance) / sqrt(n()),  # Lower 95% CI
            upper_ci = mean(importance) + 1.96 * sd(importance) / sqrt(n()),  # Upper 95% CI
            .by = c(var))




#ggplot(clusterings_nophylo, aes(k, tot.withinss)) +
# geom_line() +
# geom_point() +
#theme_minimal() +
#scale_y_continuous('Total within-cluster sum of squares') +
#scale_x_continuous('K', breaks = seq(1,20, by = 1))

####################################### SUMMARY DATA PLOTS ########################################
riskgroup_colour <- read_csv('./colour_schemes/riskgroup_cols.csv') %>%
  mutate(kclust = c(3,1, 2),
         group2 = c('minor', 'moderate', 'major'))


# Elbow Plot for summary data 
plt_1a <- ggplot(clusterings_nophylo, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 8) +
  scale_y_continuous('Total within-cluster sum of squares') +
  scale_x_continuous('K', breaks = seq(0,20, by = 5)) + 
  geom_vline(xintercept = 3, colour= 'darkgreen', linetype = 'dashed')


plt_1b <- assignments_nophylo %>%
  filter(k == 3) %>%
  select(-c(kclust, tidied, glanced)) %>%
  recipe(~ .) %>%
  step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
  prep() %>%
  bake(NULL)%>%
  left_join(meta %>% select(cluster_profile, cluster_label) %>% distinct() %>% drop_na()) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = .cluster, alpha = ifelse(cluster_label %in% c('H5N1/2022/R7_NAmerica', 
                                                                        'H5N1/2020/R1_Europe', 
                                                                        'H5N8/2019/R7_Africa', 
                                                                        'H5N1/2021/R3_Europe', 
                                                                        'H5N1/2022/R12_Europe',
                                                                        'H5N1/2021/R1_Europe',
                                                                        'H5N1/2023/R29_NAmerica'), '1', '0')), size = 2) + 
  geom_text(aes(label=ifelse(cluster_label %in% c('H5N1/2022/R7_NAmerica', 
                                                  'H5N1/2020/R1_Europe', 
                                                  'H5N8/2019/R7_Africa', 
                                                  'H5N1/2021/R3_Europe', 
                                                  'H5N1/2022/R12_Europe',
                                                  'H5N1/2021/R1_Europe',
                                                  'H5N1/2023/R29_NAmerica'), gsub('_.*', '', cluster_label),""), 
                colour = .cluster), hjust = 0, nudge_x = 0.15, nudge_y = -0.1, size = 2) +

  scale_colour_manual('K-Means Clusters',
                      values = c('3' = '#FF0000',
                      '2' = '#2ca02c',
                      '1' = '#1f77b4')) +
  scale_alpha_manual(values = c('1' = 1, '0' = 0.3)) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))
        



plt_1c <- ggplot(ari_summary) +
  geom_point(aes(x = var, y = mean_importance)) +
  geom_linerange(aes(x = var, ymin = lower_ci, ymax = upper_ci)) + 
  
  scale_x_discrete('Variable',
                   labels = c('host_richness' = 'Host richness',
                              'if_mammal' = 'Mammals' ,
                              'Length_Between_First_Last_Sample' = 'Time between samples',
                              'max_distance_km' = 'Maximum distance',
                              'Num_Sequence' = 'Number of sequences',
                              'number_Conti' = 'Number of Continents') %>%
                     str_wrap(., width = 15)) +
  scale_y_continuous(
    'Variable Importance') +
  coord_cartesian(ylim = c(0,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


plt_1 <- cowplot::plot_grid( plt_1b, plt_1c, ncol = 1, align = 'hv', axis = 'tb', labels = 'AUTO', label_size = 10)

ggsave(
       '~/Downloads/flu_plots/summarydata_clustering.jpeg', 
       plt_1,
       height = 20,
       width = 15,
       units =  "cm",
       dpi = 360)

####################################### PHYLO + SUMMARY DATA KCLUST ########################################

# Now repeat for phylo data

kclusts <- tibble(k = 1:20) %>%
  mutate(
    kclust = purrr::map(k, ~kmeans(kclust_updated_data %>% select(-cluster_profile) %>%drop_na(), .x)),
    tidied = purrr::map(kclust, tidy),
    glanced = purrr::map(kclust, glance),
    augmented = purrr::map(kclust, augment, kclust_updated_data %>%drop_na()))

clusters <- kclusts %>%
  unnest(cols = c(tidied))

assignments <- kclusts %>% 
  unnest(cols = c(augmented))

clusterings <-  kclusts %>%
  unnest(cols = c(glanced))

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 8) +
  scale_y_continuous('Total within-cluster sum of squares') +
  scale_x_continuous('K', breaks = seq(0,20, by = 5)) + 
  geom_vline(xintercept = 4, colour= 'darkgreen', linetype = 'dashed')

####################################### PHYLO + SUMMARY DATA PERMUTATIONS #######################################
# data points used for clustering and cluster_profiles
labelled_phylo <- assignments %>%
  filter(k == 3) %>%
  select(-c( kclust,
            tidied,
            glanced,
            k,
            .cluster)) 



# generate permutations for all columns (except 1 - the cluster profile label)
kclust_permutations <- list()

for (i in 2:ncol(labelled_phylo)){
  kclust_permutations[[i]] <- labelled_phylo %>%
    permutations(permute = all_of(i), times = 10)
}

names(kclust_permutations) <- colnames(labelled_phylo)

# bind together permutaed data in a single dataframe, and nest
kclust_permutations %<>% 
  bind_rows(., .id = 'permuted_var') %>%
  unite(id, permuted_var, id)  %>%
  mutate(data = purrr::map(splits, ~ analysis(.x))) %>%
  select(-splits) %>%
  unnest(data) 



# execute kmeans on nested data (ie one run of kmeans per permuted dataset)
permuted_phylo <- kclust_permutations %>%
  #select(-cluster_profile) %>%
  nest(., .by = id) %>%
  mutate(
    kclust = purrr::map(data,  ~select(.x , -cluster_profile) %>% kmeans(3)),
    tidied = purrr::map(kclust, tidy),
    glanced = purrr::map(kclust, glance),
    augmented = map2(kclust, data, augment)
  ) 


# extract results
permuted_clusters_phylo <- permuted_phylo  %>%
  unnest(cols = c(tidied)) 

permuted_assignments_phylo <- permuted_phylo  %>% 
  unnest(cols = c(augmented)) 

permuted_clusterings_phylo <- permuted_phylo  %>%
  unnest(cols = c(glanced))



####################################### PHYLO + SUMMARY DATA ARI ########################################
# was cluster assignment correct according to baseline kmeans?
lookup_clusters <- assignments %>%
  filter(k == 3) %>%
  select(c(cluster_profile,
           .cluster)) %>%
  rename(original_cluster = .cluster)


ari_phylo <- permuted_assignments_phylo %>%
  select(c(cluster_profile,
           .cluster,
           id)) %>%
  left_join(lookup_clusters, 
            by = join_by(cluster_profile)) %>%
  mutate(across(ends_with('cluster'), .fns = ~as.integer(.x))) %>%
  separate_wider_delim(id, delim = regex('_(?!.*_)'), 
                       names = c('var', 'permutation')) %>%
  group_by(var, permutation) %>%
  mutate(ari = mclust::adjustedRandIndex(original_cluster, .cluster) )%>%
  summarise(ari = mean(ari)) %>%
  ungroup()

ari_phylo_summary <- ari_phylo %>%
  mutate(importance = 1 - ari) %>% # Higher impact means lower ARI score
  summarise(mean_importance = mean(importance),
            lower_ci = mean(importance) - 1.96 * sd(importance) / sqrt(n()),  # Lower 95% CI
            upper_ci = mean(importance) + 1.96 * sd(importance) / sqrt(n()),  # Upper 95% CI
            .by = c(var))


####################################### PHYLO + SUMMARY DATA PLOTS ########################################
riskgroup_colour <- read_csv('./colour_schemes/riskgroup_cols.csv') %>%
  mutate(kclust = c(3,1,2),
         group2 = c( 'minor', 'major','dominant'))

####################################################################################################
# Elbow Plot for summary data 
plt_2a <- ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 8) +
  scale_y_continuous('Total within-cluster sum of squares') +
  scale_x_continuous('K', breaks = seq(0,20, by = 5)) + 
  geom_vline(xintercept = 3, colour= 'darkgreen', linetype = 'dashed')


marked_clusters <- reassortant_offspring %>% 
  filter(offspring >= 5) %>%
  pull(cluster_label)

marked_clusters_2 <- reassortant_ancestral_changes %>% filter(cluster_class == 'major') %>% pull(cluster_label)
plt_2b <-  
  assignments %>%
  filter(k == 3) %>%
  select(-c(kclust, tidied, glanced)) %>%
  recipe(~ .) %>%
  step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
  prep() %>%
  bake(NULL) %>%
  left_join(meta %>% select(cluster_profile, cluster_label) %>% distinct() %>% drop_na()) %>%
  ggplot(., aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = .cluster, alpha = ifelse(cluster_label %in% marked_clusters_2, '1', '0')), size = 2) + 
  geom_text(aes(label=ifelse(cluster_label %in% marked_clusters_2, gsub('_.*', '', cluster_label),""), 
                colour = .cluster), hjust = 1.1, nudge_y = -0.1, size = 2) +
  
  scale_colour_brewer('K-Means Clusters',
                      palette = 'Set1')+
  scale_alpha_manual(values = c('1' = 1, '0' = 0.3)) + 
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


cluster.stats(d = NULL, clustering, al.clustering = NULL)



plt_2c <- ari_phylo_summary %>%
  filter(!var %in% ari_summary$var) %>%
  ggplot()+
  geom_point(aes(x = var, y = mean_importance)) +
  geom_linerange(aes(x = var, ymin = lower_ci, ymax = upper_ci)) + 
  
  scale_x_discrete('Variable',
                   labels = c('persist.time_ha' = 'HA Persistence Time' ,
                              'persist.time_pb2' = 'PB2 Persistence Time' ,
                              'persist.time_nx' = 'Nx Persistence Time',
                              'weighted_diff_coeff_ha' = 'HA diffusion coefficient',
                              'weighted_diff_coeff_pb2' = 'PB2 diffusion coefficient',
                              'weighted_diff_coeff_nx' = 'Nx diffusion coefficient',
                              'evoRate_pb2' = 'PB2 evo rate',
                              'evoRate_ha' = 'HA evo rate',
                              'evoRate_nx' = 'Nx evo rate',
                              'count_cross_species_pb2' = 'PB2 species jumps',
                              'count_cross_species_ha' = 'HA species jumps',
                              'count_cross_species_nx' = 'Nx species jumps',
                              'count_to_mammal_pb2' = 'PB2 mammal jumps',
                              'count_to_mammal_ha' = 'HA mammal jumps',
                              'count_to_mammal_nx' = 'Nx mammal jumps') %>%
                     str_wrap(., width = 15)
  ) +
  scale_y_continuous(
    'Variable Importance') +
  coord_cartesian(ylim = c(0,1)) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8))


plt_2 <- cowplot::plot_grid( plt_2b, plt_2c, ncol = 1, align = 'hv', axis = 'tb', labels = 'AUTO', label_size = 10)

ggsave(
  '~/Downloads/flu_plots/phylodata_clustering.jpeg', 
  plt_2,
  height = 20,
  width = 15,
  units =  "cm",
  dpi = 360)



####################################### ALL DATA CONTINGENCY TBL ########################################
assignments %>%
  filter(k == 4) %>%
  select(c(cluster_profile, .cluster)) %>%
  left_join(assignments_nophylo %>%
              filter(k == 3) %>%
              select(c(cluster_profile, .cluster)), by = join_by(cluster_profile)) %>%
  select(-cluster_profile) %>%
  ftable()


####################################### Silhouettes ########################################
# Calculate the mean intra-cluster distance (a) and the mean nearest-cluster 
# distance (b) for each sample. The Silhouette Coefficient for a sample is (b - a) / max(a, b)

sil_nophlyo <- cluster::silhouette(assignments_nophylo %>%
                             filter(k == 4) %>% pull(.cluster) %>%
                             as.integer(),
                           assignments_nophylo %>%
                             filter(k == 4) %>%
                             select(-c(k, kclust, tidied, glanced, cluster_profile, .cluster,group2)) %>% as.matrix() %>% dist()) %>%
  .[, 1:3] %>%
  as_tibble() %>%
  mutate(data = 'Summary Data') %>%
  mutate(cluster = case_when(cluster == 1 ~ 'minor', 
                             cluster == 2 ~ 'major a',
                             cluster == 3 ~ 'major b',
                             cluster == 4 ~ 'dominant',
                             ))




sil <- cluster::silhouette(assignments %>%
                    filter(k == 4) %>% pull(.cluster) %>%
                                              as.integer(),
                  assignments %>%
                    filter(k == 4) %>%
                    select(-c(k, kclust, tidied, glanced, cluster_profile, .cluster)) %>% as.matrix() %>% dist()
                    ) %>%
  .[, 1:3] %>%
  as_tibble() %>%
  mutate(data = 'All Data') %>%
  mutate(cluster = case_when(cluster == 3 ~ 'minor',
                             cluster == 4 ~ 'major b',
                             cluster == 2 ~  'major a',
                             cluster == 1 ~ 'minor'))


riskgroup_colour <- read_csv('./colour_schemes/riskgroup_cols.csv') %>%
  mutate(kclust = c(3,1,2),
         group2 = c( 'minor', 'major','dominant'))

#a negative silhouette score for a point does not necessarily mean a "wrong" assignment, it 
# just means that the point is on average closer to points in another cluster than to points 
# in its own cluster

bind_rows(sil,
          sil_nophlyo) %>%
  ggplot(aes(x = sil_width, fill = cluster)) +
  geom_density( colour = NA, alpha = 0.8) +
  scale_fill_brewer(palette = 'Dark2') +
  facet_grid(cols = vars(data)) +
  scale_x_continuous('Silhouette Score')+
  scale_y_continuous('Probability Density') + 
  theme_minimal(base_size = 18) +
  theme(legend.position = 'bottom')

### facet pca plot for different clusters
#assignments %>%
  #filter(k %in% seq(from = 3, to = 6)) %>%
  #select(-c(kclust, tidied, glanced)) %>%
  #group_split(k, .keep = TRUE) %>%
  #map(function(df) {
    # 2. Apply PCA to each group separately
   # df %>%
     # recipe(~ .) %>%
     # step_pca(all_numeric(), -k, num_comp = 2) %>%  # Perform PCA (e.g., 2 components)
     # prep() %>%
     # bake(NULL)
 # }) %>%
 # bind_rows() %>%
 # ggplot(., aes(x = PC1, y = PC2)) +
  #geom_point(aes(colour = .cluster, shape = group2), size = 4, alpha = 0.8) + 
  #theme_minimal(base_size = 20) +
 # facet_wrap(~ k)



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