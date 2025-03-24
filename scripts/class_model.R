# class model

# variables required
# - class of reassortant /
# - immediately preceeding reassortant /
# - time since last ressortant /
# - continent /
# (- diversity)
# - proxy incidence
# - number of genomes


combined_data <- read_csv('./2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')
summary_data <- read_csv('./2024Aug18/treedata_extractions/summary_reassortant_metadata_20240904.csv') %>%
  select(-c(cluster_label,
            clade)) 

final_clusters <- read_csv('./final_clusters_2025Mar21.csv')

class_data <- combined_data %>% 
  select(-c(group2, col2)) %>%
  left_join(final_clusters) %>%
  
  select(cluster_profile, 
          segment,
          TMRCA,
          collection_regionname,
          class) %>%
  
  mutate(collection_regionname = case_when(grepl('europe', collection_regionname) ~ 'europe',
                                           grepl('africa', collection_regionname) ~ 'africa',
                                           grepl('asia', collection_regionname) ~ 'asia',
                                           grepl('(central|northern) america', collection_regionname) ~ 'central & northern america',
                                           .default = collection_regionname
  )) %>%
  
  mutate(TMRCA = date_decimal(TMRCA) %>%
           round_date(unit = 'day')%>%
           as_date(.)) %>%
  
  filter(segment== 'ha') %>%
  arrange(TMRCA) %>%
  group_by(collection_regionname) %>%
  
  # class of most recent prior reassortant
  mutate(previous_class = lag(class),
         previous_class = replace_na(previous_class, 'none')) %>%
  
  # time since last reassortant
  # calculate as interval, then convert to decimal of years
  mutate(time_since_previous = interval(lag(TMRCA), TMRCA) %>%
           as.numeric('years')) %>%
  
  # revert formatting of dates for model
  # TMRCA to decimal, then subtract 2016 for scale
  mutate(tmrca_decimal = decimal_date(TMRCA),
         tmrca_decimal_adjusted = tmrca_decimal - 2016) %>%
  
  # select vars of interest
  select(collection_regionname,
         class, 
         previous_class,
         tmrca_decimal_adjusted,
         time_since_previous) %>%
  
  # exclude the first reassortant from each continent
  drop_na(time_since_previous) %>%
  drop_na(class)


  
ggplot(class_data) + geom_bar(aes(x = class, fill = previous_class), position = 'stack') + facet_wrap( ~ collection_regionname)
ggplot(class_data, aes(x = time_since_previous)) + geom_histogram() + 
  #geom_smooth(method = "lm", se = FALSE) +
  facet_wrap( class ~ collection_regionname)
    


class ~ 0 + collection_regionname + time_since_previous 

time_since_previous + previous_class + (1 | collection_regionname:previous_class)


class_formula_priors <- get_prior(class_formula,
                                  data = class_data,
                                  family = categorical(link ='logit')) 


class_formula_priors$prior <- "normal(0,5)"

class_formula_priors$prior[c(1:8, 10:17)] <- "normal(0,5)"

basic_model <- brm(
  class_formula,
  prior = class_formula_priors,
  data = class_data,
  family = categorical(link ='logit'),
  chains = 2, iter = 10000, cores = 2
)

pp_check(basic_model, ndraws = 100) 


new_data = class_data %>% 
  select(collection_regionname) %>% 
  distinct() %>% 
  mutate(time_since_previous = median(class_data$time_since_previous))

# expectation of the posterior
basic_model %>%
  epred_draws(newdata = new_data) %>%
  #filter(.epred >0.01) %>%
  ggplot(aes(x  = .epred, y = collection_regionname,
             slab_colour = collection_regionname,
             slab_fill = collection_regionname)) +
  stat_halfeye(p_limits = c(0.001, 0.999),
               point_interval = "median_hdi",
              # expand = T,
               #density = 'histogram',
               .width =  0.95,
               normalize = "xy") +
  facet_grid(cols = vars(.category)) 

basic_model %>%
  avg_slopes(newdata = 'balanced')

basic_model %>%
  predicted_draws(newdata = new_data %>% filter(previous_class == 'minor')) %>%
  count(predicted_class = .prediction) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = factor(predicted_class), y = n / sum(n))) +
  geom_col() + 
  facet_grid(cols = vars(collection_regionname))




# transition probabilities
test <- basic_model %>%
  epred_draws(newdata = new_data)

tidy_dag <- test %>%
  filter(previous_reassortant_class != 'none') %>%
  filter(.category != 'none') %>%
  median_hdci() %>%
  rename(name = previous_reassortant_class, 
         to = .category) %>% 
  as_tidy_dagitty() 

tidy_dag %>%
  mutate(name = case_when(name == 'major'~ 'Mj',
                          name == 'moderate' ~ "Md",
                          name == 'minor' ~ 'Mn')) %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_point(aes(alpha = name, colour = collection_regionname)) +
  geom_dag_text(size = 3) +
  geom_dag_edges_fan(aes(edge_width = .epred*2,
                         label = paste0(round(.epred, digits = 3), 
                                        ' (', round(.lower, digits = 3),
                                        '-', round(.upper, digits = 3), ')')), 
                     spread = 5,
                     label_size = 3,
                     label_dodge	= 0.5) + 
  #geom_dag_edges_arc(aes(edge_width = .epred))+ 
  facet_grid(cols = vars(collection_regionname),
             labeller = labeller(collection_regionname= str_to_title)) + 
  scale_alpha_manual(values = c('D' = 1, 'Mj' = 0.7, 'Mn' = 0.5)) +
  scale_colour_manual(values = region_colours) +
  theme_dag(legend.position = 'none', base_size = 8)