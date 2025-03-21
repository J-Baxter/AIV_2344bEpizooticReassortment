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
         time_since_previous)


  
ggplot(class_data) + geom_bar(aes(x = class, fill = previous_class), position = 'stack') + facet_wrap( ~ collection_regionname)
ggplot(class_data) + geom_point(aes(x = class, y = time_since_previous), position = 'jitter') + facet_wrap( ~ collection_regionname)


class_formula <-bf(class ~  collection_regionname + previous_class + time_since_previous)

class_formula_priors <- get_prior(class_formula,
                                  data = class_data,
                                  family = categorical(link ='logit')) 



diffusionmodel1_priors$prior[1:8] <- "normal(0,5)"

basic_model <- brm(
  class_formula
  data = class_data,
  family = categorical(link ='logit'),
  chains = 2, iter = 10000, cores = 2
)

  