################################################################################
## Script Name:        Tree Data Extraction
## Purpose:            Update combined data with persistence, tmrca and evo.rate
##                    estimates from current analysis
## Author:             James Baxter
## Date Created:       2025-08-02
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(beastio)

################################### DATA #######################################
# Read and inspect data
csv_paths <- list.files("./2024Aug18/treedata_extractions",
                        pattern = ".csv",
                        full.names = TRUE
)

csv_names <- gsub(
  ".*treedata_extractions/|.csv$",
  "",
  csv_paths
)


# Import Tree Data Extractions
csv_data <- lapply(csv_paths, read_csv) %>%
  setNames(csv_names)

# Expand to separate objects - this includes those called later:
# >summary_reassortant_metadata_20240904,
# >reassortant_stratifiedpersistence,
# >RegionalTre_mrca_stats_disp_jumps_combined
#>2024-09-20_combined_data
#>domainantTre_mrca_stats_disp_jumps_combined
for (i in 1:length(csv_data)) {
  assign(names(csv_data)[i], as_tibble(csv_data[[i]]))
}

reassortant_log_path <- list.files("./2024Nov18/plain_tree",
                                   pattern = '.log',
                                   full.names = T)


################################### MAIN #######################################
# Main analysis or transformation steps

## Import Evolutionary Rates from Reassortant log files
combined_data <- RegionalTre_mrca_stats_disp_jumps_combined %>%
  # remove dominant cluster profiles
  filter(!Rcode %in% c(
    "3_2_3_1_3_2_1_2",
    "2_1_1_1_1_1_1_1",
    "1_1_1_1_1_1_1_1",
    "2_1_2_1_1_1_1_1",
    "1_1_2_1_1_1_1_1",
    "1_6_2_1_1_1_1_1",
    "1_1_4_1_4_1_1_4",
    "2_6_1_1_6_1_1_1",
    "2_1_6_1_1_4_1_1",
    "7_1_5_2_1_3_1_2",
    "4_3_1_1_2_1_1_3",
    "5_1_1_1_2_1_1_3",
    "5_4_9_1_2_1_1_1"
  )) %>%
  # remove deprecated columns
  select(-c(
    cluster_number,
    collection_countryname,
    
    prob_cluster_number,
    prob_collection_countryname,
    
  )) %>%
  # bind results from dominant trees
  bind_rows(domainantTre_mrca_stats_disp_jumps_combined) %>%
  # bind summary table
  left_join(
    summary_reassortant_metadata_20240904 %>% ###
      select(c(
        cluster_profile,
        col2,
        group2,
        Length_Between_First_Last_Sample
      )),
    by = join_by(Rcode == cluster_profile)
  ) %>%
  
  # split labels to extract segments
  separate_wider_delim(data_name, "_", names = c("segment", "deprecated")) %>%
  mutate(cluster_profile = case_when(deprecated == "11111111A" & Rcode == "1_1_1_1_1_1_1_1" ~ "1_1_1_1_1_1_1_1A",
                                     .default = Rcode
  )) %>%
  select(-c(deprecated, Rcode)) %>%
  mutate(segment = gsub("n[:0-9:]", "nx", segment)) %>%
  relocate(cluster_profile, segment)


# Import Evolutionary Rates from Reassortant log files 
imported_logs <- lapply(reassortant_log_path, readLog) %>%
  lapply(., getHPDMedian) %>%
  lapply(., function(x) as_tibble(x, rownames = 'parameter') %>% filter(parameter %in% c("ucld.mean", "age.root."))) %>%
  setNames(gsub(".*plain_tree/|_subsampled.*", "", reassortant_log_path)) %>%
  bind_rows(.id = 'key') %>%
  mutate(key = gsub("n[:0-9:]", "nx", key),
         parameter = gsub('age\\.root\\.', 'TMRCA', parameter) %>%
           gsub('ucld\\.mean', 'evoRate', .)) %>%
  rename(`95lower` = lower,
         `95upper` = upper) %>%
  pivot_wider(names_from = parameter,
              values_from = c('95lower', 'med', '95upper'),
              names_glue =  "{parameter}_{.value}") %>%
  rename_with(.fn = ~gsub('_med', '', .x) , .cols = contains('med'))


reassortant_stratifiedpersistence %<>%
  rename_with(.fn = ~gsub('persistence_host_|_simplifiedhost', '', .x)) %>%
  pivot_wider(names_from = host,
              values_from = c(median, max, sum),
              names_sep = '_',
              values_fn = median) 


################################### OUTPUT #####################################
# Save output files, plots, or results
# Combine all data
combined_data <- combined_data %>%
  unite(key, segment, cluster_profile, remove = F) %>%
  mutate(key = str_remove_all(key, "(?<=_\\d{1,2})_")) %>%
  
  
  rows_update(.,
              imported_logs %>% dplyr::select(- c(TMRCA_95lower , TMRCA_95upper)),
              by = "key",
              unmatched = "ignore"
  ) %>%
  # join reassortant-stratified host persistence time
  
  left_join(.,
            y = reassortant_stratifiedpersistence) %>%
  select(-contains('+'))


write_csv(combined_data, './2024Aug18/treedata_extractions/2024-09-20_combined_data.csv')


#################################### END #######################################
################################################################################