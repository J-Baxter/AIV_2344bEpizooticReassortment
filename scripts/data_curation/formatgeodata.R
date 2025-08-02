# format geographical data to csv
library(tidyverse)
library(V8)
library(countrycode)
# dataset from https://github.com/oodavid/iso-3166-2/blob/master/iso_3166_2.js
cx <- v8()
cx$source("iso_3166_2.js")
data <- cx$get("iso_3166_2") %>%
  bind_rows() %>%
  
  # Remove higher-level labels
  filter(is.na(parent) | parent != "EARTH") %>%
  filter(is.na(division) | division != "planet") %>%
  
  # Extract county codes from subdivision codes
  mutate(countrycode = case_when(
    grepl("-| ", parent) ~ sub("-.*| .*", "", parent),
    is.na(parent) ~ code,
    .default = parent
  )) %>%
  
  # map subdivisions to countries
  left_join(cx$get("iso_3166_2") %>%
    bind_rows(),
    by = join_by(countrycode == code)) %>%
  
  mutate(across(contains("x"), .fns = ~ case_when(
    grepl("country", division.x) ~ NA,
    .default = .x
  ))) %>%
  
  mutate(code = case_when(
    str_count(code) == 2 ~ NA,
    .default = code
  )) %>%
  
  # map countries to regions
  left_join(countrycode::codelist %>% 
              select(c(iso2c, region23)), 
            by = join_by(countrycode == iso2c)) %>%
  
  # interim re-formatting
  select(-contains("division")) %>%
  rename(
    iso2_country = countrycode,
    name_country = name.y,
    lat_country = lat.y,
    long_country = lng.y,
    region = region23
  ) %>%
  
  # pivot wider subdivisions
  group_by(iso2_country) %>%
  
  mutate(subdivision_tier = case_when(
    str_count(parent.x) == 2 ~ "subdiv1",
    str_count(parent.x) > 2 ~ "subdiv2",
    .default = NA
  )) %>%
  
  mutate(code.x = code) %>%
  
  pivot_wider(
    names_from = subdivision_tier,
    values_from = contains("x"),
    names_sep = "_"
  ) %>%
  ungroup() %>%
  
  # remove deprecated columns
  select(-c(contains("_NA"), code)) %>%
  
  # tidy renaming and reformatting
  rename_with(~ gsub(".x", "", .x)) %>%
  rename_with(~ gsub("lng_", "long_", .x)) %>%
  relocate(region,
           contains("country"), 
           contains("subdiv1"), 
           contains("subdiv2")) %>%
  
  # map subdivisions
  left_join(cx$get("iso_3166_2") %>% 
              bind_rows(), 
            by = join_by(parent_subdiv2 == code)) %>%
  
  mutate(name_subdiv1 = coalesce(name_subdiv1, name)) %>%
  mutate(code_subdiv1 = coalesce(code_subdiv1, parent_subdiv2)) %>%
  mutate(lat_subdiv1 = coalesce(lat_subdiv1, lat)) %>%
  mutate(long_subdiv1 = coalesce(long_subdiv1, lng)) %>%
  
  # drop deprecated rows and rename
  select(-c(contains("parent"), lat, lng, division, name)) %>%
  rename_with(~ gsub("code", "iso2", .x))  %>%
  mutate(row_id = row_number())
  
  # generate query column
query <- data %>%
  mutate(query = case_when(
    is.na(name_subdiv1) &  is.na(name_subdiv2) ~ name_country,
    !is.na(name_subdiv1) &  is.na(name_subdiv2) ~ name_subdiv1,
    !is.na(name_subdiv1) &  !is.na(name_subdiv2) ~ name_subdiv2)) %>%
  distinct() %>%
  mutate(across(-contains('iso2'), .fns = ~ tolower(.x)))
  
#t_geo <- apply(query,  1,
 # function(x) paste0('atlasTable[', '"', x['query'], '"]', '= {', '"', x['region'], '", ', '"', 
                   #  x['iso2_country'], '", ', '"', x['name_country'], '", ', '"', 
                    # str_trim(x['lat_country']), '", ', '"', str_trim(x['long_country']), '", ', '"',
                    # x['iso2_subdiv1'], '", ', '"', x['name_subdiv1'], '", ', '"', 
                   #  str_trim(x['lat_subdiv1']), '", ', '"', str_trim(x['long_subdiv1']),  '", ', '"',
                    # x['iso2_subdiv2'], '", ', '"', x['name_subdiv2'], '", ', '"',
                    # str_trim(x['lat_subdiv2']), '", ', '"', str_trim(x['long_subdiv2']), '" };'))

#write_lines(t_geo, 'test_geo.txt')  


write_csv(data, "geographicalmetadata.csv")

# END #
