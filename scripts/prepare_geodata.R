# prepare geodata
library(sf)
country <- read_sf('./spatial_data/gadm_410-levels.gpkg', layer = 'ADM_0')
subdivision_one <- read_sf('./spatial_data/gadm_410-levels.gpkg', layer = 'ADM_1')
subdivision_two <- read_sf('./spatial_data/gadm_410-levels.gpkg', layer = 'ADM_2')

country_centroid <- country %>% 
  mutate(is.valid = st_is_valid(geom)) %>%
  mutate(geom = case_when(is.valid != TRUE ~ st_make_valid(geom), .default = geom)) %>%
  mutate(centroid = st_centroid(geom, of_largest_polgon = TRUE))

#saveRDS(country_centroid, 'temp1.rdata')

subdiv1_centroid <- subdivision_one %>% 
  mutate(is.valid = st_is_valid(geom)) %>%
  mutate(geom = case_when(is.valid != TRUE ~ st_make_valid(geom), .default = geom)) %>%
  mutate(centroid = st_centroid(geom, of_largest_polgon = TRUE))

subdiv2_centroid <- subdivision_two %>% 
  mutate(is.valid = st_is_valid(geom)) %>%
  mutate(geom = case_when(is.valid != TRUE ~ st_make_valid(geom), .default = geom)) %>%
  rowwise() %>% 
  mutate(centroid = possibly(st_centroid, otherwise = NA)(geom, of_largest_polgon = TRUE)) %>% 
  as_tibble()

#saveRDS(subdiv1_centroid, 'temp2.rdata')
#saveRDS(subdiv1_centroid, 'temp2.rdata')

subdiv1_withcountry <- subdiv1_centroid %>% as_tibble() %>%
  left_join(as_tibble(country_centroid), by = join_by(COUNTRY, GID_0)) %>%
  rename(adm0_centroid = centroid.y,
         adm0_geom = geom.y,
         adm1_centroid = centroid.x,
         adm1_geom = geom.x) subdiv1_withcountry %>%
 unite(., 'match', COUNTRY, NAME_1, sep = '_', remove = FALSE)


subdiv2_withcountry <- subdiv2_centroid %>% as_tibble() %>%
  left_join(as_tibble(subdiv1_withcountry), by = join_by(COUNTRY, GID_0, GID_1, NAME_1, NL_NAME_1)) %>%
  rename(adm2_centroid = centroid,
         adm2_geom = geom) %>%
  unite(., 'match', COUNTRY, NAME_1, NAME_2, sep = '_', remove = FALSE)


country_centroid <- country_centroid %>%
  as_tibble() %>%
  rename(adm0_centroid = centroid,
         adm0_geom = geom) %>%
  mutate(match = COUNTRY)

#load('temp2.rdata')

all_geom <- bind_rows(subdiv2_withcountry, subdiv1_withcountry, country_centroid)
saveRDS(all_geom, 'temp3.rdata')

necessary_geom <- all_geom %>%
  select(c(GID_0, 
           COUNTRY,
           NAME_1,
           HASC_1,
           ISO_1,
           NAME_2,
           HASC_2,
           ends_with('centroid'),
           match)) %>%
  mutate(across(ends_with('centroid'), .fns = ~ as.character(.x))) %>%
  mutate(across(ends_with('centroid'), .fns = ~ str_remove_all(.x, "[c(),]"))) %>%
  separate(adm0_centroid, into = c("adm0_long", "adm0_lat"), sep = " ") %>%
  separate(adm1_centroid, into = c("adm1_long", "adm1_lat"), sep = " ") %>%
  separate(adm2_centroid, into = c("adm2_long", "adm2_lat"), sep = " ") %>%
  rowwise() %>%
  mutate(across(starts_with('adm2'), .fns = ~ case_when(is.list(.x) ~ NA, .default = .x))) %>%
  as_tibble() %>%
  mutate(across(where(is.character), .fns = ~ case_when(!grepl('^NA$', .x) ~ tolower(.x), .default = NA))) %>%
  rename_with(., ~ tolower(.x)) %>%
  left_join(., region) %>% # from annotated_geodata.csv %>% region <- old_geo %>% select(c(gid_0, country, region)) %>% mutate(region = case_when(country == 'russia' ~ 'eastern europe', .default = region)) %>% distinct()
  relocate(region,
           gid_0, 
           country,
           adm0_lat,
           adm0_long,
           name_1,
           hasc_1,
           iso_1,
           adm1_lat,
           adm1_long,
           name_2,
           hasc_2,
           adm2_lat,
           adm2_long,
           match) 

issues <- necessary_geom%>%filter(is.na(adm2_lat))
write_csv(necessary_geom, './spatial_data/updated_geodata.csv')

