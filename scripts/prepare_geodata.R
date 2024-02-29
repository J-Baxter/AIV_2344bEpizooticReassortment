# prepare geodata

country <- read_sf('./spatial_data/gadm_410-levels.gpkg', layer = 'ADM_0')
subdivision_one <- read_sf('./spatial_data/gadm_410-levels.gpkg', layer = 'ADM_1')

country_centroid <- country %>% 
  mutate(is.valid = st_is_valid(geom)) %>%
  mutate(geom = case_when(is.valid != TRUE ~ st_make_valid(geom), .default = geom)) %>%
  mutate(centroid = st_centroid(geom, of_largest_polgon = TRUE))

saveRDS(country_centroid, 'temp1.rdata')

subdiv1_centroid <- subdivision_one %>% 
  mutate(is.valid = st_is_valid(geom)) %>%
  mutate(geom = case_when(is.valid != TRUE ~ st_make_valid(geom), .default = geom)) %>%
  mutate(centroid = st_centroid(geom, of_largest_polgon = TRUE))
saveRDS(subdiv1_centroid, 'temp2.rdata')

subdiv1_withcountry <- subdiv1_centroid %>% as_tibble() %>%
  left_join(as_tibble(country_centroid), by = join_by(COUNTRY, GID_0)) %>%
  rename(adm0_centroid = centroid.y,
         adm0_geom = geom.y,
         adm1_centroid = centroid.x,
         adm1_geom = geom.x)

country_centroid <- country_centroid %>%
  as_tibble() %>%
  rename(adm0_centroid = centroid,
         adm0_geom = geom)


all_geom <- bind_rows(subdiv1_withcountry, country_centroid)
saveRDS(all_geom, 'temp3.rdata')

necessary_geom <- all_geom %>%
  select(c(GID_0, 
           COUNTRY,
           NAME_1,
           VARNAME_1,
           HASC_1,
           ISO_1,
           ends_with('centroid'))) %>%
  mutate(across(ends_with('centroid'), .fns = ~ as.character(.x))) %>%
  mutate(across(ends_with('centroid'), .fns = ~ str_remove_all(.x, "[c(),]"))) %>%
  separate(adm0_centroid, into = c("adm0_long", "adm0_lat"), sep = " ") %>%
  separate(adm1_centroid, into = c("adm1_long", "adm1_lat"), sep = " ") %>%
  mutate(across(where(is.character), .fns = ~ tolower(.x))) %>%
  rename_with(., ~ tolower(.x))

write_csv(necessary_geom, './spatial_data/gdam_geodata.csv')

