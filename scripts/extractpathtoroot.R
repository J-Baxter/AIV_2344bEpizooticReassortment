library(treeio)

test_tree <- read.beast('./2024Aug18/reassortant_subsampled_outputs/traits_mcc/ha_11111111_subsampled_traits_mcc.tree')
test_tree_tidy <- as_tibble(test_tree) %>%
  mutate(tipdate =  str_extract(label,  "(?<=\\|)\\d{4}(?![[:lower:]]).*$") %>%
           ymd(.)) %>%
  select(parent, 
         node,
         height_median,
         collection_regionname,
         host_simplifiedhost) %>%
  left_join( x = . ,
             y = .,
             by = c('parent' = 'node'), suffix = c('', '_parent')) %>%
  select(-c(parent_parent)) %>%
  
  rename(fromNode = parent,
         toNode = node,
         fromHeight = height_median_parent,
         toHeight = height_median,
         fromRegion = collection_regionname_parent,
         toRegion = collection_regionname,
         fromHost = host_simplifiedhost_parent,
         toHost = host_simplifiedhost
         ) %>%
  
  relocate(fromNode,
           toNode,
           fromHeight,
           toHeight,
           fromRegion,
           toRegion,
           fromHost,
           toHost) %>%
  
  mutate(persistencetime_host = case_when(toHost == fromHost ~ as.numeric(fromHeight) - as.numeric(toHeight),
                                          .default = 0),
         persistencetime_region = case_when( toRegion == fromRegion ~ as.numeric(fromHeight) - as.numeric(toHeight),
                                            .default = 0)
         )


