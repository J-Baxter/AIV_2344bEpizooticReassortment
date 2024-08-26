library(beastio)
library(tidyverse)


logfiles <- list.files('./2024Jul12/region_beastlog',
                       pattern = 'africa',
                       full.names = T)


imported_logs <- lapply(logfiles, function(x) try(readLog(x)))


probs <- cbind.data.frame('filename' = rep(NA, length(imported_logs)),
                          'complete' = rep(NA, length(imported_logs)),
                          'problems' = rep(NA, length(imported_logs)))

for( i in 1:length(imported_logs)){
  print(i)
  probs[i,1] <- logfiles[[i]]
  probs[i,2] <- nrow(imported_logs[[i]]) > 9000
  
  lowESS <- try(checkESS(imported_logs[[i]]))
  if(length(lowESS)==0){
    probs[i,3] <- NA
  }else{
    lowESS_char <- paste0(lowESS %>% names(), ' (', lowESS, ")")
    probs[i,3] <- paste(lowESS_char, collapse = ', ')
  }
  
}

probs_edited <- probs %>%
  mutate(runtype = 'traits') %>%
  mutate(region = str_split(probs$filename, '_')%>% lapply(., `[`, 2) %>% unlist()) %>%
  relocate(runtype,
          region) %>%
  as_tibble()

#write_csv(probs_edited, 'checklogfiles_apr18.csv')

treefiles <- list.files('./2024Jul12/region_beasttreefile',
                        pattern = 'africa',
                        full.names = T)

# Patterns to match against the file paths
patterns <- str_split(treefiles,  '/') %>% 
  lapply(., tail, n = 1) %>% 
  #lapply(., `[`, c(1,2)) %>%
  lapply(., function(x) gsub("^([^_]*_[^_]*).*", "\\1",x)) %>%
  unlist() %>%
  tolower() %>%
  gsub('na', 'n', .) %>%
  unique()

# Initialize an empty list to store the results
group_list <- list()

# Iterate over each pattern
for (pattern in patterns) {
  # Filter the file paths based on the current pattern
  matching_files <- treefiles[grepl(pattern, treefiles)]
  
  all_files <- c(matching_files, gsub('[^_]*$', 'combined1800.trees', matching_files[1]))
  # Store the filtered file paths in the result list
  group_list[[pattern]] <- all_files
}


logcombiner<- lapply(group_list, function(x) paste('/Applications/bin/logcombiner -burnin 10000000 -trees -resample 100000', x[1], x[2], x[3])) %>% 

  lapply(., function(x) gsub('./2024Jul12/region_beasttreefile/' , '', x)) %>%
  unlist()

write_lines(logcombiner,
            './2024Jul12/region_beasttreefile/logcombiner.sh'
)









