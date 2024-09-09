# Dependencies
library(tidyverse)

# Import Tree Data Extractions
csv_paths <- list.files('./2024Aug18/treedata_extractions',
                        pattern = '.csv',
                        full.names = TRUE)

csv_names <- gsub('.*treedata_extractions/|.csv$',
                  '',
                  csv_paths)

csv_data <- lapply(csv_paths, read_csv) %>%
  setNames(csv_names)

# Expand to separate objects
for(i in 1:length(csv_data)){
  assign(names(csv_data)[i],as.data.frame(csv_data[[i]]))
}


# Persistence time of reassortants

