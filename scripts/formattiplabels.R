library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)

pb2_tree <- read.newick('./data/phylo_ml/RAxML_bipartitions.H5_large_pb2')
ggtree(pb2_tree)

metadata <- read_csv('./data/metadata/h5_metadata.csv')
will_metadata <- read_csv('./data/metadata/Re_H5_clade2344b_Asia_blast_PB2_meta.will.2.csv')

t <- left_join(metadata, will_metadata )
tip_labels <- TipLabels(pb2_tree) 


joined_data <- tip_labels %>%
  sapply(., str_split_4, pattern = '\\|') %>%
  sapply(., `[`, 4) %>%
  bind_cols() %>%
  setNames('isolate_id') %>%
  left_join(metadata)

prob_seqs <- joined_data %>% filter(is.na(isolate_name))

tree <- pb2_tree
n<-length(tree$tip.label)
ee<-setNames(tree$edge.length[sapply(4:n,function(x,y)   which(y==x),y=tree$edge[,2])],tree$tip.label)


data_from_tipnames <- tip_labels %>%
  sapply(., str_split_4, pattern = '\\|') %>%
  sapply(., function(x) ifelse(!grepl('\\.', tail(x, n = 4)), return(x), NA)) %>%
  .[!is.na(.)] %>%
  
  do.call(rbind.data.frame,. )


#Isolate id = starts with EPI_ISL
#subtype = HxNx
#isolate name = 
#segment
#date

MyFunc4 <- function(x){
  
  isolate.id = x[grep('EPI_ISL_*' ,x)]
  subtype = x[grep('H[[:digit:]]*N[[:digit:]]*' ,x)]
  segment = x[grep('HA|PB4|PB2|PA|NP|NA|NS' ,x)]
  date = x[grep('[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]' ,x)]
  isolate.name = x[grep('/[^/]*/[^/]*/' ,x)]
  
  if(identical(isolate.id, character(0))){
    isolate.id <- NA
  }
  
  if(identical(subtype, character(0))){
    subtype <- NA
    print(subtype)
  }
  if(identical(segment, character(0))){
    segment <- NA
  } 
  if(identical(date, character(0))){
    date <- NA
  } 
  if(identical(isolate.name, character(0))){
    isolate.name <- NA
  }
  
  out <- cbind.data.frame('isolate.id' = isolate.id,
                    'subtype' = subtype,
                    'segment' = segment,
                    'date' = date,
                    'isolate.name' = isolate.name)
  
  
  return(out)
}

df <- lapply(data_from_tipnames, MyFunc4) %>% bind_rows()

data.frame('isolate.id' = data_from_tipnames[[4]][grep('EPI_ISL_*' ,data_from_tipnames[[4]])],
           'subtype' = data_from_tipnames[[4]][grep('H[[:digit:]]N[[:digit:]]|H[[:digit:]][[:digit:]]N[[:digit:]][[:digit:]]' ,data_from_tipnames[[4]])],
           'segment' = data_from_tipnames[[4]][grep('HA|PB4|PB2|PA|NP|NA|NS' ,data_from_tipnames[[4]])],
           'date' = data_from_tipnames[[4]][grep('[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]' ,data_from_tipnames[[4]])],
           'isolate.name' = data_from_tipnames[[4]][grep('/[^/]*/[^/]*/' ,data_from_tipnames[[4]])])

