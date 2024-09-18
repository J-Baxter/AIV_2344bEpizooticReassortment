####################################################################################################
# Infer global alignments
library(tidyverse)
library(treeio)
library(ape)


####################################################################################################
# Import reference tree
will_tree <- read.beast('./global_subsample/h52344b_ha_s1.mcc.trees')
metadata <- read_csv('./2024-09-09_meta.csv')
included_isolates <- str_extract(will_tree@phylo$tip.label,  "EPI_ISL_(china_){0,1}\\d+[^.|]*") %>%
  .[!. =='EPI_ISL_13370925'] #dropped as no NP

# Import all region alignments (as list)
region_alignmentfiles <- list.files('./global_subsample/region_alignments', 
                                    pattern = '.fasta',
                                    full.names = T)

region_alignments <- lapply(region_alignmentfiles, read.dna,
                            format = 'fasta',
                            as.matrix = FALSE) 

names(region_alignments) <- gsub('.*\\/h5_|.fasta$', '', region_alignmentfiles)


# concatenate all alignments (within segment)
segments <- gsub('_.*', '', names(region_alignments)) %>%
  gsub('n[:0-9:]', 'n[:0-9:]', .) %>%
  unique()


####################################################################################################
# initialise empty list
selected_genomes <- rep_along(segments, list())

# select isolates

for (i in 1:length(selected_genomes)){
  segment <- segments[i]
  
  relevant_alns <- region_alignments[grep(segment, names(region_alignments))] 
  
  temp <-  do.call(c, relevant_alns) 
  names(temp) <- gsub('.*\\.', '', names(temp))
  
  selected_genomes[[i]] <- temp[str_extract(names(temp), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% included_isolates]
  
}


selected_genomes_unique <- selected_genomes %>%
  lapply(., function(x) x[unique(names(x))]) 


# confirm all isolates in all fasta files
length(included_isolates) == sapply(selected_genomes_unique, length, simplify = F) 


####################################################################################################
# Stratify dataset by a) reassortant and b) sampling through time
# aim to restrict n_sequences to ~ 500

# Sample first and last sequence from reassortant
metadata_subsampled <- metadata %>%
  
  # include only those data present in Will's tree
  filter(isolate_id %in% included_isolates) %>%
  group_by(cluster_profile,
           collection_datemonth) %>%
  
  # Sample first and last sequence from reassortant
  slice(which(collection_date == max(collection_date, na.rm = T) | collection_date == min(collection_date, na.rm = T)), 
        .preserve = T) %>%
  group_by(collection_date, .add = TRUE) %>%
  slice_sample(n = 1) %>%
  ungroup() 


####################################################################################################
# Update alignments
alignments_subsampled <- lapply(selected_genomes_unique,
                                function(x) x[str_extract(names(x), "EPI_ISL_(china_){0,1}\\d+[^.|]*") %in% metadata_subsampled$isolate_id]) 



####################################################################################################
# export all fasta for alignment and manual editing
alignmentfiles_subsampled <- paste('./2024Sept16/global_strictdownsample',
                                   paste(segments, 'global.fasta', sep = '_'),
                                   sep = '/' )

mapply(ape::write.dna, 
       alignments_subsampled, 
       alignmentfiles_subsampled ,
       format = 'fasta')

