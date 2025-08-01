FindIdenticalSeqs<- function(alignment, snp_threshold = 1){
  
  # calculate pairwise distances
  distances <- ape::dist.dna(alignment, model = "N")
  
  # Identify identical sequences
  identical_pairs <- which(as.matrix(distances) < snp_threshold, arr.ind = TRUE) %>%
    dplyr::as_tibble(rownames = 'tipnames') %>% 
    filter(row !=col) 
  
  # infer adjcency network 
  sequence_net <- igraph::graph_from_data_frame(identical_pairs %>% 
                                                  select(!tipnames), 
                                                directed = F)
  
  # returns a dataframe of tipnames and group number (int)
  groupings <- merge(identical_pairs,
                 stack(igraph::clusters(sequence_net)$membership),
                 by.x = "row", 
                 by.y = "ind", all.x = TRUE) %>%
    dplyr::select(c(tipnames, values)) %>%
    dplyr::distinct() %>%
    dplyr::rename(group = values)
  
  return(groupings)
}

#identical_seqs <- lapply(subsetted_alignments, FindIdenticalSeqs)

