# Basic filtering of maximum likelihood phylogenies
# A function to create an initial subsample of nodes based on the following broad criteria:
# 1. Genetic distance - select a single example of trait/location/date from a clade in which 
#     all descending branches are less than 0
# 2. Remove outliers (>2sd from the root-to-tip regression)


# Filter phylogenetic tree according to branch length
FilterBL <- function(phylo,  bl_threshold = 0){
  
  return(tiplist)
}

# filter phylogenetic tree accoding to root-to-tip regression
FilterR2T <- function(phylo,  sd_threshold = 0){
  
  return(tiplist)
}


SubSamplePhylo <- function(phylo, 
                           bl_threshold = 0, 
                           sd_threshold = 2){
  
  require(ape)
  
  if(class(phylo) == 'treedata'){
    phylo = phylo@phylo
    
  }else if(class(phylo) != 'phylo'){
    warning('Object is not a phylogenetic tree of class phylo.')
    stop()
    
  }
  
  
  return(tiplist)
}
