# Basic filtering of maximum likelihood phylogenies
# A function to create an initial subsample of nodes based on the following broad criteria:
# 1. Genetic distance - select a single example of trait/location/date from a clade in which 
#     all descending branches are less than 0
# 2. Remove outliers (>2sd from the root-to-tip regression)


# Filter phylogenetic tree according to branch length
FilterBL <- function(phylo,  bl_threshold = 0){
  require(ape)
  
  return(tiplist)
}

# filter phylogenetic tree accoding to root-to-tip regression
FilterR2T <- function(phylo,  sd_threshold = 0){
  require(ape)
  
  return(tiplist)
}

# root to tip regression plots
pb2_tree_dropped <- drop.tip(pb2_tree, TipLabels(pb2_tree)[grep('NA$', TipLabels(pb2_tree))])
rttdist <- adephylo::distRoot(r2t)

test_df <-  cbind.data.frame('dist' = rttdist, 'tipnames' = names(rttdist)) %>% 
  left_join(df )

my_model <- lm(dist ~ collection.datedecimal, data =test_df)
rate = coef(my_model)[2]
x_intercept <- -coef(my_model)[1] / coef(my_model)[2] #y = mx+c to y = m(x-d) i.e d = (c\times-1)/m

ggplot(test_df, aes(x = collection.datedecimal, y = dist)) +
  geom_point(aes(colour = host.order)) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  geom_point(aes(x = x_intercept, y = 0)) + 
  geom_smooth(method = 'lm', fullrange = T) +
  coord_cartesian(xlim = c(floor(x_intercept/10) * 10, 2030), ylim = c(0, NA)) +
  annotate('text', x = x_intercept, y = 0.19, label = paste('x intercept =', round(x_intercept, 2))) +
  annotate('text', x = x_intercept, y = 0.18, label = paste('rate =', round(rate, 6))) + 
  theme_classic()
SubSamplePhylo <- function(phylo, 
                           bl_threshold = 0, 
                           sd_threshold = 2){
  
  require(ape)
  # If n>2
  # IF (all lengths to all tips from node x are less than y ) = rmultinom
  if(class(phylo) == 'treedata'){
    phylo = phylo@phylo
    
  }else if(class(phylo) != 'phylo'){
    warning('Object is not a phylogenetic tree of class phylo.')
    stop()
    
  }
  
  
  return(tiplist)
}
