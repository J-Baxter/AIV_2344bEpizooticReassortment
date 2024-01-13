####################################################################################################
# Basic temporal signal analysis
# NB: inferior analysis to TEMPEST as doesn't account for variable precision in dates and only uses
#     RSM (not heuristic RSM) to infer best-fitting root

# Requires a phylogeny (branch lengths inferred from genetic distance only) with tipdates


####################################################################################################
# Package dependencies
library(ape)
library(tidyverse)
library(treeio)
library(TreeTools)

# Functions

# Filter out sequences with no date
CleanPhylo <- function(phylo){
  out <- drop.tip(phylo,
                  TipLabels(phylo)[grep('NA$', TipLabels(phylo))])
  
  return(out)
}


CleanandRoot <- function(phylo, metadata){
  require(ape)
  require(TreeTools)
  require(adephylo)
  
  stopifnot(Ntip.phylo(phylo) == nrow(metadata))
  
  phylo_dropped <- CleanPhylo(phylo = phylo)
  
  r2t <- ape::rtt(phylo_dropped,
                  objective = 'rms', 
                  tip.dates = metadata$collection.datedecimal[!is.na(metadata$collection.datedecimal)])
  
  rttdist <- adephylo::distRoot(r2t)
  
  df <-  cbind.data.frame('dist' = rttdist, 'tipnames' = names(rttdist)) %>% 
    left_join(metadata)
  
  return(df)
  
}


RTTReg <- function(x){
  my_model <- lm(dist ~ collection.datedecimal, data = x)
  rate <-  coef(my_model)[2]
  x_intercept <- -coef(my_model)[1] / coef(my_model)[2] 
  
  out <- cbind.data.frame(rate, x_intercept) %>%
    as_tibble()
  return(out)
}

test_treedata <- tidytree::as_tibble(trees[[1]]) %>%
  left_join(., HA_strict, by = join_by(label == tipnames)) %>%
  as.treedata()
####################################################################################################
# Import Phylogenies
# Import tree data
treefiles <- list.files(path = './data/alignments/subsampled_alignments/2024Jan10_strict/iqtrees',
                        recursive = FALSE,
                        include.dirs = FALSE, 
                        full.names = TRUE) %>%
  .[grep('treefile$', .)]

trees <- lapply(treefiles, read.newick) 

# Create treedata object (joins dataframe and phylo into one object)


####################################################################################################
# Infer best-fitting root and get rtt distances (returns list of dataframes)

rtt_dist <- mapply(CleanandRoot,
                   phylo = renamed_phylos,
                   metadata = formatted_metadata,
                   SIMPLIFY = FALSE) %>%
  setNames(segnames) 


t <- rtt_dist %>% bind_rows(., .id = 'segment') %>%
  left_join(ratesandintercept)

# Calculate gradient (evolutionary rate) and intercept (tmrca)
ratesandintercept <- lapply(rtt_dist, rtt_reg) %>%
  bind_rows(.id = 'segment') %>%
  left_join(ratesandintercept)

# Plot
ggplot(t, aes(x = collection.datedecimal, y = dist)) +
  geom_point(aes(colour = !is.na(cluster.profile))) + 
  scale_y_continuous('Genetic Distance', 
                     expand = expansion(mult = c(0,0.05))) + 
  scale_x_continuous('Time') +
  scale_colour_brewer(
    'Reassortment',
    palette = 'Paired', 
    labels = c('FALSE' = 'Non Reassortment',
               'TRUE' = 'Reassortment'))+
  geom_point(aes(x = x_intercept, y = 0)) + 
  geom_smooth(method = 'lm', 
              fullrange = T,
              colour = 'black') +
  coord_cartesian(ylim = c(0, NA)) +
  geom_text(aes(x = x_intercept,
                y = Inf,
                label = paste('x intercept =', round(x_intercept, 2))),
            size = 3,
            vjust= +3, 
            hjust = -0.2, 
            data = ratesandintercept) +
  geom_text(aes(x = x_intercept,
                y = Inf,
                label = paste('rate =', formatC(rate, format = "e", digits = 2))), 
            size = 3,
            vjust= +5,
            hjust = -0.27,
            data = ratesandintercept) + 
  facet_wrap(.~segment, scales = 'free') +
  my_theme + 
  theme(legend.position = 'bottom')


###################
#t <- rooted_phylos %>%
# lapply(., as_tibble) %>%
# {lapply(seq_along(along.with = .), 
#        function(i)  left_join( x = .[[i]], 
#                                y = formatted_metadata[[i]], 
#                                join_by(label == tipnames)))} %>% 
#lapply(., as.treedata) %>%
# lapply(., function(x) ggtree(x) + geom_tippoint(aes(colour = is.na(cluster.genome))) +
#         scale_colour_brewer(
#           'Reassortment',
#           palette = 'Paired', 
#           direction = -1,
#           labels = c('TRUE' = 'Non Reassortment',
#                      'FALSE' = 'Reassortment')) + theme(legend.position = 'none')) %>%
# cowplot::plot_grid(plotlist = ., align = 'hv', labels = segnames, ncol = 3)


#t_mat <- rooted_phylos %>%
#  lapply(., as_tibble) %>%
#  {lapply(seq_along(along.with = .), 
#          function(i)  left_join( x = .[[i]], 
#                                 y = formatted_metadata[[i]], 
#                                 join_by(label == tipnames)))} %>% 
#  lapply(., as.treedata) %>%
#  lapply(., function(x) ggtree(x) + 
#          theme(legend.position = 'none') %>%
#           gheatmap(., is.na(x@data$cluster.genome),
#                   offset=8, 
#                    width=0.6, 
#                    colnames=FALSE) ) %>%
#  cowplot::plot_grid(plotlist = ., align = 'hv', labels = segnames, ncol = 3)
###################
# Remove tip names identified as problems



