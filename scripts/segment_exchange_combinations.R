################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-06-10
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)

test <-  function(x) {
  parent_col = gsub('cluster_', 'parent_', cur_column())
  if_else(x == .data[[parent_col]], "0", "1")
}

# Function to compute Jaccard between 2 columns
jaccard <- function(x, y) {
  sum(x & y) / sum(x | y)
}


################################### DATA #######################################
# Read and inspect data
segments <- c('PB2', 'PB1', 'PA', 'HA', 'NP',  'NA','M', 'NS')

################################### MAIN #######################################
# Main analysis or transformation steps

segments_change <- updated %>%
  separate_wider_delim(cluster_profile, delim = '_', names = paste0('cluster_',segments), cols_remove = F) %>%
  separate_wider_delim(parent_profile, delim = '_', names = paste0('parent_',segments), cols_remove = F) %>%
  separate_wider_delim(last_major_profile, delim = '_', names = paste0('last_major_',segments), cols_remove = F) %>%
  
  # changes since last reassortant
  mutate(PB2_switch = if_else(cluster_PB2 == parent_PB2, '0', '1'),
         PB1_switch = if_else(cluster_PB1 == parent_PB1, '0', '1'),
         PA_switch = if_else(cluster_PA == parent_PA, '0', '1'),
         HA_switch = if_else(cluster_HA == parent_HA, '0', '1'),
         NP_switch = if_else(cluster_NP == parent_NP, '0', '1'),
         NA_switch = if_else(cluster_NA == parent_NA, '0', '1'),
         M_switch = if_else(cluster_M == parent_M, '0', '1'),
         NS_switch = if_else(cluster_NS == parent_NS, '0', '1')) %>%
  
  # changes since last major reassortant
  mutate(PB2_switch_mj = if_else(cluster_PB2 == last_major_PB2, '0', '1'),
         PB1_switch_mj = if_else(cluster_PB1 == last_major_PB1, '0', '1'),
         PA_switch_mj = if_else(cluster_PA == last_major_PA, '0', '1'),
         HA_switch_mj = if_else(cluster_HA == last_major_HA, '0', '1'),
         NP_switch_mj = if_else(cluster_NP == last_major_NP, '0', '1'),
         NA_switch_mj = if_else(cluster_NA == last_major_NA, '0', '1'),
         M_switch_mj = if_else(cluster_M == last_major_M, '0', '1'),
         NS_switch_mj = if_else(cluster_NS == last_major_NS, '0', '1')) %>%
  
  select(-ends_with(c('PB2', 'PB1', 'PA', 'HA', 'NP',  'NA','M', 'NS')))


# total number of changes
segments_change %>%
  select(ends_with('switch')) %>%
  pivot_longer(everything(), values_to = 'switch', names_to = 'segment') %>%
  mutate(segment = gsub('_switch', '', segment),
         switch = as.numeric(switch)) %>%
  summarise(switch = sum(switch), .by = segment) %>%
  ggplot(aes(x = segment, y = switch)) +
  geom_bar(stat = 'identity')


# Linkage between reassortments (ie which seg)
tmp <- segments_change %>%
  select(ends_with('switch')) %>%
  mutate(across(everything(), .fns = ~as.numeric(.x))) %>%
  
  mutate(row_id = row_number()) %>%                              # Add row ID
  pivot_longer(-row_id, names_to = "segment", values_to = "selected") %>%
  mutate(segment = gsub('_switch', '', segment)) %>%
  filter(selected == 1) %>%                                     # Keep only selected = 1
  inner_join(., ., by = "row_id") %>%                           # Self-join on row_id
  filter(segment.x < segment.y) %>%                                     # Avoid duplicates and self-pairs
  count(segment.x, segment.y, name = "n") %>%                           # Count co-occurrences
  pivot_wider(names_from = segment.y, values_from = n,              # Make it wide
              values_fill = 0)

totals <- segments_change %>%
  select(ends_with('switch')) %>%
  mutate(across(everything(), .fns = ~as.numeric(.x))) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(everything(), names_to = "segment",  values_to = "selected_count") %>%
  mutate(segment = gsub('_switch', '', segment)) 


## Calculate Jaccard
df <- segments_change %>%
select(ends_with('switch')) %>%
  mutate(across(everything(), .fns = ~as.numeric(.x))) 

segments <- colnames(df) 

observed_jaccards <- expand.grid(seg1 = segments, seg2 = segments) %>%
  rowwise() %>%
  mutate(jaccard = jaccard(df[[seg1]], df[[seg2]])) %>%
  ungroup() %>%
  mutate(jaccard = if_else(jaccard == 1, NaN, jaccard),
         across(starts_with('seg'), ~ gsub('_switch', '', .x)))

# Permutation test
permute_jaccards_long <- function(df, segments, n_perm = 1000) {
  combs <- expand.grid(seg1 = segments, seg2 = segments) %>% 
    filter(seg1 != seg2)
  
  replicate(n_perm, {
    df_perm <- df %>% mutate(across(all_of(segments), ~ sample(.x)))
    combs %>%
      rowwise() %>%
      mutate(jaccard = jaccard(df_perm[[seg1]], df_perm[[seg2]])) %>%
      pull(jaccard)
  }) %>% t()  # matrix: rows = permutations, cols = segment pairs
}

perm_matrix <- permute_jaccards_long(df, segments, n_perm = 1000)

valid_pairs <- observed_jaccards %>% filter(seg1 != seg2)

p_vals <- map2_dbl(
  seq_len(nrow(valid_pairs)),
  valid_pairs$jaccard,
  ~ mean(perm_matrix[, .x] >= .y, na.rm = TRUE)
)

observed_jaccards %<>%
  left_join(
    bind_cols(valid_pairs[, c("seg1", "seg2")], 
              p_value = p.adjust(p_vals, method = "BH")),
    by = c("seg1", "seg2")
  )


## permutation test for polymerase complex
observed_rdrp <- sum(df$PB2_switch & df$PB1_switch & df$PA_switch)
n <- nrow(df)

# Estimate marginal probabilities
p_pb2 <- mean(df$PB2_switch)
p_pb1 <- mean(df$PB1_switch)
p_pa  <- mean(df$PA_switch)
expected_rdrp <- p_pb2 * p_pb1 * p_pa

# Binomial test
binom.test(observed_rdrp, n, expected_rdrp, alternative = "greater")


## permutation test for polymerase complex
observed_rdrp <- sum(df$PB2_switch & df$PB1_switch & df$PA_switch & df$NP_switch)
n <- nrow(df)

# Estimate marginal probabilities
p_pb2 <- mean(df$PB2_switch)
p_pb1 <- mean(df$PB1_switch)
p_pa  <- mean(df$PA_switch)
p_np <- mean(df$NP_switch)
expected_rdrp <- p_pb2 * p_pb1 * p_pa * p_np

# Binomial test
binom.test(observed_rdrp, n, expected_rdrp, alternative = "greater")

################################### OUTPUT #####################################
# Save output files, plots, or results
plt_a <- segments_change %>%
  select(ends_with('switch')) %>%
  pivot_longer(everything(), values_to = 'switch', names_to = 'segment') %>%
  mutate(segment = gsub('_switch', '', segment),
         switch = as.numeric(switch)) %>%
  summarise(switch = sum(switch), .by = segment) %>%
  ggplot(aes(x = segment, y = switch, fill = switch)) +
  scale_x_discrete('Segment' ,expand = c(0,0)) + 
  scale_y_continuous('Reassortants (n)', expand = c(0,0)) + 
  scale_fill_distiller(palette = 'YlOrRd', na.value = NA, direction = 1) +
  geom_bar(stat = 'identity', colour = 'black') + 
  theme_classic() + 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))


plt_b <- observed_jaccards %>%
  ggplot() + 
  geom_point(aes(x = seg1, y = seg2, fill = jaccard, size = jaccard), shape = 21) +
  scale_size( range = c(0.5, 10)) + 
  #geom_tile(aes(x = seg1, y = seg2, fill = jaccard)) +
  scale_fill_distiller(palette = 'YlOrRd', na.value = NA, direction = 1) +
  scale_x_discrete('Segment' ) + 
  scale_y_discrete('Segment' ) + 
  theme_classic()+ 
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.position = 'none')

cowplot::plot_grid(plt_a + theme(legend.position = 'none'), 
                   plt_b ,
                   labels = 'AUTO',
                   label_size = 9, 
                   align = 'hv', 
                   axis = 'rltb',
                   ncol = 2)

ggsave('~/Downloads/flu_plots/figure_segment_exchange.jpeg', height = 10, width = 20, units = 'cm', dpi = 360)
#################################### END #######################################
################################################################################