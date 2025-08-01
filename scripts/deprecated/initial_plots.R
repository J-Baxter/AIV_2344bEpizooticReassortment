library(tidyverse)


dominant_persistencetime <- domainantTre_mrca_stats_disp_combined %>%
  select(c(data_name,
           youngestTip.time,
           TMRCA,
           Rcode)) %>%
  left_join(summary_reassortant_metadata_20240902,
            by = join_by(Rcode == cluster_profile)) %>%
  select(c(data_name,
           TMRCA,
           Rcode,
           Length_Between_First_Last_Sample,
           youngestTip.time)) %>%
  separate_wider_delim(data_name, '_', names = c('segment', 'cluster_profile'))

ggplot(dominant_persistencetime) +
  geom_linerange(aes(x = cluster_profile, ymin = youngestTip.time-Length_Between_First_Last_Sample, ymax = youngestTip.time), linewidth = 0.5, colour = 'red') +
  geom_linerange(aes(x = cluster_profile, ymin = TMRCA, ymax = youngestTip.time), linewidth = 0.5, alpha = 0.5, colour = 'red') +
  
  coord_flip() + 
  facet_grid(rows = vars(segment)) + 
  theme_bw()


ggplot(dominant_persistencetime) +
  geom_histogram(aes(x = youngestTip.time-Length_Between_First_Last_Sample-TMRCA)) +
  facet_grid(rows = vars(cluster_profile))