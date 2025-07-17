################################################################################
## Script Name:        Downsampled diversity plot
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2025-07-17
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


################################### DATA #######################################
# Read and inspect data
down_sampled <- read_csv('./profile_diversity/04h5_diversity_downsample/20250605_h5_diversity_downsample.csv')

################################### MAIN #######################################
# Main analysis or transformation steps
down_sampled %>%
  filter(! continent %in%  c('South America', 'Antarctica')) %>%
  mutate(continent = str_to_lower(continent) %>%
           case_when(grepl('north america', .) ~ 'central & northern america',
                     .default = .)) %>%
  mutate(midpoint = (start_point + end_point)/2) %>%
  mutate(midpoint_date = date_decimal(midpoint)) %>%
  #filter(random == 1) %>%
  ggplot(aes(y = diversity, 
             x = midpoint_date,
             group = interaction(continent, as.factor(random)), 
             colour = continent)) +
  #geom_point() + 
  geom_line(stat="smooth",
            method = "gam",
            formula = y ~ s(x, bs = "cs"),
            size = 0.8,
            se = FALSE,
            alpha = 0.1) + 
  #geom_labelsmooth(text_smoothing = 30, 
  # fill = "white",
  #formula = y ~ s(x, bs = 'cs'), 
  #method = "gam",
  #size = 4, 
  #linewidth = 0.1,
  #boxlinewidth = 0.3) +
  #scale_fill_manual('Continent', values = region_colours, labels = str_to_title) + 
  scale_colour_manual('Continent', values = region_colours, labels = str_to_title) +
  scale_y_continuous('Numbers-Equivalent Shanon Entropy', expand = c(0.01, 0)) + 
  scale_x_datetime(limits = as_datetime(c('2020-06-01', '2024-02-01')),
                   breaks = '1 year', 
                   date_labels = "%Y", 'Date (Years)',
                   expand = c(0,0)) + 
  
  theme_classic()+ 
  theme(legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

################################### OUTPUT #####################################
# Save output files, plots, or results
ggsave('~/Downloads/flu_plots/diversity_downsample.jpeg', 
       height = 10,
       width = 12, 
       units = 'cm', 
       dpi = 360)

#################################### END #######################################
################################################################################