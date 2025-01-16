####################################################################################################
####################################################################################################
## Script name: Plot Standards
##
## Purpose of script: Script specifies colour schemes and graphical settings for this work
##
## Date created: 2025-01-16
##
##
########################################## SYSTEM OPTIONS ##########################################
options(scipen = 6, digits = 7) 
memory.limit(30000000) 

  
########################################## DEPENDENCIES ############################################
# Packages
library(tidyverse)
library(magrittr)


# User functions


############################################## DATA ################################################



############################################## MAIN ################################################
global_theme <- theme_classic(base_family = "LM Sans 10")+
  theme(
    #text = element_text(size=10),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 10),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 10),
    axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    legend.text = element_text(size = 8),
    legend.position = 'none', 
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  )


# Host Palette
host_colours <- c(
  'anseriformes-domestic' = '#a6cee3',
  'anseriformes-wild' = '#1f78b4',
  'galliformes-domestic' = '#b2df8a',
  'galliformes-wild' = '#33a02c',
  'mammal' = '#fb9a99',
  'human' = '#e31a1c',
  'charadriiformes-wild' = '#fdbf6f',
  'other-bird' = '#ff7f00',
  'unknown' = '#cab2d6',
  'environment' = '#6a3d9a')


#E8E1E9FF #C0A5AAFF #4D3944FF #7083A4FF #B3A2B4FF #C9CCEAFF #3B3960FF  

#host_colours <- c(
# 'anseriformes-domestic' = '#E8E1E9FF',
# 'anseriformes-wild' = '#C0A5AAFF',
# 'galliformes-domestic' = '#4D3944FF',
# 'galliformes-wild' = '#7083A4FF',
# 'mammal' = '#B3A2B4FF',
# 'human' = '#C9CCEAFF',
# 'charadriiformes-wild' = '#3B3960FF',
# 'other-bird' = '#1E2142FF',
# 'unknown' = '#586085FF',
# 'environment' = '#F6E0D2FF')


# Region Palette
region_colours <- c('europe' = '#1b9e77',
                    'asia' ='#d95f02',
                    'africa' ='#7570b3',
                    'australasia' = '#e7298a',
                    'central & northern america' ='#66a61e',
                    'south america' ='#e6ab02')


############################################## WRITE ###############################################




############################################## END #################################################
####################################################################################################
####################################################################################################