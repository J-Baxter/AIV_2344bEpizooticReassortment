require(extrafont)


#https://www.fontsquirrel.com/fonts/latin-modern-sans
#font_import(pattern = "lmsans10*") 
#loadfonts()

#my_theme <- theme_classic(base_family = "LM Sans 10")+
 # theme(
    #text = element_text(size=10),
    #axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
    #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
    #axis.text = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 7),
    #strip.text  = element_text(size = 8),
    #legend.text = element_text(size = 7),
    #legend.title = element_text(size = 8),
    #legend.position = 'none', 
    #panel.spacing = unit(2, "lines"), 
    #strip.background = element_blank()
 # )

global_theme <- theme_classic()+
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

region_colours <- c('europe' = '#1b9e77',
                    'asia' ='#d95f02',
                    'africa' ='#7570b3',
                    'australasia' = '#e7298a',
                    'central & northern america' ='#66a61e',
                    'south america' ='#e6ab02')

host_colours <- c(
  'anseriformes-domestic' = '#a6cee3',
  'anseriformes-wild' = '#1f78b4',
  'galliformes-domestic' = '#b2df8a',
  'galliformes' = '#33a02c',
  'mammal' = '#fb9a99',
  'human' = '#e31a1c',
  'charadriiformes-wild' = '#fdbf6f',
  'other-bird' = '#ff7f00',
  'unknown' = '#cab2d6',
  'environment' = '#6a3d9a')


class_colours <- c('major' = '#A50104', 
                   'moderate' = '#D17A22', 
                   'minor' = '#B4C292')