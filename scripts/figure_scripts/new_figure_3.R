class_colours <- c('major'=        '#A92523'
                   ,'moderate'= '#FA9F42',
                   'minor' = '#2B4162')

g_maj <- make_star(n = 6, mode = "undirected", center = 1)  
theta   <- seq(0, 2*pi, length.out = 6)[-1]               
layout_maj  <- data.frame(
  x = c(0,  cos(theta)),                                 
  y = c(0,  sin(theta))
  )


g_mod <- make_star(n = 3, mode = "undirected", center = 1)  
layout_mod  <- data.frame(
  x = c(0,  0.5, -0.5),                                 
  y = c(0,  0.8660254, 0.8660254)
)


g_min <- make_star(n = 2, mode = "undirected", center = 1)  
layout_min  <- data.frame(
  x = c(-0.5, 0.5),                                 
  y = c(0,  0)
)


panel1 <- ggraph(g_maj, layout = "manual", x = layout_maj[,1], y = layout_maj[,2]) +
  geom_edge_link(colour = "black") +
  geom_node_point(size = 4, colour = "#A92523") +
  theme_void() +
  theme(legend.position = 'none',
       # panel.background = element_rect(fill = alpha('#A92523', 0.5) , colour = NA)
          panel.background = element_rect(colour = '#A92523' , linewidth = 4 )) + 
  annotate("text", x=0, y=4,   label = "Major",    size=4, fontface=2) +
  annotate("text", x=0, y=-3.5, label = "Persist ~37 months", size=3, fontface=1) +
  annotate("text", x=0, y=-4.5, label = "~ 16 `offspring' reassortants",      size=3, fontface=1) + 
  annotate("text", x=0, y=-5.5, label = "Frequent host switches, including mammals", size=3, fontface=1) + 
  scale_y_continuous(limits=c(-7,5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-1.1,1.1), expand = c(0,0)) 

panel2 <- ggraph(g_mod, layout = "manual", x = layout_mod[,1], y = layout_mod[,2]) +
  geom_edge_link(colour = "black") +
  geom_node_point(size = 4, colour = "#FA9F42") +
  theme_void() +
  theme(legend.position = 'none',
        #panel.background = element_rect(fill = alpha('#FA9F42', 0.5), colour =NA)
        panel.background = element_rect(colour = '#FA9F42' , linewidth = 4 )) + 
  annotate("text", x=0, y=4,   label = "Moderate",    size=4, fontface=2) +
  annotate("text", x=0, y=-3.5, label = "Persist ~15 months", size=3, fontface=1) +
  annotate("text", x=0, y=-4.5, label = "~ 2 `offspring' reassortants",      size=3, fontface=1) + 
  annotate("text", x=0, y=-5.5, label = "Periodic host switches, but rarely mammals", size=3, fontface=1) + 
  scale_y_continuous(limits=c(-7,5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-1.1,1.1), expand = c(0,0)) 

panel3 <- ggraph(g_min, layout = "manual", x = layout_min[,1], y = layout_min[,2]) +
  geom_edge_link(colour = "black") +
  geom_node_point(size = 4, colour = "#2B4162") +
  theme_void() +
  theme(legend.position = 'none',
        #panel.background = element_rect(fill = alpha('#2B4162', 0.5) , colour = NA)
        panel.background = element_rect(colour = '#2B4162' , linewidth = 4 )) + 
  annotate("text", x=0, y=4,   label = "Minor",    size=4, fontface=2) +
  annotate("text", x=0, y=-3.5, label = "Persist ~3 months", size=3, fontface=1) +
  annotate("text", x=0, y=-4.5, label = "1 `offspring' reassortant",      size=3, fontface=1) + 
  annotate("text", x=0, y=-5.5, label = "Single host, no mammal infections", size=3, fontface=1) + 
  
  scale_y_continuous(limits=c(-7,5), expand = c(0,0)) +
  scale_x_continuous(limits=c(-1.1,1.1), expand = c(0,0)) 

plt_5a <- grid.arrange(panel1, panel2, panel3, ncol=3)


plt_5b <-avg_predictions(ordinal_model, by = 'cluster_region') %>%
  get_draws(shape = "rvar") |>
  #filter(group == 'major') %>%
  ggplot(aes(x = cluster_region, ydist = rvar, colour = group)) + 
  stat_interval(aes(color_ramp = after_stat(level)), 
                position = position_dodgejust(width = 0.5, justification = 0),
                point_interval = "median_hdci")+
  
  scale_colour_manual(values = class_colours) + 
  
  scale_y_continuous('Pr(Class = k)' ,
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     expand = c(0,0)) +
  
  scale_x_discrete('Continent', expand = c(0,0), labels =  function(x) str_to_title(x) ) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.text = element_text(face = 'bold', size = 10),
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))



### F (Segment)
plt_5c <-avg_predictions(ordinal_model,  variables = list('segments_changed' = 1:7)) %>%
  as_tibble() %>%
  ggplot(aes(ymin = conf.low, ymax = conf.high, y = estimate, x = segments_changed, fill = group,
             colour = group)) + 
  
  geom_ribbon(aes(colour = NULL), alpha = 0.2) + 
  geom_line(linewidth = 1) + 
  
  scale_fill_manual(values = class_colours) + 
  scale_colour_manual(values = class_colours) + 
  
  scale_y_continuous('Pr(Class = k)' ,
                     breaks = seq(0,1,by=0.25),
                     labels = seq(0,1,by=0.25),
                     expand = c(0,0)) +
  
  scale_x_continuous('Segments changed from previous (N)', expand = c(0,0)) + 
  global_theme + 
  theme(strip.placement  = 'inside',
        strip.background = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8))

logo_file <- system.file('~/Downloads/Picture 2.png')

ggdraw() +
  draw_plot(plt_5a, x = 0, y = 0.667, width = 1, height = 0.334) +
  draw_image(
    '~/Downloads/Picture 2.png',x = 0, y = 0.334, width = 1,height = 0.334
  ) + 
  draw_plot(plt_5b, x = 0, y = 0, width = 0.5, height = 0.334) +
  draw_plot(plt_5c, x = 0.5, y = 0, width = 0.5, height = 0.334) +
  draw_plot_label(label = c("A", "B", "C", 'D'), size = 10,
                  x = c(0, 0, 0, 0.5), y = c(1, 0.667, 0.337, 0.337))

plt5_left <- align_plots(plt_5a, plt_5b, align = 'v', axis = 'l')
plt5_bttm <- plot_grid(plt5_left[[2]], plt_5c, nrow = 1, align = 'h', axis = 'tb')
plt5 <- plot_grid(plt5_left[[1]], plt5_bttm, nrow = 2, align = 'v', axis = 'lr')

ggsave('~/Downloads/flu_plots/figure3_new.pdf', height = 257, width = 206, units = 'mm', dpi = 360)
