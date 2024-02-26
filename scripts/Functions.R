#### CREATION OF MY THEME ####
my_theme <- theme(
  plot.title = element_text(size = rel(2)),
  panel.grid.major.y = element_line(color = 'gray90', linewidth = 0.3),
  panel.grid.minor.y = element_line(color = 'gray90', linewidth = 0.1),
  panel.grid.major.x = element_line(color = 'gray90', linewidth = 0.3),
  panel.grid.minor.x = element_line(color = 'gray90', linewidth = 0.1),
  panel.border = element_rect(color = 'gray50', fill = NA, linewidth = 1),
  plot.background = element_rect(fill = NULL),
  panel.background = element_rect(fill = 'gray99'),
  axis.line = element_line(color = 'gray50'),
  axis.text = element_text(color = 'gray40', face = 'bold'),
  axis.text.x = element_text(size = rel(1.2)),
  axis.text.y = element_text(size = rel(1.2)),
  axis.title = element_text(size = rel(1), face = 'bold'),
  axis.ticks = element_line(color = 'gray50'),
  legend.position = "right",
  legend.margin = margin(5, 5, 5, 5),
  legend.title = element_text(face = 'bold'),
  legend.background = element_blank(),
  legend.key = element_rect(fill = 'gray98'),
  legend.box.background = element_rect(fill = 'gray98',color = "gray50")
)