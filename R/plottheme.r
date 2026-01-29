library(ggplot2)
theme_Publication <- function(base_size = 14, base_family = "Helvetica") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(
        face = "bold",
        size = rel(1.2),
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.text = element_text(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "grey85"),
      panel.grid.minor = element_line(colour = "grey92"),
      panel.border = element_blank(),
      panel.background = element_blank(),
      strip.background = element_rect(colour = "grey85", fill = "grey85"),
      strip.text = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", colour = "white"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = rel(0.8)),
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      complete = TRUE
    )
}