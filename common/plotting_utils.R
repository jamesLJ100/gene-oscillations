library(ggplot2)

#' Shared minimal ggplot theme for score-comparison violin/boxplots
#'
#' @param base_size Base font size (default 16)
#' @param angle_x_labels Logical, rotate x-axis labels 45 degrees (useful when
#'   cluster/condition names are long, e.g. nfkb's "wt_tnf"/"ss_pic")
#' @return A ggplot2 theme object
theme_common <- function(base_size = 16, angle_x_labels = FALSE) {
  t <- theme_bw(base_size = base_size) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.ticks.length = unit(-0.25, "cm"),
      legend.position   = "none",
      strip.background  = element_rect(fill = "transparent", colour = NA)
    )
  if (angle_x_labels) {
    t <- t + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  t
}
