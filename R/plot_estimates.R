#' Visualize Estimates of Treatment Effects
#'
#' @param est_df \code{data.frame} containing estimates of treatment effects with estimates' type in columns and
#' contrasts of interest in rows.
#'
#' @export
#'
#' @examples
#' Please see \code{\link{bcalibrate}} or \code{\link{gcalibrate}}

plot_estimates <- function(est_df) {
  bound_df <- data.frame(x1 = 1:nrow(est_df),
                        y1 = apply(est_df, 1, min),
                        x2 = 1:nrow(est_df),
                        y2 = apply(est_df, 1, max))
  data.frame(est_df, case = 1:nrow(est_df)) %>%
    gather(key = "Type", value = "effect", - case) %>%
    ggplot() +
    geom_point(aes(x = case, y = effect, col = Type), size = 1.5)  +
    geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.3) +
    scale_color_discrete(name = expression(R^2~":"),
                         labels = colnames(est_df)) +
    labs(y = "estimates") +
    scale_x_continuous(breaks = 1:nrow(est_df), labels = as.character(1:nrow(est_df)),
                       limits = c(0.5, nrow(est_df) + 0.5)) +
    theme_bw(base_size = 12)
}



