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
  if (rownames(est_df)[1] == "lower") {
    case <- rep(1, 2)
    bound_df <- data.frame(x1 = 1,
                           y1 = min(est_df["lower",]),
                           x2 = 1,
                           y2 = max(est_df["upper",]))
    } else {
    case <- 1:nrow(est_df)
    bound_df <- data.frame(x1 = 1:nrow(est_df),
                           y1 = apply(est_df, 1, min),
                           x2 = 1:nrow(est_df),
                           y2 = apply(est_df, 1, max))
    }
  df <- data.frame(est_df, case = case) %>%
    gather(key = "Type", value = "effect", - case)
  ggplot(df) +
  geom_point(aes(x = case, y = effect, col = Type), size = 1.5)  +
  geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.3) +
  scale_color_discrete(name = expression(R^2~":"),
                       labels = paste0(as.numeric(sub('...', '', levels(factor(df$Type))))*100, "%")) +
  labs(y = "estimates") +
  scale_x_continuous(breaks = unique(case), labels = as.character(unique(case)),
                     limits = c(0.5, max(case) + 0.5)) +
  theme_bw(base_size = 12)
}



