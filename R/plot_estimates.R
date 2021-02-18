#' Visualize Estimates of Treatment Effects
#'
#' @param est an return object from \code{\link{gcalibrate}} or \code{\link{bcalibrate}}, or
#' \code{data.frame} containing estimates of treatment effects with estimates' type in columns and
#' contrasts of interest in rows.
#'
#' @export
#'
#' @examples
#' Please see \code{\link{bcalibrate}} or \code{\link{gcalibrate}}

plot_estimates <- function(est) {
  if(is.list(est) & (!is.data.frame(est))) {
    est_df <- est$est_df
    R2 <- round(c(0, sort(est$R2)), 2)
  } else if (is.data.frame(est)){
    est_df <- est
    R2 <- as.numeric(sub('...', '', colnames(est))) %>%
      sort() %>% round(2)
  }
  case <- 1:nrow(est_df)
  bound_df <- data.frame(x1 = 1:nrow(est_df),
                         y1 = apply(est_df, 1, min),
                         x2 = 1:nrow(est_df),
                         y2 = apply(est_df, 1, max))
  if (length(R2) < ncol(est_df)) {
    df_wc <- unlist(cbind(est_df[,1], est_df)) %>%
      matrix(nrow = 2*nrow(est_df))
    colnames(df_wc) <- c(0, est$R2)
    df <- data.frame(df_wc, case = rep(case, 2)) %>%
      gather(key = "Type", value = "effect", - case)
  } else {
    df <- data.frame(est_df, case = case) %>%
      gather(key = "Type", value = "effect", - case)
  }



  # labels <- as.character(c(0, results$R2))
  # label_init <- levels(factor(df$Type))
  # if (length(gregexpr('_', label_init[2])[[1]]) > 1) {
  #   labels <- rep(0, length(label_init))
  #   for (i in 2:length(label_init)) {
  #     # extract the initial position
  #     initial.position = gregexpr('_', label_init[i])[[1]][1] + 1
  #     # extract the final position
  #     final.position = gregexpr('_', label_init[i])[[1]][2] - 1
  #     # extract the substring between the initial and final positions, inclusively
  #     labels[i] = substr(label_init[i], initial.position, final.position)
  #   }
  # } else {
  #   labels <- sub('...', '', label_init)
  # }

  ggplot(df) +
    ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.5, size = 1.5)  +
    geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.3) +
    labs(y = "estimates") +
    scale_x_continuous(breaks = unique(case), labels = as.character(unique(case)),
                       limits = c(0.5, max(case) + 0.5)) +
    scale_colour_brewer(name = expression(R^2~":"),
                        labels = paste0(R2*100, "%"),
                        palette = "Spectral", direction = -1) +
    theme_bw(base_size = 12)
}



