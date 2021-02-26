#' Visualize Estimates of Treatment Effects
#'
#' @param est an return object from \code{\link{gcalibrate}} or \code{\link{bcalibrate}}, or
#' \code{data.frame} containing estimates of treatment effects with estimates' type in columns and
#' contrasts of interest in rows.
#'
#' @param show_rv logical. Whether robustness values should be printed in the plot or not?
#' Available only for the "worstcase" calibration.
#'
#' @export
#'
#' @examples
#' Please see \code{\link{bcalibrate}} or \code{\link{gcalibrate}}

plot_estimates <- function(est, show_rv = TRUE,
                           order = "naive", labels=NULL) {
  
  if(is.list(est) & (!is.data.frame(est))) {
    est_df <- est$est_df
    R2 <- round(c(0, sort(est$R2)), 2)
  } else if (is.data.frame(est)){
    est_df <- est
    R2 <- as.numeric(sub('...', '', colnames(est))) %>%
      sort() %>% round(2)
  }


  case <- 1:nrow(est_df)
  if(is.null(labels))
    labels <- case

  ## order according to naive estimates
  ord <- case
  if(order == "naive") {
    ord <- order(est_df[, "R2_0"])

  } else if(order == "worstcase") {
    rnk <- 1:length(est$rv)
    rob_neg <- which(is.na(est$rv) & est_df[, "R2_1_upr"] < 0)
    rnk[rob_neg] <- rank(est_df[rob_neg, "R2_1_upr"])

    not_rob <- which(!is.na(est$rv))
    rnk[not_rob] <- rank(est_df[not_rob, "R2_0"]) + length(rob_neg)

    rob_pos <- which(is.na(est$rv) & est_df[, "R2_1_lwr"] > 0)
    rnk[rob_pos] <- rank(est_df[rob_pos, "R2_1_lwr"]) + length(c(rob_neg, not_rob))
    ord <- order(rnk)

  }
  est_df <- est_df[ord, ]
  est$rv <- est$rv[ord]
  labels <- labels[ord]

  bound_df <- data.frame(x1 = case,
                         y1 = apply(est_df, 1, min),
                         x2 = case,
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

  plot <- ggplot(df) +
    ungeviz::geom_hpline(aes(x = case, y = effect, col = Type), width = 0.025*length(case), size = 1.5)  +
    geom_segment(data = bound_df, aes(x=x1,y=y1,xend=x2,yend=y2), size = 0.3) +
    labs(y = "Causal Effect", x = "Treatment Contrast") +
    scale_x_continuous(breaks = unique(case), labels = labels,
                       limits = c(0.5, max(case) + 0.5)) +
    scale_colour_viridis_d(name = expression(R^2~":"),
                           labels = paste0(R2*100, "%"), direction=-1) +
    theme_bw(base_size = 12)

  # add robustness value #
  if(!is.null(est$rv) & show_rv) {
    est$rv[!is.na(est$rv)] <- paste0(round(est$rv[!is.na(est$rv)]), "%")
    est$rv[is.na(est$rv)] <- "R"

    plot + annotate(geom = "text", x = case + 0.4, y = est_df[,'R2_0'],
                    size = 3, label = est$rv, col="dark red")
  } else {
    plot
  }
}



