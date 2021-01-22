#' Calibration for Gaussian Outcomes
#'
#' @description Calibrates the naive estimates to account for unobserved confounding when outcome
#' variables are Gaussian. The calibration can be done with user-specific sensitivity parameters or
#' with our pre-provided calibration methods, the worst-case calibration for a single contrast
#' or multivariate calibration for multiple contrasts.
#'
#' @param y \code{data.frame}, \code{matrix} or \code{vector}. Gaussian outcome variable.
#' @param tr \code{data.frame}. Treatment variables with rows corresponding to observations and columns
#' to variables.
#' @param t1 \code{data.frame}. First treatment arms of interest.
#' May contain a single or multiple treatments in rows.
#' @param t2 \code{data.frame}. Second treatment arms of interest,
#' which has same number of row as \code{t1}.
#' @param calitype character. The calibration method to be applied. Can be one of: \cr
#' "worstcase" - apply worst-case calibration when considering a single contrast. \cr
#' "multicali" - apply mutlivariate calibration when considering multiple contrasts.\cr
#' "null" - apply calibration with user-specified sensitivity parameter, \eqn{\gamma}. \cr
#' @param mu_y_dt an optional scalar or vector that contains naive estimates of treatment effects
#' ignoring confounding.
#' @param sigma_y_t an optional scalar of the standard deviation of outcome conditional on treatments.
#' @param mu_u_dt an optional matrix of difference in conditional confounder means, \eqn{E(U \mid t1) - E(U \mid t2)},
#' with latent variables in columns.
#' @param cov_u_t an optional covariance matrix of confounders conditional on treatments.
#' @param R2 an optional scalar or vector specifying the proportion of residual variance in outcome given the
#' treatment that can be explained by confounders.
#' @param gamma sensitivity parameter vector. Must be given when \code{calitype = "null"}.
#' @param penalty_weight an optional scalar or vector specifying the penalty weight for \eqn{R^2} to place
#' constraints on its magnitude. By default, \code{penalty_weight = 0}.
#' @param ... further arguments passed to \code{\link{pcaMethods::kEstimate}}, \code{\link{pcaMethods::pca}} or
#' \code{\link{get_opt_gamma}}.
#'
#'
#' @return A \code{data.frame} with naive and calibrated estimates of average treatment effects.
#' @export
#'
#' @examples
#' # load the example data #
#' y <- GaussianT_GaussianY$y
#' tr <- subset(GaussianT_GaussianY, select = -c(y))
#' # worst-case calibration #
#' est_df1 <- gcalibrate(y = y, tr = tr, t1 = tr[1,], t2 = tr[2,],
#'                       calitype = "worstcase", R2 = c(0.6, 1))
#' plot_estimates(est_df1)
#' # multivariate calibration #
#' est_df2 <- gcalibrate(y = y, tr = tr, t1 = tr[1:10,], t2 = tr[11:20,],
#'                       calitype = "multicali")
#' plot_estimates(est_df2)
#' # user-specified calibration #
#' est_df3 <- gcalibrate(y = y, tr = tr, t1 = tr[1:2,], t2 = tr[3:4,],
#'                       calitype = "null", gamma = c(0.96, -0.29, 0),
#'                       R2 = c(0.3, 0.7, 1))
#' plot_estimates(est_df3)

gcalibrate <- function(y, tr, t1, t2, calitype = c("worstcase", "multicali", "null"),
                      mu_y_dt = NULL, sigma_y_t = NULL,
                      mu_u_dt = NULL, cov_u_t = NULL,
                      R2 = 1, gamma = NULL,
                      penalty_weight = 0, ...) {
  # by default, fitting latent confounder model by PPCA #
  if (is.null(mu_u_dt) | is.null(cov_u_t)) {
    message("Fitting the latent confounder model by PPCA with default.")
    ut_cv <- pcaMethods::kEstimate(tr, method = "ppca", allVariables = TRUE, ...)
    ut_ppca <- pcaMethods::pca(tr, method = "ppca", center = TRUE,
                               nPcs = ut_cv$bestNPcs, ...)
    W = pcaMethods::loadings(ut_ppca)
    tr_hat <- pcaMethods::scores(ut_ppca) %*% t(W)
    sig2est <- sum((tr - tr_hat)^2)/(nrow(tr)*ncol(tr))
    ## cov(U|t) = sigma2*M^{-1}, M = W'W+ sigma2*I
    if (is.null(cov_u_t)) {
      cov_u_t <- sig2est * solve(t(W) %*% W + sig2est*diag(ut_cv$bestNPcs))
    }
    if (is.null(mu_u_dt)) {
      mu_u_dt <- predict(ut_ppca, newdata = t1)$scores - predict(ut_ppca, newdata = t2)$scores
    }
  }
  if (ncol(mu_u_dt) == 1) {
    cov_halfinv <- 1/sqrt(cov_u_t)
  } else {
    eigen_cov <- eigen(cov_u_t)
    cov_halfinv <- eigen_cov$vectors %*% diag(eigen_cov$values^{-1/2}) %*% t(eigen_cov$vectors)
  }
  # by default, fitting the outcome by ordinary linear regression model #
  if (is.null(mu_y_dt) | is.null(sigma_y_t)) {
    message("Observed outcome model fitted by simple linear regression with default.")
    lm_y_t <- lm(y ~., data = data.frame(y,tr))
    if (is.null(mu_y_dt)) {
      mu_y_dt <- predict(lm_y_t, newdata = t1) - predict(lm_y_t, newdata = t2)
    }
    if (is.null(sigma_y_t)) {
      sigma_y_t <- sigma(lm_y_t)
    }
  }
  if (calitype == "worstcase") {
    message("Worst-case calibration executed.")
    bias <- sqrt(R2) * sigma_y_t * sqrt(sum((cov_halfinv %*% c(mu_u_dt))^2))
    results <- data.frame(cbind(rep(mu_y_dt, 2),
                     rbind(mu_y_dt - bias, mu_y_dt + bias)))
    colnames(results) <- paste0("R2_", c(0, R2))
    rownames(results) <- c("lower", "upper")
    results
  } else if (calitype == "multicali" | calitype == "null") {
    if (calitype == "multicali") {
      message("Multivariate calibration executed.")
      cali <- matrix(NA, nrow = length(mu_y_dt), ncol = length(penalty_weight))
      R2 <- rep(NA, length(penalty_weight))
      for (i in 1:length(penalty_weight)) {
        gamma_opt <- get_opt_gamma(mu_y_dt, mu_u_dt, cov_u_t, sigma_y_t,
                                   penalty_weight = penalty_weight[i], ...)
        cali[,i] <- mu_y_dt - mu_u_dt %*% gamma_opt
        R2[i] <- t(gamma_opt) %*% cov_u_t %*% gamma_opt / sigma_y_t^2
      }
    } else if (calitype == "null" & is.null(gamma) == FALSE) {
      # eq (33) in terms of d = gamma #
      message("User-specified calibration executed.")
      cali <- matrix(NA, nrow = length(mu_y_dt), ncol = length(R2))
      for (i in 1:length(R2)) {
        cali[,i] <- mu_y_dt - sqrt(R2[i]) * sigma_y_t * mu_u_dt %*% cov_halfinv %*% gamma
      }
    }
    cat("\n")
    results <- data.frame(cbind(mu_y_dt, cali))
    colnames(results) <- paste0("R2_", round(c(0, R2), digits = 2))
    results
  } else {
    stop("Please specify a valid calibration type or gamma.")
  }
}






