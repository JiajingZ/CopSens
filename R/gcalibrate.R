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
#' @param R2_constr an optional scalar or vector specifying the upper limit constraint on \eqn{R^2} .
#' By default, \code{R2_constr = 1}.
#' @param ... further arguments passed to \code{\link{pcaMethods::kEstimate}}, \code{\link{pcaMethods::pca}} or
#' \code{\link{get_opt_gamma}}.
#'
#' @return \code{gcalibrate} returns a list containing the following components:
#' \describe{
#'   \item{\code{est_df}}{a \code{data.frame} with naive and calibrated estimates of average treatment effects.}
#'   \item{\code{R2}}{a vector of \eqn{R^2} with elements corresponding to columns of \code{est_df}.}
#'   \item{\code{gamma}}{a matrix returned when \code{calitype = "multicali"} or \code{"worstcase"}.
#'   If \code{calitype = "multicali"}, optimized gamma are in columns,
#'   respectively resulting in estimates in columns of \code{est_df}.
#'   If \code{calitype = "worstcase"}, gamma are in rows,
#'   which respectively lead to the worstcase ignorance region for each contrast of interest.}
#'   \item{\code{rv}}{a \code{numeric vector} returned when \code{calitype = "worstcase"},
#'   with elements being the robustness value or \code{NA} if the ignorance region doesn't
#'   contains 0 for each contrast of interest.}
#' }
#
#' @export
#'
#' @examples
#' # load the example data #
#' y <- GaussianT_GaussianY$y
#' tr <- subset(GaussianT_GaussianY, select = -c(y))
#'
#' # worst-case calibration #
#' est_g1 <- gcalibrate(y = y, tr = tr, t1 = tr[1:2,], t2 = tr[3:4,],
#'                      calitype = "worstcase", R2 = c(0.3, 1))
#' plot_estimates(est_g1)
#'
#' # multivariate calibration #
#' est_g2 <- gcalibrate(y = y, tr = tr, t1 = tr[1:10,], t2 = tr[11:20,],
#'                      calitype = "multicali", R2_constr = c(1, 0.15))
#' plot_estimates(est_g2)
#'
#' # user-specified calibration #
#' est_g3 <- gcalibrate(y = y, tr = tr, t1 = tr[1:2,], t2 = tr[3:4,],
#'                      calitype = "null", gamma = c(0.96, -0.29, 0),
#'                      R2 = c(0.2, 0.6, 1))
#' plot_estimates(est_g3)
#' # apply gamma that maximizes the bias for the first contrast considered in est_g1 #
#' est_g4 <- gcalibrate(y = y, tr = tr, t1 = tr[1:2,], t2 = tr[3:4,],
#'                      calitype = "null", gamma = est_g1$gamma[1,],
#'                      R2 = c(0.2, 0.6, 1))
#' plot_estimates(est_g4)

gcalibrate <- function(y, tr, t1, t2, calitype = c("worstcase", "multicali", "null"),
                      mu_y_dt = NULL, sigma_y_t = NULL,
                      mu_u_dt = NULL, cov_u_t = NULL,
                      R2 = 1, gamma = NULL,
                      R2_constr = 1, ...) {
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
    bias <- sigma_y_t * apply(mu_u_dt %*% cov_halfinv, 1, function(x) sqrt(sum(x^2))) %o%
            c(0, rep(sqrt(R2), each = 2)*rep(c(-1, 1), times = length(R2)))
    est_df <- matrix(rep(mu_y_dt, times = 2*length(R2)+1), nrow = length(mu_y_dt)) + bias
    colnames(est_df) <- paste0("R2_", c(0, paste0(rep(R2, each = 2),
                                                   rep(c('_lwr', '_upr'), times = length(R2)))))
    gamma <- t(cov_halfinv %*%
      apply(mu_u_dt %*% cov_halfinv , 1, function(x) sigma_y_t*x/sqrt(sum(x^2))))
    rv_init <- (c(mu_y_dt^2) / apply(mu_u_dt %*% cov_halfinv, 1, function(x) sum(x^2)) /
             sigma_y_t^2) %>% round(digits = 4)
    rv <- rv_init * 100
    #rv[rv_init <= 1] <- paste0(rv_init[rv_init <= 1]*100, "%")
    rv[rv > 100] <- NA # "robust"
    list(est_df = data.frame(est_df), R2 = R2, gamma = gamma, rv = rv)
  } else if (calitype == "multicali" | calitype == "null") {
    if (calitype == "multicali") {
      message("Multivariate calibration executed.\n")
      cali <- matrix(NA, nrow = length(mu_y_dt), ncol = length(R2_constr))
      R2 <- rep(NA, length(R2_constr))
      gamma_mat <- matrix(NA, nrow = ncol(mu_u_dt), ncol = length(R2_constr))
      cat("Calibrating with R2_constr = ")
      for (i in 1:length(R2_constr)) {
        cat(R2_constr[i], " ")
        gamma_opt <- get_opt_gamma(mu_y_dt, mu_u_dt, cov_u_t, sigma_y_t,
                                   R2_constr = R2_constr[i], ...)
        cali[,i] <- mu_y_dt - mu_u_dt %*% gamma_opt
        R2[i] <- t(gamma_opt) %*% cov_u_t %*% gamma_opt / sigma_y_t^2
        gamma_mat[,i] <- gamma_opt
      }
      cat("\n")
      est_df <- data.frame(cbind(mu_y_dt, cali))
      colnames(est_df) <- paste0("R2_", round(c(0, R2), digits = 2))
      colnames(gamma_mat) <- colnames(est_df)[-1]
      list(est_df = est_df, R2 = R2, gamma = gamma_mat)
    } else if (calitype == "null" & is.null(gamma) == FALSE) {
      # eq (33) in terms of d = gamma #
      message("User-specified calibration executed.")
      cali <- matrix(NA, nrow = length(mu_y_dt), ncol = length(R2))
      for (i in 1:length(R2)) {
        cali[,i] <- mu_y_dt - sqrt(R2[i]) * sigma_y_t * mu_u_dt %*% cov_halfinv %*% gamma
      }
      cat("\n")
      est_df <- data.frame(cbind(mu_y_dt, cali))
      colnames(est_df) <- paste0("R2_", round(c(0, R2), digits = 2))
      list(est_df = est_df, R2 = R2)
    }
  } else {
    stop("Please specify a valid calibration type or gamma.")
  }
}






