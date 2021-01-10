#' Obtain Optimized Sensitivity Parameters Using Multivariate Calibration Criterion
#'
#' @param mu_y_dt Scalar or vector that contains naive estimates of treatment effects
#' ignoring confounding.
#' @param mu_u_dt Matrix of difference in conditional confounder means, \eqn{E(U \mid t1) - E(U \mid t2)},
#' with latent variables in columns.
#' @param cov_u_t Covariance matrix of confounders conditional on treatments.
#' @param sigma_y_t Scalar of the standard deviation of outcome conditional on treatments.
#' @param penalty_weight scalar specifying the penalty weight when constraints on \eqn{R^2} are placed.
#' @param gamma0 An optional vector with initial values of the sensitivity parameters to be optimized over.
#' @param n_iter Number of times of optimization to execute.
#' @param normtype character. Optional function \code{m} for the multivariate calibration criterion.
#' By default, the L2 norm will be applied.\cr
#' "L1" - apply the L1 norm, \code{sum(abs(x))}. \cr
#' "L2" - apply the L2 norm, \code{sqrt(sum(x^2))}.\cr
#' "Inf" - apply the infinity norm, \code{max(abs(x))}. \cr
#'
#' @return Optimized sensitivity parameters.
#'

get_opt_gamma <- function(mu_y_dt, mu_u_dt, cov_u_t, sigma_y_t,
                          penalty_weight = 0, gamma0 = NULL, n_iter = 100,
                          normtype = "L2") {
  # Objective function for Optimization #
  objective <- function(gamma) {
    cali <- mu_y_dt - mu_u_dt %*% gamma
    R2 <- t(gamma) %*% cov_u_t %*% gamma / sigma_y_t^2
    if (normtype == "L1") {
      sum(abs(cali)) + penalty_weight * R2
    } else if (normtype == "L2") {
      sqrt(sum(cali^2)) + penalty_weight * R2
    } else if (normtype == "Inf") {
      max(abs(cali)) + penalty_weight * R2
    } else {stop("Please specify a norm type.")}
  }
  latent_dim <- ncol(mu_u_dt)
  if (is.null(gamma0)){
    gamma0 <- rep(0, latent_dim)
  }
  obj_min <- objective(gamma0)
  gamma_opt <- gamma0
  for (i in 1:n_iter) {
    solution <- optim(par = gamma0, fn = objective, method = "BFGS")
    if (solution$value < obj_min) {
      obj_min <- solution$value
      gamma_opt <- solution$par
    }
    gamma0 <- rnorm(latent_dim)
  }
  gamma_opt
}


