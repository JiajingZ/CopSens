# Generating data with Gaussian-Gaussian Linear Model #
k <- 10
s <- 3
set.seed(123)
B <- rstiefel::rustiefel(10, 3) %*% diag(c(2, -1, 0.5))
tau <- runif(k, min = -1, max = 1)
gamma <- 5 * svd(B)$v[, 1] + 1.5 * svd(B)$v[, 2]
sigma2_t <- 0.75
sigma2_y <- 0.2

set.seed(234)
n <- 10000
u <- MASS::mvrnorm(n = n, mu = rep(0, s), Sigma = diag(s))
colnames(u) <- paste0(rep('u', s), 1:s)
tr <- u%*%t(B) + MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = sigma2_t*diag(k))
colnames(tr) <- paste0(rep('t', k), 1:k)
y_gaussian <-  tr%*%tau + u%*%gamma + rnorm(n, mean = 0, sd = sqrt(sigma2_y))
y_binary <- ifelse(y_gaussian > 0, 1, 0)

GaussianT_GaussianY <- data.frame(y = y_gaussian, tr)
GaussianT_BinaryY <- data.frame(y = y_binary, tr)

# parameters value in theory #
Cov_u_t <- diag(s) - t(B) %*% solve(B %*% t(B) + sigma2_t * diag(k)) %*% B
sigma2_y_t <- t(gamma) %*% Cov_u_t %*% gamma + sigma2_y
R2 <- t(gamma) %*% Cov_u_t %*% gamma / sigma2_y_t

usethis::use_data(GaussianT_GaussianY, GaussianT_BinaryY, overwrite = TRUE)
# usethis::use_data(GaussianT_GaussianY, GaussianT_BinaryY, internal = TRUE, overwrite = TRUE)





