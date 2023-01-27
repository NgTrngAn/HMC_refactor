library(Rcpp)
library(RcppArmadillo)

sourceCpp("ctest.cpp")

means <- list(c(1, 1), c(10, 10))
covs <- list(diag(c(1, 1)), diag(c(2, 2)))

params <- list("means" = means, "covs" = covs, "num_steps" = 10, "step_size"=0.05, "alpha"=1)

compute_mixture_pdf(c(0, 0), params, c(0.5, 0.5))
0.5 * dmvnorm(c(0, 0), c(1, 1), diag(c(1, 1))) + 0.5 * dmvnorm(c(0, 0), c(10, 10), diag(c(2, 2)))

params <- list("means" = means, "covs" = covs, "num_steps" = 200, "step_size"=0.3, "alpha"=1.04)
path <- get_path(x=c(1, 1), p = c(1, -1), params, weights = c(0.5, 0.5))
