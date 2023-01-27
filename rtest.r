library(Rcpp)
library(RcppArmadillo)

sourceCpp("ctest.cpp")

means <- list(c(1, 1), c(10, 10))
covs <- list(diag(c(1, 1)), diag(c(2, 2)))

params <- list("means" = means, "covs" = covs, "num_steps" = 10, "step_size"=0.05, "alpha"=1)

compute_mixture_pdf(c(0, 0), params, c(0.5, 0.5))
0.5 * dmvnorm(c(0, 0), c(1, 1), diag(c(1, 1))) + 0.5 * dmvnorm(c(0, 0), c(10, 10), diag(c(2, 2)))

params <- list("means" = means, "covs" = covs, "num_steps" = 200, "step_size"=0.3, "alpha"=1.04)
path <- get_path(x=c(-0.4, -0.9), p = c(0.7, -0.9), params, weights = c(0.5, 0.5))

# turn path into df
path <- data.frame(x = unlist(path)[seq(1, 600, 3)], y = unlist(path)[seq(2, 600, 3)], temp = unlist(path)[seq(3, 600, 3)])

path %>%
    ggplot(aes(x = x, y = y, color = temp)) + geom_path()

params <- list("means" = means, "covs" = covs, "num_steps" = 200, "step_size"=0.3, "alpha"=1.05)
another_path <- get_path(c(5, 5), c(-1, 1), params, c(0.5, 0.5))

another_path <- data.frame(x = unlist(another_path)[seq(1, 600, 3)], y = unlist(another_path)[seq(2, 600, 3)], temp = unlist(another_path)[seq(3, 600, 3)])

another_path %>%
    ggplot(aes(x = x, y = y, color = temp)) + geom_path()


samples <- get_samples(c(2, 2), params, c(1/2, 1/2), 10000)
samples <- unlist(samples)
x <- samples[seq(1, 2 * 10000 - 1, 2)]
y <- samples[seq(2, 2 * 10000, 2)]

df <- as.data.frame(cbind(x, y))
df %>%
    ggplot(aes(x = x, y = y)) + geom_point()

T <- 16
cov1 <- diag(c(1, 1)) / T
cov2 <- diag(c(2, 2)) / T
covs <- list(cov1, cov2)

means <- list(c(1, 1), c(10, 10))
weights <- c(0.5, 0.5)

params <- list("means" = means, "covs" = covs,
               "num_steps" = 200,
               "step_size" = 1.5,
               "alpha" = 1.0,
               "mass" = diag(c(3, 3)))

samples <- get_samples(c(2, 2), params, c(1/2, 1/2), 10000)
samples <- unlist(samples)
x <- samples[seq(1, 2 * 10000 - 1, 2)]
y <- samples[seq(2, 2 * 10000, 2)]

df <- as.data.frame(cbind(x, y))

df %>% mutate(temp = T) -> df3

df %>%
    ggplot(aes(x = x, y = y)) + geom_point()
