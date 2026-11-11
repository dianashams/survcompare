cox_test_data <- function() {
  set.seed(714)
  x1 <- stats::rnorm(120)
  x2 <- stats::rnorm(120)
  failure <- stats::rexp(120, exp(0.7 * x1 - 0.3 * x2))
  censor <- stats::rexp(120, 0.1)
  data.frame(x1 = x1, x2 = x2,
             time = pmin(failure, censor) + 0.01,
             event = as.numeric(failure <= censor))
}
