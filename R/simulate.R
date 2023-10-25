
# These functions return simulated data with exponential or weibull distributions
# 0.5 / 0.75 / 0.9 - by time = 10 (scaled) expected number of events
# distr = "Exp" or "Weibull"
# (rho_w =0.5, lambda = 0.447) hazard slopes down
# (rho_w=1.5, lambda 0.027) hazard up and down;
# drop-out - additional drop out before the end of study (expected, independent from event)

#' Simulated sample with survival outcomes with linear dependencies
#' @description
#' Simulated sample with exponentially or Weibull distributed time-to-event;
#' log-hazard (lambda parameter) depends linearly on risk factors.
#'
#' @param N sample size, 300 by default
#' @param observe_time study's observation time, 10 by default
#' @param percentcensored expected number of non-events by observe_time, 0.75 by default (i.e. event rate is 0.25)
#' @param drop_out expected rate of drop out before observe_time, 0.3 by default
#' @param randomseed random seed for replication
#' @param distr time-to-event distribution, "Exp" for exponential (default), "W" for Weibull
#' @param lambda baseline hazard rate, 0.1 by default
#' @param rho_w shape parameter for Weibull distribution, 0.3 by default
#' @examples
#' mydata <- simulate_linear()
#' head(mydata)
#' @return data frame; "time" and "event" columns describe survival outcome; predictors are "age", "sex", "hyp", "bmi"
#' @export
simulate_linear <- function(N = 300, observe_time = 10,
                                percentcensored = 0.75,
                                randomseed = NULL, lambda = 0.1,
                                distr = "Exp", rho_w = 1,
                                drop_out = 0.3) {

  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}

  # simulate the data
  df <- simulate_population(N, randomseed)
  # calculate betas
  exp_beta <- linear_beta(df)
  df["exp_beta"] <- exp_beta
  # simulate censored and event times
  df <- df_event_times(
    exp_beta = exp_beta, df = df, N = N, observe_time = observe_time,
    percentcensored = percentcensored, randomseed = randomseed + 1,
    lambda = lambda, distr = distr, rho_w = rho_w, drop_out = drop_out
  )
  return(df)
}

#' Simulated sample with survival outcomes with non-linear dependencies
#' @description
#' Simulated sample with exponentially or Weibull distributed time-to-event;
#' log-hazard (lambda parameter) depends non-linearly on risk factors.
#' @param N sample size, 300 by default
#' @param observe_time study's observation time, 10 by default
#' @param percentcensored expected number of non-events by observe_time, 0.75 by default (i.e. event rate is 0.25)
#' @param drop_out expected rate of drop out before observe_time, 0.3 by default
#' @param randomseed random seed for replication
#' @param distr time-to-event distribution, "Exp" for exponential (default), "W" for Weibull
#' @param lambda baseline hazard rate, 0.1 by default
#' @param rho_w shape parameter for Weibull distribution, 0.3 by default
#' @examples
#' mydata <- simulate_nonlinear()
#' head(mydata)
#' @return data frame; "time" and "event" columns describe survival outcome; predictors are "age", "sex", "hyp", "bmi"
#' @export
simulate_nonlinear <- function(N = 300, observe_time = 10,
                                   percentcensored = 0.75,
                                   randomseed = NULL, lambda = 0.1,
                                   distr = "Exp", rho_w = 1,
                                   drop_out = 0.3) {
  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}
  # simulate the data
  df <- simulate_population(N, randomseed)
  # calculate betas
  exp_beta <- nonlinear_beta(df)
  # simulate censored and event times
  df <- df_event_times(
    exp_beta = exp_beta, df = df, N = N, observe_time = observe_time,
    percentcensored = percentcensored, randomseed = randomseed + 1,
    lambda = lambda, distr = distr, rho_w = rho_w, drop_out = drop_out
  )
  return(df)
}

#' Simulated sample with survival outcomes with non-linear and cross-term dependencies
#' @description
#' Simulated sample with exponentially or Weibull distributed time-to-event;
#' log-hazard depends non-linearly on risk factors, and includes cross-terms.
#' @param N sample size, 300 by default
#' @param observe_time study's observation time, 10 by default
#' @param percentcensored expected number of non-events by observe_time, 0.75 by default (i.e. event rate is 0.25)
#' @param drop_out expected rate of drop out before observe_time, 0.3 by default
#' @param randomseed random seed for replication
#' @param distr time-to-event distribution, "Exp" for exponential (default), "W" for Weibull
#' @param lambda baseline hazard rate, 0.1 by default
#' @param rho_w shape parameter for Weibull distribution, 0.3 by default
#' @examples
#' mydata <- simulate_crossterms()
#' head(mydata)
#' @return data frame; "time" and "event" columns describe survival outcome; predictors are "age", "sex", "hyp", "bmi"
#' @export
simulate_crossterms <- function(N = 300, observe_time = 10,
                                    percentcensored = 0.75,
                                    randomseed = NULL, lambda = 0.1,
                                    distr = "Exp", rho_w = 1,
                                    drop_out = 0.3) {
  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}
  # simulate the data
  df <- simulate_population(N, randomseed)
  # calculate betas
  exp_beta <- xt_beta(df)
  # simulate censored and event times
  df <- df_event_times(
    exp_beta = exp_beta, df = df, N = N, observe_time = observe_time,
    percentcensored = percentcensored, randomseed = randomseed + 1,
    lambda = lambda, distr = distr, rho_w = rho_w, drop_out = drop_out
  )
  return(df)
}

# Simulated sample with survival outcomes with non-linear and cross-term dependencies and non-PH survival times.
# @description
# Simulated sample with Weibull distributed time-to-event;
# log-hazard are non-linear and with cross-terms. Survival curves for different sexes intersect (non-PH).
#
# @param N sample size, 300 by default
# @param observe_time study's observation time, 10 by default
# @param percentcensored expected number of non-events by observe_time, 0.75 by default (i.e. event rate is 0.25)
# @param drop_out expected rate of drop out before observe_time, 0.3 by default
# @param randomseed random seed for replication
# @param lambda baseline hazard rate, 0.1 by default
# @examples
# mydata <- simulate_lin_nonPH()
# head(mydata)
# @return data frame; "time" and "event" columns describe survival outcome; predictors are "age", "sex", "hyp", "bmi"
# simsurv_lin_nonPH <- function(N = 300,
#                                    observe_time = 10,
#                                    percentcensored = 0.75,
#                                    randomseed = NULL,
#                                    lambda = 0.1,
#                                    drop_out = 0.3) {
#   # TODO check why Cox model is performing so well for these simulations
#   # simulate the data
#   if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}
#   df <- simulate_population(N, randomseed)
#   # calculate betas
#   exp_beta <-
#     exp(0.4 * df$age + 1.0 * df$bmi + 1 * df$hyp + 0.5 * df$sex)
#   df["exp_beta"] <- exp_beta
#   df["shape_rho"] <- ifelse(df$sex == 1, 2.5, 0.5)
#   # simulate event times
#   {
#     set.seed(randomseed)
#     v <- stats::runif(n = N)
#   }
#   # Weibull density
#   event_time <-
#     (-log(v) / (lambda * df$exp_beta))^(1 / df$shape_rho)
#
#   # re-scale the time to have 1-percentcensored of events=1
#   # by the observe_time
#   final_time <- quantile(event_time, 1 - percentcensored)
#   # scale to observe_time:
#   df$event_time <-
#     0.001 + pmin(round(event_time / final_time * observe_time, 3), observe_time)
#   # generate drop-out times for random drop_out % observations
#   if (drop_out > 0) {
#     set.seed(randomseed + 1)
#     randcentime <-
#       stats::runif(round(N * drop_out, 0), 0, observe_time)
#     cens_obs <- sample.int(nrow(df), round(N * drop_out, 0))
#     df[cens_obs, "cens_time"] <- randcentime
#     df[-cens_obs, "cens_time"] <- observe_time
#   } else {
#     df[, "cens_time"] <- observe_time
#   }
#
#   # final time and event definition
#   # event =1 if event time < cens_time and observe_time
#   df$time <- pmin(df$event_time, df$cens_time, observe_time)
#   df$event <- ifelse(df$event_time == df$time, 1, 0)
#
#   return(df)
# }


#' Auxiliary function for simulatedata functions
#' @param df data
linear_beta <- function(df) {
  return(exp(0.4 * df$age + 1.0 * df$bmi + 0.7 * df$hyp))
}

nonlinear_beta <- function(df) {
  # BMI impact is 2 for very low and high levels,
  # 1 for high/ low level, 0 for normal range
  bmi_beta <-
    ifelse((df$bmi < -1.5) |
             (df$bmi > 1.5), 2, ifelse((df$bmi < -1) | (df$bmi > 1), 1, 0))
  # Age impact is 1 for age>=55; linear age impact is also present,
  # but is smaller than in linear simulation
  age_beta <- ifelse((df$age >= 1), 1, 0)
  return(exp(bmi_beta + (df$hyp * 0.7) + df$age * 0.2 + age_beta))
}

xt_beta <- function(df) {
  # BMI impact is 2 for very low and high levels,
  # 1 for high/ low level, 0 for normal range
  bmi_beta <- ifelse((df$bmi < -1.5) | (df$bmi > 1.5), 2,
                     ifelse((df$bmi < -1) | (df$bmi > 1), 1, 0)
  )

  # hypertension x age interaction
  hyp_beta <- ifelse((df$age >= 1 & df$hyp == 1), 2,
                     ifelse((df$age < 1 & df$hyp == 1), 1, 0)
  )
  # simulating event time
  return(exp(bmi_beta + hyp_beta + df$age * 0.2))
}


simulate_population <- function(N = 300, randomseed = 42) {
  # Auxiliary function for simulatedata
  # Simulates age, bmi, hyp and sex for a generated population
  # param N sample size
  # param randomseed random seed, 42 by default
  # return data frame
  set.seed(randomseed)
  df <- data.frame(
    age = round(stats::runif(N, -1.73, 1.73), 1),
    bmi = round(stats::rnorm(N, 0, 1), 1),
    hyp = stats::rbinom(N, 1, 0.20),
    sex = stats::rbinom(N, 1, 0.5)
  )
  return(df)
}

# Auxiliary function for simsurve_linear, non_linear, crossterms
# param exp_beta exponential distribution beta param
# param df data
# param N sample size
# param observe_time observation time
# param percentcensored percent of censored observations (with no event) by observe_time
# param randomseed random seed
# param lambda exponential/Weibull distribution lambda param
# param distr distribution for times, "Exp" or "W" for exponential or Weibull
# param rho_w Rho parameter for Weibull distribution
# param drop_out drop out rate of those censored before observe_time
# return data frame with added columns "time", "event", "event_time"
df_event_times <- function(exp_beta, df, N, observe_time,
                           percentcensored, randomseed,
                           lambda, distr, rho_w, drop_out) {

  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}

  # simulate event times
  if (distr == "Exp") {
    # Exponential time distribution, h0(t)=??=0.1  with shape ??>0 and scale ??>0
    {
      set.seed(randomseed)
      event_time <- stats::rexp(N, lambda * exp_beta)
    }
  } else {
    # simulate Weibull distribution
    # h0(t)=lambda*rho*t^(rho???1), shape rho_w>0,  scale lambda>0.
    # If rho=1, it is exponential
    # rho>1 => upward sloping, rho<1 => downward, rho=1 constant (exponential)
    {
      set.seed(randomseed)
      v <- stats::runif(n = N)
    }
    # Weibull density
    event_time <- (-log(v) / (lambda * exp_beta))^(1 / rho_w)
  }

  # re-scale the time to have 1-percentcensored of events=1
  # by the observe_time
  final_time <- quantile(event_time, 1 - percentcensored)
  # scale to observe_time
  df$event_time <- 0.001 + pmin(round(event_time / final_time * observe_time, 3), observe_time)
  # generate drop-out times for random drop_out % observations
  if (drop_out > 0) {
    set.seed(randomseed + 1)
    randcentime <- stats::runif(round(N * drop_out, 0), 0, observe_time)
    cens_obs <- sample.int(nrow(df), round(N * drop_out, 0))
    df[cens_obs, "cens_time"] <- randcentime
    df[-cens_obs, "cens_time"] <- observe_time
  } else {
    df[, "cens_time"] <- observe_time
  }

  # final time and event definition
  # event =1 if event time < cens_time and observe_time
  df$time <- pmin(df$event_time, df$cens_time, observe_time)
  df$event <- ifelse(df$event_time == df$time, 1, 0)

  return(df)
}

