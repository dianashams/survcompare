#' Calculates time-dependent Brier Score
#' @description
#' Calculates time-dependent Brier Scores for a vector of times. Calculations are similar to that in:
#' https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.brier_score.html#sksurv.metrics.brier_score
#' https://github.com/sebp/scikit-survival/blob/v0.19.0.post1/sksurv/metrics.py#L524-L644
#' The function uses IPCW (inverse probability of censoring weights), computed using the Kaplan-Meier
#' survival function, where events are censored events from train data
#'
#' @param y_predicted_newdata computed event probabilities
#' @param df_brier_train train data
#' @param df_newdata test data for which brier score is computed
#' @param time_points times at which BS calculated
#' @param weighted TRUE/FALSE for IPWC to use or not
#' @return vector of time-dependent Brier Scores for all time_points
surv_brierscore <-
  function(y_predicted_newdata,
           df_brier_train,
           df_newdata,
           time_points,
           weighted = TRUE) {
    # compute K-M probabilities of censoring for each observation till its individual time
    df_newdata$p_km_obs <-
      survival_prob_km(df_brier_train, df_newdata$time, estimate_censoring = TRUE)
    df_newdata$p_km_obs <-
      pmax(pmin(df_newdata$p_km_obs, 0.9999), 0.0001)

    # ! impute with mean observations if can't estimate !
    df_newdata[is.na(df_newdata$p_km_obs), "p_km_obs"] <-
      mean(df_newdata$p_km_obs, na.rm = 1)

    p_km_t <-
      survival_prob_km(df_brier_train, time_points, estimate_censoring = TRUE)
    p_km_t <- pmax(pmin(p_km_t, 0.9999), 0.0001)
    p_km_t

    bs <- c()
    for (j in 1:length(time_points)) {
      # assign t_point and predicted  probabilities
      if (length(time_points) == 1) {
        # only 1 time
        t_point <- time_points
        ppp <- pmax(pmin(y_predicted_newdata, 0.99999), 0.00001)
      } else {
        # many times
        t_point <- time_points[j]
        ppp <-
          pmax(pmin(y_predicted_newdata[, j], 0.99999), 0.00001)
      }
      # cases and controls by time t_point
      id_case <-
        ((df_newdata$time <= t_point) & (df_newdata$event == 1))
      id_control <- (df_newdata$time > t_point)

      # compute BS with weights which are 1/G(t) for controls and 1/G(obs_i) for cases
      # if weights == false, use w=1 for all
      if (weighted == TRUE) {
        # brier score is average of weighted squared errors
        bs[j] <-
          (
            sum(
              as.numeric(id_case) * (1 - ppp) ^ 2 * 1 / df_newdata$p_km_obs,
              na.rm = 1
            ) +
              sum(id_control * (0 - ppp) ^ 2 * 1 / p_km_t[j], na.rm = 1)
          ) / dim(df_newdata)[1]
      } else {
        # unweighted BS
        bs[j] <-
          sum(id_case * (1 - ppp) ^ 2 + id_control * (0 - ppp) ^ 2, na.rm = 1) / dim(df_newdata)[1]
      }
    }
    names(bs) <- round(time_points, 6)
    return(bs)
  }


#' Calculates survival probability estimated by Kaplan-Meier survival curve
#' Uses polynomial extrapolation in survival function space, using poly(n=3)
#' @param df_km_train event probabilities (!not survival)
#' @param times times at which survival is estimated
#' @param estimate_censoring FALSE by default, if TRUE, event and censoring is reversed (for IPCW calculations)
#' @return vector of survival probabilities for time_points
survival_prob_km <-
  function(df_km_train, times, estimate_censoring = FALSE) {
    if (estimate_censoring == FALSE) {
      km <-
        survival::survfit(survival::Surv(time, event) ~ 1, data = df_km_train)
      kmf <- stats::approxfun(km$time, km$surv, method = "constant")
      kmdf <- data.frame(cbind("time" = km$time, "surv" = km$surv))
      extrap <-
        stats::lm(surv ~ poly(time, 3, raw = TRUE), data = kmdf)
      km_extrap <- function(x) {
        (cbind(1, x, x ** 2, x ** 3) %*% extrap$coefficients[1:4])
      }
    } else {
      df_km_train$censor_as_event <- 1 - df_km_train$event
      km <-
        survival::survfit(survival::Surv(time, censor_as_event) ~ 1, data = df_km_train)
      kmdf <- data.frame(cbind("time" = km$time, "surv" = km$surv))
      kmf <- stats::approxfun(km$time, km$surv, method = "constant")
      extrap <-
        stats::lm(surv ~ poly(time, 3, raw = TRUE), data = kmdf)
      km_extrap <- function(x) {
        (cbind(1, x, x ** 2, x ** 3) %*% extrap$coefficients[1:4])
      }
    }
    return(km_extrap(times))
  }


#' Computes performance statistics for a survival data given the predicted event probabilities
#'
#' @param y_predict probabilities of event by predict_time (matrix=observations x times)
#' @param predict_time times for which event probabilities are given
#' @param df_train train data, data frame
#' @param df_test test data, data frame
#' @param weighted TRUE/FALSE, for IPWC
#' @param alpha calibration alpha as mean difference or from logistic regression
#' @return  data.frame(T, AUCROC, Brier Score, Scaled Brier Score, C_score, Calib slope, Calib alpha)
#' @export
surv_validate <- function(y_predict,
                          predict_time,
                          df_train,
                          df_test,
                          weighted = TRUE,
                          alpha = "logit") {
  # The function computes auc, brier score, c-index,
  # calibration slope and alpha for df_test
  # for apparent statistics use test  = train
  auc_score <- c()
  brier_score <- c()
  brier_score_scaled <- c()
  c_score <- c()
  calibration_slope <- c()
  calibration_alpha <- c()

  for (i in 1:length(predict_time)) {
    t_i <- predict_time[i]
    if (length(predict_time) > 1) {
      y_hat <- y_predict[, i]
    } else {
      y_hat <- unlist(y_predict)
    }

    temp <-
      try(timeROC::timeROC(
        T = df_test$time,
        delta = df_test$event,
        marker = y_hat,
        times = t_i * 0.9999999,
        cause = 1,
        weighting = "marginal"
      ),
      silent = TRUE)
    # time dependent AUC
    auc_score[i] <-
      ifelse(inherits(temp, "try-error"), NaN, temp$AUC[2])

    # compute time-dependent Brier score:
    temp <-
      try(surv_brierscore(y_hat, df_train, df_test, t_i, weighted = weighted),
          silent = TRUE)
    brier_score_scaled[i] <- NaN
    brier_score[i] <- NaN

    if (!inherits(temp, "try-error")) {
      brier_score[i] <- temp
      bs_base <-
        surv_brierscore(rep(mean((df_test$event) * (df_test$time <= t_i)
        ), dim(df_test)[1]),
        df_train, df_test, t_i, weighted = weighted)
      brier_score_scaled[i] <- 1 - brier_score[i] / bs_base
    }

    remove(temp)

    # compute concordance - time-dependent in a sense that a predictor is
    # event probability at t_i. For Cox model it is the same for each time
    # (probabilities are always ordered according to linear predictors for any t)
    temp <-
      try(survival::concordancefit(survival::Surv(df_test$time, df_test$event),-1 * y_hat),
          silent = TRUE)

    c_score[i] <-
      ifelse((inherits(temp, "try-error")) |
               is.null(temp$concordance),
             NaN,
             temp$concordance)

    # compute calibration slope and alpha:
    # 1/0 by t_i:
    df_test$event_ti <-
      ifelse(df_test$time <= t_i & df_test$event == 1, 1, 0)

    # cut 0 and 1 predicted probabilities for the logit to work:
    df_test$predict_ti <- pmax(pmin(y_hat, 0.99999), 0.00001)

    # Excluding censored observations before t_i, leave those with known state
    df_test_in_scope <-
      df_test[(df_test$time >= t_i) |
                (df_test$time < t_i & df_test$event == 1),]

    # Calibration slope and alpha.
    y_hat_hat <-
      log(df_test_in_scope$predict_ti / (1 - df_test_in_scope$predict_ti))
    y_actual_i <- df_test_in_scope$event_ti

    temp <-
      try(stats::glm(y_actual_i ~ y_hat_hat, family = binomial(link = "logit")),
          silent = TRUE)

    if (inherits(temp, "try-error")) {
      calibration_slope[i] <- NaN
      calibration_alpha[i] <- NaN
    } else {
      calibration_slope[i] <- temp$coefficients[2]
      if (alpha == "logit") {
        # take alpha from alpha: logit(y)~ logit(y_hat) + alpha
        calibration_alpha[i] <-
          stats::glm(y_actual_i ~ offset(y_hat_hat),
                     family = binomial(link = "logit"))$coefficients[1]
      } else {
        # take alpha as alpha= mean(y) - mean(y_hat)
        calibration_alpha[i] <-
          mean(y_actual_i) - mean(df_test_in_scope$predict_ti)
      }
    } # end "else"
  } # end "for"

  output <- data.frame(
    "T" = predict_time,
    "AUCROC" = auc_score,
    "BS" = brier_score,
    "BS_scaled" = brier_score_scaled,
    "C_score" = c_score,
    "Calib_slope" = calibration_slope,
    "Calib_alpha" = calibration_alpha
  )
  return(output)
}


#' Calibration stats of a fitted Cox PH model
#' @description
#' Computes calibration alpha and slope for a fitted coxph model
#' in the data.
#'
#' Crowson, C. S., Atkinson, E. J., & Therneau, T. M. (2016).
#' Assessing calibration of prognostic risk scores.
#' Statistical methods in medical research, 25(4), 1692-1706.
#'
#' https://journals.sagepub.com/doi/pdf/10.1177/0962280213497434
#'
#' @param cox_model fitted cox model, namely, coxph() object
#' @param test_data test data, should be a data frame with "time" and "event" columns for survival outcome
#' @return c(calibration alpha, calibration slope)
#' @export
cox_calibration_stats <-  function(cox_model, test_data) {
  if (!inherits(cox_model, "coxph")) {
    stop("The model should be a coxph object.")
  }
  if (!inherits(test_data, "data.frame")) {
    stop("The test data should be a dataframe.")
  }

  temp <-
    try(predict(cox_model, newdata = test_data, type = "lp"), silent = TRUE)

  if (inherits(temp, "try-error")) {
    stop ("Predictions can not be made for test data using the model provided.")
  }

  p <-
    log(predict(cox_model, newdata = test_data, type = "expected"))
  lp <- predict(cox_model, newdata = test_data, type = "lp")
  logbase <- p - lp

  fit1 <-
    try(stats::glm(event ~ offset(p), family = poisson, data = test_data),
        silent = TRUE)
  fit2 <-
    try(stats::glm(event ~ lp + offset(logbase),
                   family = poisson,
                   data = test_data),
        silent = TRUE)

  if (inherits(fit1, "try-error") |
      (inherits(fit2, "try-error"))) {
    stop("Stats computations failed.")
  }

  #group <- cut(lp, c(-Inf, quantile(lp, (1:9) / 10), Inf))
  #fit3 <- stats::glm(event ~ -1 + group + offset(p),family = poisson,data = test_data)

  calib_alpha <- as.numeric(fit1$coefficients[1])
  calib_slope <- as.numeric(fit2$coefficients[2])

  return(c("calib_alpha" = calib_alpha, "calib_slope" = calib_slope))
}


eligible_params <- function(params, df) {
  # This function checks eligible predictors from params list for split
  # It deletes those which are
  # 1) not in df and
  # 2) taking only 1 value (constants)
  # TODOLater may delete collinear factors
  if (length(params) == 0) {
    return(NULL)
  }
  # take only columns which are in df
  z <- params %in% names(df)
  if (sum(!z) == length(params)) {
    return(NULL) # no eligible params
  } else {
    params <- params[z] # there are some potentially eligible
  }
  params_eligible <- params
  for (i in 1:length(params)) {
    if (length(unique(df[, params[i]])) < 2) {
      params_eligible <- params_eligible[params_eligible != params[i]]
    }
  }
  return(params_eligible)
}


check_call <- function(inputs, inputclass, call_to_check) {
  #checks if inputs correspond to the right class
  indx <-
    match(names(call_to_check), names(inputclass), nomatch = 0)
  indx <- indx[indx != 0]
  classok <- rep(0, length(inputs))
  names(classok) = names(inputclass)
  classok[-indx] <- 1 #defaults are ok
  msg <- vector("character", length(inputs))
  for (i in indx) {
    obj <- inputs[[i]]
    class_required <- inputclass[[i]]
    if (class(obj) == class_required) {
      classok[i] = 1
    } else{
      classok[i] = 0
      msg[i] <-
        paste(
          names(inputclass)[i],
          " should be of class ",
          class_required,
          " (",
          class(obj),
          " is supplied).",
          sep = ""
        )
    }
  }
  anyerror <- any(classok==0)
  return(list("anyerror"=anyerror, "classok" = classok, "msg"= msg))
}


