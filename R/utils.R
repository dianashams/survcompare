#' Calculates time-dependent Brier Score
#' @description
#' Calculates time-dependent Brier Scores for a vector of times. Calculations are similar to that in:
#' https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.brier_score.html#sksurv.metrics.brier_score
#' https://github.com/sebp/scikit-survival/blob/v0.19.0.post1/sksurv/metrics.py#L524-L644
#' The function uses IPCW (inverse probability of censoring weights), computed using the Kaplan-Meier
#' survival function, where events are censored events from train data
#' @param y_predicted_newdata computed event probabilities (! not survival probabilities)
#' @param df_brier_train train data
#' @param df_newdata test data for which brier score is computed
#' @param time_point times at which BS calculated
#' @param weighted TRUE/FALSE for IPWC to use or not
#' @return vector of time-dependent Brier Scores for all time_point
surv_brierscore <-
  function(y_predicted_newdata,
           df_brier_train,
           df_newdata,
           time_point,
           weighted = TRUE) {

    cut01 = function(x) {return(pmax(pmin(x,0.9999),0.0001))}
    # compute K-M probabilities of censoring for each observation till its individual time
    if (weighted) {
      df_newdata$p_km_obs <-
        cut01(survival_prob_km(df_brier_train, df_newdata$time, estimate_censoring = TRUE))
      p_km_t <-
        cut01(survival_prob_km(df_brier_train, time_point, estimate_censoring = TRUE))
    }else{
      df_newdata$p_km_obs = 1
      p_km_t=1
      }
    ppp <- cut01(y_predicted_newdata)
    #! impute with mean observations if can't estimate !
    df_newdata[is.na(df_newdata$p_km_obs), "p_km_obs"] <- mean(df_newdata$p_km_obs, na.rm = TRUE)

    # cases and controls by time_point
    id_case <- ((df_newdata$time <= time_point) & (df_newdata$event == 1))
    id_control <- (df_newdata$time > time_point)|
      ((df_newdata$time == time_point)&(df_newdata$event == 0))

    # compute BS with weights which are 1/G(t) for controls and 1/G(obs_i) for cases
    bs <-
      (sum(as.numeric(id_case) * (1 - ppp) ^ 2 * 1 / df_newdata$p_km_obs,na.rm = FALSE) +
         sum(id_control * (0 - ppp) ^ 2 * 1 / p_km_t, na.rm = FALSE)
       ) / dim(df_newdata)[1]
    return(bs)
}


#' Calculates survival probability estimated by Kaplan-Meier survival curve
#' Uses polynomial extrapolation in survival function space, using poly(n=3)
#' @param df_km_train event probabilities (!not survival)
#' @param times times at which survival is estimated
#' @param estimate_censoring FALSE by default, if TRUE, event and censoring are reversed (for IPCW calculations)
#' @return vector of survival probabilities for time_point
survival_prob_km <-
  function(df_km_train, times, estimate_censoring = FALSE) {
    if (estimate_censoring == FALSE) {
      km <- survival::survfit(survival::Surv(time, event) ~ 1, data = df_km_train)
    } else {
      times[times==max(df_km_train$time)] = 0.9999*max(df_km_train$time)
      km <- survival::survfit(survival::Surv(time, 1- event) ~ 1, data = df_km_train)
    }
    kmf <- stats::approxfun(km$time, km$surv, method = "constant")
    return(kmf(times))
  }


#' Computes performance statistics for a survival data given the predicted event probabilities
#'
#' @param y_predict probabilities of event by predict_time (matrix=observations x times)
#' @param predict_time times for which event probabilities are given
#' @param df_train train data, data frame
#' @param df_test test data, data frame
#' @param weighted TRUE/FALSE, for IPWC
#' @param alpha calibration alpha as mean difference in probabilities, or in log-odds (from logistic regression, default)
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

  #if all NaNs, return NaNs
  if (sum(is.na(y_predict))==length(y_predict)){
    output <- data.frame(
      "T" = predict_time,
      "AUCROC" = NaN,
      "BS" = NaN,
      "BS_scaled" = NaN,
      "C_score" = NaN,
      "Calib_slope" = NaN,
      "Calib_alpha" = NaN
    )
    return(output)
  }

  # 1) Concordance
  # time-dependent in a sense that a predictor is event prob @ predict_time.
  # For Cox model it is the same for each time
  # (probabilities are always ordered as linear predictors for each t)
  temp <- try(survival::concordancefit(
      survival::Surv(df_test$time, df_test$event),-1 * y_predict),
      silent = TRUE)
  c_score <- ifelse((inherits(temp, "try-error")) |
                           inherits(try(temp$concordance, silent = TRUE), "try-error"),
                         NaN, temp$concordance)

  # 2) time dependent AUC
  # this gives NA for the final time, so we move it to t-
  predict_time_auc <- predict_time
  if (predict_time == max(df_test$time)) {predict_time_auc <- 0.9999 * max(df_test$time)}
  temp <-  try(timeROC::timeROC(
      T = df_test$time,delta = df_test$event,
      marker = y_predict,times = predict_time_auc,cause = 1),
      silent = TRUE)
  auc_score <- ifelse(inherits(temp, "try-error"), NaN, temp$AUC[2])

  # 3) time-dependent Brier score:
  temp <- try(
    surv_brierscore(y_predict, df_train, df_test, predict_time, weighted = weighted),
    silent = TRUE)
  brier_score_scaled <- NaN
  brier_score <- NaN
  if (!inherits(temp, "try-error")) {
    brier_score <- temp
    bs_base <-
      surv_brierscore(
        y_predicted_newdata =
          rep(mean((df_test$event) * (df_test$time <= predict_time)),
              dim(df_test)[1]),
        df_brier_train = df_train,
        df_newdata = df_test,
        time_point = predict_time,
        weighted = weighted
      )
    brier_score_scaled <- 1 - brier_score / bs_base
  }
  # 4) Calibration slope and alpha:
  # 1/0 by predict_time:
  df_test$event_ti <- ifelse(df_test$time <= predict_time & df_test$event == 1, 1, 0)
  # cut 0 and 1 predicted probabilities for the logit to work:
  df_test$predict_ti <- pmax(pmin(y_predict, 0.9999), 0.0001)

  # Exclude censored observations before predict_time, leave those with known state
  df_test_in_scope <-
    df_test[(df_test$time >= predict_time) |
              (df_test$time < predict_time & df_test$event == 1), ]

  # 4) Calibration slope and alpha.
  y_predict_hat <-
    log(df_test_in_scope$predict_ti / (1 - df_test_in_scope$predict_ti))
  y_actual_i <- df_test_in_scope$event_ti
  temp <- try(stats::glm(y_actual_i ~ y_predict_hat,
                         family = binomial(link = "logit")),
              silent = TRUE)
  calibration_slope <- NaN
  calibration_alpha <- NaN
  if (!inherits(temp, "try-error")) {
    calibration_slope <- temp$coefficients[2]
    if (alpha == "logit") {
      # take alpha from alpha: logit(y)~ logit(y_predict) + alpha
      temp2 <- try(stats::glm(y_actual_i ~ y_predict_hat,
                             family = binomial(link = "logit")),
                   silent=TRUE)
      if(!inherits(temp2, "try-error")){
        calibration_alpha <-
        stats::glm(y_actual_i ~ offset(y_predict_hat),
                   family = binomial(link = "logit"))$coefficients[1]}
    } else {
      # take alpha as alpha= mean(y) - mean(y_predict)
      calibration_alpha <-
        mean(y_actual_i) - mean(df_test_in_scope$predict_ti)
    }# end "else"
  } # end if try-error

  remove(temp)
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


