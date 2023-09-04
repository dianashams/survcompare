############# Basic Cox Model functions ##########

#' Trains CoxPH using survival package, or trains CoxLasso (cv.glmnet, lambda.min),
#'  and then re-trains survival:coxph on non-zero predictors
#'
#' @param df_train  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of the column names to be used as predictors
#' @param useCoxLasso  TRUE or FALSE
#' @param fixed_time  not used here, to re-align with other methods
#' @param retrain_cox if useCoxLasso is TRUE, whether to re-train coxph on non-zero predictors, FALSE by default
#' @examples
#' df<- simsurv_nonlinear()
#' @return fitted CoxPH or CoxLasso model
survcox_train <- function(df_train,
                             predict.factors,
                             useCoxLasso = FALSE,
                             fixed_time = NaN,
                             retrain_cox = FALSE) {
  stopifnot(expr = {
    is.data.frame(df_train)
    predict.factors %in% colnames(df_train)
  })

  # wrapper for coxph() function returning a trained Cox model
  if (useCoxLasso == FALSE) {
    cox.m <- NULL
    try(
      {
        cox.m <- survival::coxph(
          as.formula(
            paste(
              "survival::Surv(df_train$time, df_train$event) ~",
              paste(predict.factors, collapse = "+")
            )
          ),
          data = df_train, x = TRUE
        )
        # replace NAwith 0 i.e. ignore params that Cox couldn't estimate
        cox.m$coefficients[is.na(cox.m$coefficients)] <- 0
      },
      silent = TRUE
    )
    if (is.null(cox.m)) {
      print(paste(
        "Warning: cox.m == NULL, N/Events=",
        dim(df_train)[1],
        sum(df_train$event == 1)
      ))
    }
    return(cox.m)
  } else {
    return(survcoxlasso_train(df_train, predict.factors))
  }
}


#' Trains CoxLasso, using cv.glmnet(s="lambda.min")
#'
#' @param df_train  data frame with the data, "time" and "event" should describe survival outcome
#' @param predict.factors list of the column names to be used as predictors
#' @param fixed_time  not used here, to re-align with other methods
#' @param retrain_cox whether to re-train coxph on non-zero predictors; FALSE by default
#' @param verbose TRUE/FALSE prints warnings if no predictors in Lasso
#' @examples
#' df<- simsurv_nonlinear()
#' @return fitted CoxPH model
survcoxlasso_train <- function(df_train,
                                  predict.factors,
                                  fixed_time = NaN,
                                  retrain_cox = FALSE,
                                  verbose = FALSE) {
  stopifnot(expr = {
    is.data.frame(df_train)
    predict.factors %in% colnames(df_train)
  })

  cox.m <- NULL
  try(
    {
      cv10 <- glmnet::cv.glmnet(
        as.matrix(df_train[predict.factors]),
        survival::Surv(df_train$time, df_train$event),
        family = "cox",
        nfold = 5,
        alpha = 1
      )
      new.predictors <-
        rownames(coef(cv10, s = "lambda.min"))[as.matrix(coef(cv10, s = "lambda.min")) != 0]

      if (length(new.predictors) == 0) {
        if (verbose) {print("0 predictors in lasso!")}
        cox.m <-
          survival::coxph(
            survival::Surv(df_train$time, df_train$event) ~ 1,
            data = df_train,
            x = TRUE
          )
      } else {
        # check if non-regularised model to be re-trained
        if (retrain_cox) {
          # re-train cox on new.predictors
          f <-
            as.formula(paste(
              "survival::Surv(df_train$time, df_train$event) ~",
              paste(new.predictors, collapse = "+")
            ))
          cox.m <- survival::coxph(f, data = df_train, x = TRUE)
          # replace NA  with 0 i.e. ignore params that Cox couldn't estimate
          cox.m$coefficients[is.na(cox.m$coefficients)] <- 0
        } else {
          # return coxlasso in coxph object
          f <-
            as.formula(paste(
              "survival::Surv(df_train$time, df_train$event) ~",
              paste(predict.factors, collapse = "+")
            ))
          cox.m <- survival::coxph(f, data = df_train, x = TRUE)
          cox.m$coefficients <-
            as.numeric(coef(cv10, s = "lambda.min"))
        }
      }
    },
    silent = TRUE
  )
  return(cox.m)
}


#' Computes event probabilities from a trained cox model
#'
#' @param model_cox  pre-trained cox model
#' @param newdata data to compute event probabilities for
#' @param times  at which event probabilities are computed
#' @examples
#' df<- simsurv_nonlinear()
#' @return returns matrix(nrow = length(newdata), ncol = length(times))
survcox_predict <- function(model_cox,
                               newdata,
                               times) {
  # returns event probability from trained cox model model_cox

  # define bh - baseline hazard as dataframe with "time" and "hazard"
  # if baseline hazard can't be calibrated, # return mean(y) for all times
  # we take baseline hazard from K-M estimate and lp from Cox !!!! :((
  if (class(try(survival::basehaz(model_cox), silent = TRUE)) == "try-error") {
    bh <- summary(survival::survfit(model_cox$y ~ 1), times)$cumhaz
    predicted_event_prob <-
      matrix(nrow = dim(newdata)[1], ncol = length(times))
    for (i in seq(length(times))) {
      predicted_event_prob[, i] <- 1 -
        exp(-bh[i] * exp(
          predict(
            model_cox,
            newdata = newdata,
            type = "lp",
            reference = "zero"
          )
        ))
    }
    colnames(predicted_event_prob) <- round(times, 6)
    return(predicted_event_prob)
  } else {
    bh <- survival::basehaz(model_cox)
  }

  # define bh as function to compute bh for any time
  bh_approx <-
    stats::approxfun(bh[, "time"], bh[, "hazard"], method = "constant")

  # define bh_extrap how to extrapolate outside of the times in the training data
  if (class(try(stats::lm(hazard ~ poly(time, 3, raw = TRUE),
                   data = bh
  ), silent = TRUE)) != "try-error") {
    extrap <- stats::lm(hazard ~ poly(time, 3, raw = TRUE), data = bh)
    bh_extrap <- function(x) {
      sum(c(1, x, x**2, x**3) * extrap$coefficients[1:4])
    }
  } else {
    min_bh <- min(bh[, "hazard"], na.rm = 1)
    max_bh <- max(bh[, "hazard"], na.rm = 1)
    l <- dim(bh)[1]
    bh[1, c("hazard", "time")] <- c(0.0000001, min_bh)
    bh[l + 1, c("hazard", "time")] <-
      c(bh[l, "time"] + 100000, max_bh)
    bh_extrap <-
      stats::approxfun(bh[, "time"], bh[, "hazard"], method = "constant")
  }
  # compute event probability for times:
  # create placeholder
  predicted_event_prob <-
    matrix(nrow = dim(newdata)[1], ncol = length(times))
  # go over each time in times
  for (i in seq(length(times))) {
    if (is.na(bh_approx(times[i]))) {
      # if interpolation doesn't work, take extrapolated value
      bh_time <- bh_extrap(times[i])
      if (is.na(bh_time)) {
        bh_time <- mean(bh[, "hazard"], na.rm = TRUE)
      }
    } else {
      bh_time <- bh_approx(times[i])
    }
    # if baseline hazard is infinite, event probability is 1
    if (bh_time == Inf) {
      predicted_event_prob[, i] <- 1
      # if baseline hazard is ==0, event prob is 0 for all with survival==1
      # (somehow "survival" calculates even with bh==0)
    } else if (bh_time == 0) {
      predicted_event_prob[, i] <- 0
      # if baseline hazard is a number, use the survival formula
    } else {
      predicted_event_prob[, i] <-
        1 - exp(-bh_time * exp(
         predict(
            model_cox,
            newdata = newdata,
            type = "lp",
            reference = "zero"
          )
        ))
    }
  }
  # name columns by the time for which it predicts event prob
  colnames(predicted_event_prob) <- round(times, 6)
  return(predicted_event_prob)
}



#' Cross-validates Cox or CoxLasso model
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which event probabilities are computed
#' @param cv_number k in k-fold CV
#' @param randomseed random seed
#' @param useCoxLasso TRUE/FALSE
#' @param parallel TRUE/FALSE
#' @param return_models TRUE/FALSE, whether to return all CV models. Default is FALSE
#' @param repeat_cv if not NULL, repeats CV repeat_cv times
#' @examples
#' df<- simsurv_nonlinear()
#' @return output list: output$train, test, testaverage, traintaverage, time,tuned_cv_models
survcox_cv <- function(df,
                          predict.factors,
                          fixed_time = NaN,
                          cv_number = 5,
                          randomseed = NULL,
                          useCoxLasso = FALSE,
                          parallel = FALSE,
                          return_models = FALSE,
                          repeat_cv = NULL) {

  time_0 <- Sys.time()
  if (is.null(randomseed)) {randomseed <- round(stats::runif(1)*1e9,0)}

  stopifnot(expr = {
    is.data.frame(df)
    predict.factors %in% colnames(df)
  })

  if (sum(is.nan(fixed_time)) > 0) {
    fixed_time <- round(quantile(df[df$event == 1, "time"], 0.9), 1)
  }

  predict.factors <- eligible_params(predict.factors, df)
  if (length(predict.factors) == 0) {
    print("No eligible params")
    return(NULL)
  }

  #defining number of repeated cv
  if (is.null(repeat_cv)) {repeat_cv = 1}
  if (is.numeric(repeat_cv) & repeat_cv > 1) {
    repeat_cv = round(repeat_cv, 0)
  }else{
    repeat_cv = 1
  }

  print(paste("Cross-validating ",
              ifelse(useCoxLasso, "Cox-Lasso", "Cox-PH"),
              " ( ",repeat_cv,
              " repeat(s), ", cv_number," loops)",sep=""))

  modelstats_train <- list()
  modelstats_test <- list()
  models_for_each_cv <- list()
  #progress bar
  pb <- utils::txtProgressBar(0, cv_number*repeat_cv, style = 3)
  utils::setTxtProgressBar(pb, cv_number*repeat_cv / 50)

  for (rep_cv in 1:repeat_cv) {
    set.seed(randomseed + rep_cv)
    if (rep_cv != 1) {
      df <- df[sample(1:nrow(df)),]
    }
    cv_folds <-
      caret::createFolds(df$event, k = cv_number, list = FALSE)

    for (cv_iteration in 1:cv_number) {
      df_train_cv <- df[cv_folds != cv_iteration, ]
      df_test_cv <- df[cv_folds == cv_iteration, ]
      cox.model <- survcox_train(
        df_train_cv,
        predict.factors,
        useCoxLasso = useCoxLasso
      )
      y_predict_test <-
        survcox_predict(cox.model, df_test_cv, fixed_time)
      y_predict_train <-
        survcox_predict(cox.model, df_train_cv, fixed_time)
      modelstats_test[[cv_iteration + (rep_cv-1)*cv_number]] <-
        survval(y_predict_test, fixed_time, df_train_cv,
                df_test_cv, weighted = 1)
      modelstats_train[[cv_iteration + (rep_cv-1)*cv_number]] <-
        survval(y_predict_train,fixed_time,df_train_cv,
                df_train_cv,weighted = 1)
      if (return_models) {
        models_for_each_cv[[cv_iteration + (rep_cv-1)*cv_number]] <-
          cox.model
        }
      utils::setTxtProgressBar(pb, cv_iteration + (rep_cv-1)*cv_number)
      # create data frame with results:
    }
  }

  df_modelstats_test <- data.frame(modelstats_test[[1]])
  df_modelstats_train <- data.frame(modelstats_train[[1]])
  for (i in 2:(cv_number*repeat_cv)) {
    df_modelstats_test[i, ] <- modelstats_test[[i]]
    df_modelstats_train[i, ] <- modelstats_train[[i]]
  }
  row.names(df_modelstats_train)<- 1:(cv_number*repeat_cv)
  row.names(df_modelstats_test)<- 1:(cv_number*repeat_cv)

  utils::setTxtProgressBar(pb, cv_number*repeat_cv)
  close(pb)
  df_modelstats_test$test <- 1
  df_modelstats_train$test <- 0

  # comprise  output object
  output <- list()
  output$test <- df_modelstats_test
  output$train <- df_modelstats_train
  output$testaverage <-sapply(df_modelstats_test[, 1:8], mean, na.rm = 1)
  output$trainaverage <-sapply(df_modelstats_train[, 1:8], mean, na.rm = 1)
  output$tuned_cv_models <- models_for_each_cv
  time_1 <- Sys.time()
  output$time <- time_1 - time_0
  return(output)
}

#' Internal function, checks eligible parameters
#' to avoid tree-building with non-existing or constant predictors
#'
#' @param params list of predictor names
#' @param df data frame
#'
#' @return params_eligible - list of params which are in the df and have more than 1 value
#'
eligible_params <- function(params,df) {
  # This function checks eligible predictors from params list for split
  # It deletes those which are
  # 1) not in df and
  # 2) taking only 1 value (constants)
  # ! Later may delete collinear factors
  
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

