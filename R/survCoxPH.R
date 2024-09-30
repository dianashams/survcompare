###################### Basic Cox Model functions #################################
#' Trains CoxPH using survival package, or trains CoxLasso (cv.glmnet, lambda.min),
#'  and then re-trains survival:coxph on non-zero predictors
#'
#' @param df_train  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of the column names to be used as predictors
#' @param fixed_time  target time, NaN by default; needed here only to re-align with other methods
#' @param useCoxLasso  TRUE or FALSE
#' @param retrain_cox if useCoxLasso is TRUE, whether to re-train coxph on non-zero predictors, FALSE by default
#' @param inner_cv k in k-fold CV for training lambda for Cox Lasso, only used for useCoxLasso = TRUE
#' @return fitted CoxPH or CoxLasso model
#' @export
survcox_train <- function(df_train,
                          predict.factors,
                          fixed_time = NaN,
                          useCoxLasso = FALSE,
                          retrain_cox = FALSE,
                          inner_cv = 5) {
  Call <- match.call()
  indx <-
    pmatch(c("df_train", "predict.factors"), names(Call), nomatch = 0)
  if (indx[1] * indx[2] == 0) {
    stop("Please supply data and predictors")
  }
  stopifnot(
    "The data is not a data frame" = inherits(df_train, "data.frame"),
    "Predictors are not found" = inherits(predict.factors, "character"),
    "Predictors are not in the data supplied" = predict.factors %in% colnames(df_train)
  )
  # if Lasso, then return survcoxlasso_train()
  if (useCoxLasso){
    return(survcoxlasso_train(df_train = df_train,
                              predict.factors = predict.factors,
                              inner_cv = inner_cv,
                              retrain_cox = retrain_cox))
  }
  # if not Lasso
  # wrapper for coxph() function returning a trained Cox model
  cox.m <- NULL
  try({
    cox.m <- survival::coxph(as.formula(
      paste(
        "survival::Surv(df_train$time, df_train$event) ~",
        paste(predict.factors, collapse = "+")
      )
    ),
    data = df_train, x = TRUE)
    # replace NA with 0 i.e. ignore params that Cox couldn't estimate
    cox.m$coefficients[is.na(cox.m$coefficients)] <- 0
  },
  silent = TRUE)
  if (is.null(cox.m)) {
    print(paste("Warning: cox.m == NULL, N/Events=",
                dim(df_train)[1],sum(df_train$event == 1)))
  }
  return(cox.m)
}


#' Trains CoxLasso, using cv.glmnet(s="lambda.min")
#'
#' @param df_train  data frame with the data, "time" and "event" should describe survival outcome
#' @param predict.factors list of the column names to be used as predictors
#' @param inner_cv k in k-fold CV for lambda tuning
#' @param fixed_time  not used here, for internal use
#' @param retrain_cox whether to re-train coxph on non-zero predictors; FALSE by default
#' @param verbose TRUE/FALSE prints warnings if no predictors in Lasso
#' @return fitted CoxPH object with coefficient of CoxLasso or re-trained CoxPH with non-zero CoxLasso if retrain_cox = FALSE or TRUE
#' @export
survcoxlasso_train <- function(df_train,
                               predict.factors,
                               inner_cv = 5,
                               fixed_time = NaN,
                               retrain_cox = FALSE,
                               verbose = FALSE) {
  stopifnot(expr = {
    is.data.frame(df_train)
    predict.factors %in% colnames(df_train)
  })

  cox.m <- NULL
  try({
    cv10 <- glmnet::cv.glmnet(
      as.matrix(df_train[predict.factors]),
      survival::Surv(df_train$time, df_train$event),
      family = "cox",
      nfold = inner_cv,
      alpha = 1,
      maxit = 200
    )
    new.predictors <-
      rownames(coef(cv10, s = "lambda.min"))[as.matrix(coef(cv10, s = "lambda.min")) != 0]

    if (length(new.predictors) == 0) {
      if (verbose) {print("Warning: No predictors are left in lasso.")}
      cox.m <-
        survival::coxph(
          survival::Surv(df_train$time, df_train$event) ~ 1,
          data = df_train,
          x = TRUE
        )
      return(cox.m)
    }

    # check if non-regularized model to be re-trained
    if (retrain_cox) {
      # re-train cox on new.predictors
      f <-as.formula(paste("survival::Surv(df_train$time, df_train$event) ~",
                           paste(new.predictors, collapse = "+")))
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
      temp<- as.numeric(coef(cv10, s = "lambda.min"))
      names(temp) = predict.factors
      cox.m$coefficients <- temp
    }
  },
  silent = TRUE)
  return(cox.m)
}

#' Computes event probabilities from a trained cox model
#'
#' @param trained_model  pre-trained cox model of coxph class
#' @param newdata data to compute event probabilities for
#' @param fixed_time  at which event probabilities are computed
#' @param interpolation "constant" by default, can also be "linear", for between times interpolation for hazard rates
#' @return returns matrix(nrow = length(newdata), ncol = length(fixed_time))
#' @export
survcox_predict <- function(trained_model,
                            newdata,
                            fixed_time,
                            interpolation = "constant") {
  # returns event probability from trained cox model trained_model

  #checks
  if (!inherits(trained_model, "coxph")) {
    stop("Supply coxph model.");    return(NULL)
  }
  if (!inherits(newdata, "data.frame"))   stop("Supply newdata as data.frame.")
  if (!inherits(fixed_time, "numeric"))
    stop("Supply fixed_time as a non-empty numeric list.")
  if (length(fixed_time) == 0)
    stop("Supply fixed_time as a non-empty numeric list.")
  if (is.null(newdata) |  dim(newdata)[1] == 0 |  dim(newdata)[2] == 0)
    stop("Empty or NULL data is supplied.")

  # compute baseline hazard function
  temp <- try(survival::basehaz(trained_model), silent = TRUE)

  #compute exponentiated linear predictors
  explp <- predict(trained_model, newdata, type = "risk") #exp(beta x X)

  # if baseline hazard can't be calibrated, return NaNs
  if (inherits(temp, "try-error")) {return(rep(NaN, dim(newdata)[1]))   }
  remove(temp)

  # define bh as function to compute bh for variable time points
  bh <- survival::basehaz(trained_model)
  bh_approx <-
    stats::approxfun(bh[, "time"], bh[, "hazard"], method = interpolation)

  # compute event probability if possible
  # this would not work for the times outside of training data, return NaN
  if (is.na(bh_approx(fixed_time))) {
    return(rep(NaN, dim(newdata)[1]))
    } else {
      if(fixed_time > max(bh[,"time"])) {return(rep(NaN, dim(newdata)[1]))}
      bh_t <- bh_approx(fixed_time)
    }
  # if baseline cumulative hazard == Inf, event probability is 1
  # if  == 0, event prob is 0 for all with survival==1
  # (somehow "survival" calculates even with bh==0)
  if (bh_t == Inf) {
      predicted_event_prob <- rep(1, dim(newdata)[1])
    } else if (bh_t == 0) {
      predicted_event_prob <- rep(0, dim(newdata)[1])
    } else {
      predicted_event_prob <- 1 - exp(-bh_t * explp)
    }

  return(predicted_event_prob)
}



#' Cross-validates Cox or CoxLasso model
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param outer_cv k in k-fold CV, default 3
#' @param repeat_cv if NULL, runs once, otherwise repeats CV
#' @param randomseed random seed
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param inner_cv k in the inner loop of k-fold CV, default is 3; only used if CoxLasso is TRUE
#' @param useCoxLasso TRUE/FALSE, FALSE by default
#' @examples \donttest{
#' df <- simulate_nonlinear()
#' coxph_cv <- survcox_cv(df, names(df)[1:4])
#' summary(coxph_cv)
#' }
#' @return list of outputs
#' @export
survcox_cv <- function(df,
                       predict.factors,
                       fixed_time = NaN,
                       outer_cv = 3,
                       repeat_cv = 2,
                       randomseed = NaN,
                       return_models = FALSE,
                       inner_cv = 3,
                       useCoxLasso = FALSE,
                       parallel = FALSE) {
  Call <- match.call()
  inputs <- list(df , predict.factors, fixed_time,
                 outer_cv,inner_cv, repeat_cv,
                 randomseed, return_models,
                 useCoxLasso)
  inputclass<- list(df = "data.frame", predict.factors = "character", fixed_time = "numeric",
                    outer_cv = "numeric",inner_cv = "numeric", repeat_cv = "numeric",
                    randomseed = "numeric",return_models = "logical",
                    useCoxLasso = "logical")
  cp<- check_call(inputs, inputclass, Call)
  if (cp$anyerror) stop (paste(cp$msg[cp$msg!=""], sep=""))

  if (sum(is.na(df[c("time", "event", predict.factors)])) > 0) {
    stop("Missing data can not be handled. Please impute first.")
  }

  output <- surv_CV(
    df = df,
    predict.factors = predict.factors,
    fixed_time = fixed_time,
    outer_cv = outer_cv,
    inner_cv = inner_cv,
    repeat_cv = repeat_cv,
    randomseed = randomseed,
    return_models = return_models,
    train_function = survcox_train,
    predict_function = survcox_predict,
    model_args = list("useCoxLasso" = useCoxLasso),
    model_name = ifelse(!useCoxLasso, "CoxPH", "CoxLasso"),
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
