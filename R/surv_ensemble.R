######################### Ensemble 1 ########################

#' Fits an ensemble of Cox-PH and Survival Random Forest (SRF)
#' with internal CV to tune SRF hyperparameters.
#'
#' Details: the function trains Cox model, then adds its out-of-the-box
#' predictions to Survival Random Forest as an additional predictor
#' to mimic stacking procedure used in Machine Learning and reduce over-fitting.
#' #' Cox model is fitted to .9 data to predict the rest .1 for each 1/10s fold;
#' these out-of-the-bag predictions are passed on to SRF
#'
#' @param df_train  data, "time" and "event"  describe survival outcome
#' @param predict.factors list of the column names to be used as predictors
#' @param fixed_time  for which the performance is maximized
#' @param inner_cv  number of inner cycles for model tuning
#' @param randomseed  random seed
#' @param srf_tuning list of mtry, nodedepth and nodesize, to use default supply empty list()
#' @param fast_version  TRUE/FALSE, TRUE by default
#' @param oob FALSE/TRUE, TRUE by default
#' @param useCoxLasso FALSE/TRUE, FALSE by default
#' @param var_importance_calc FALSE/TRUE, TRUE by default
#' @return trained object of class survensemble
#' @export
survensemble_train <- function(df_train,
                               predict.factors,
                               fixed_time = NaN,
                               inner_cv = 3,
                               randomseed = NULL,
                               srf_tuning = list(),
                               fast_version = TRUE,
                               oob = TRUE,
                               useCoxLasso = FALSE,
                               var_importance_calc = 1) {
  # the function trains Cox model, then adds its predictions
  # into Survival Random Forest model
  # to mimic stacking procedure and reduce overfitting,
  # we train Cox model on 0.9 of the data and predict
  # on the rest 0.1 for each 1/10s fold
  # so we pass out-of-the-bag prediction to SRF

  Call <- match.call()

  predict.factors <- eligible_params(predict.factors, df_train)
  if (length(predict.factors) == 0) {
    print("No eliible params")
    return(NULL)
  }

  # defining output for fixed_time
  if (sum(is.nan(fixed_time)) > 0) {
    fixed_time <-
      round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
  }
  #setting random seed
  if (is.null(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0)
  }
  set.seed(randomseed)
  #creating folds
  cv_folds <-
    caret::createFolds(df_train$event, k = 10, list = FALSE)
  cindex_train <- vector(length = 10)
  cindex_test <- vector(length = 10)
  for (cv_iteration in 1:10) {
    cox_train <- df_train[cv_folds != cv_iteration, ]
    cox_oob <- df_train[cv_folds == cv_iteration, ]
    # train cox model on cox_train
    cox_m_cv <-
      survcox_train(cox_train,
                    eligible_params(predict.factors, cox_train),
                    useCoxLasso = useCoxLasso)
    # predict for cox_oob
    cox_predict_oob <-
      survcox_predict(cox_m_cv, cox_oob, fixed_time)
    # adding Cox prediction to the df_train in the column "cox_predict"
    df_train[cv_folds == cv_iteration, "cox_predict"] <-
      cox_predict_oob
  }

  # adding Cox predictions as a new factor to tune SRF
  predict.factors.1A <- c(predict.factors, "cox_predict")
  ensemble1_model <-
    survsrf_train(
      df_train = df_train,
      predict.factors = predict.factors.1A,
      fixed_time = fixed_time,
      inner_cv = inner_cv,
      randomseed = randomseed,
      srf_tuning = srf_tuning,
      fast_version = fast_version,
      oob = oob
    )

  if (var_importance_calc) {
    v <- randomForestSRC::vimp(ensemble1_model$model,
                               importance = "permute",
                               seed = randomseed)
    var_importance <- sort(v$importance, decreasing = TRUE)
    vimp10 <- var_importance[1:min(length(var_importance), 20)]
  } else {
    vimp10 <- c(NaN)
  }

  #base cox model
  cox_base_model <-
    survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

  #output
  ensemble1_model$vimp10 <- vimp10
  ensemble1_model$model_base <- cox_base_model
  ensemble1_model$randomseed <- randomseed
  ensemble1_model$call <- Call
  class(ensemble1_model) <- "survensemble"
  return(ensemble1_model)
}


#' Predicts event probability for a fitted survensemble
#' @description
#' \link[survcompare:predict.survensemble]{predict.survensemble}
#' @param object trained survensemble model
#' @param newdata test data
#' @param fixed_time  time for which probabilities are computed
#' @param oob TRUE/FALSE , default is FALSE, if out of bag predictions are to be made from SRF
#' @param ... other parameters to pass
#' @return matrix of predictions for observations in newdata by times
#' @export
predict.survensemble <- function(object,
                                 newdata,
                                 fixed_time,
                                 oob = FALSE,
                                 ...) {
  if (!inherits(object, "survensemble")) {
    stop("Not a \"survensemble\" object")
  }
  if (is.null(newdata)) {
    stop("The data for predictions is not supplied")
  }
  if (!inherits(newdata, "data.frame")) {
    stop("The data should be a data frame")
  }

  # use model_base with the base Cox model to find cox_predict
  newdata$cox_predict <- survcox_predict(object$model_base,
                                         newdata, fixed_time)
  # now use "model" which is SRF which needs additional risk factor
  # "cox_predict" which was created in the previous row
  predicted_event_prob <-
    1 - srf_survival_prob_for_time(object$model, newdata, fixed_time, oob = oob)
  return(predicted_event_prob)
}


#' Cross-validates predictive performance for Ensemble 1
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param outer_cv k in k-fold CV, default 3
#' @param inner_cv kk in the inner look of kk-fold CV, default 3
#' @param repeat_cv if NULL, runs once (or 1), otherwise repeats CV
#' @param randomseed random seed
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param useCoxLasso TRUE/FALSE, default is FALSE
#' @param srf_tuning list of tuning parameters for random forest: 1) NULL for using a default tuning grid, or 2) a list("mtry"=c(...), "nodedepth" = c(...), "nodesize" = c(...))
#' @param oob TRUE/FALSE use out-of-bag predictions while tuning instead of cross-validation, TRUE by default
#' @examples \donttest{
#' \dontshow{rfcores_old <- options()$rf.cores; options(rf.cores=1)}
#' df <- simulate_nonlinear()
#' ens_cv <- survensemble_cv(df, names(df)[1:4])
#' summary(ens_cv)
#' \dontshow{options(rf.cores=rfcores_old)}
#' }
#' @return list of outputs
#' @export
survensemble_cv <- function(df,
                            predict.factors,
                            fixed_time = NaN,
                            outer_cv = 3,
                            inner_cv = 3,
                            repeat_cv = 2,
                            randomseed = NULL,
                            return_models = FALSE,
                            useCoxLasso = FALSE,
                            srf_tuning = list(),
                            oob = TRUE) {
  Call <- match.call()
  inputs <- list(df , predict.factors, fixed_time,
                 outer_cv,inner_cv, repeat_cv,
                 randomseed, return_models,
                 useCoxLasso,srf_tuning, oob)
  inputclass<- list(df = "data.frame", predict.factors = "character", fixed_time = "numeric",
                    outer_cv = "numeric",inner_cv = "numeric", repeat_cv = "numeric",
                    randomseed = "numeric",return_models = "logical",
                    useCoxLasso="logical", srf_tuning = "list", oob = "logical")

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
    train_function = survensemble_train,
    predict_function = predict.survensemble,
    model_args = list(
      "useCoxLasso" = useCoxLasso,
      "srf_tuning" = srf_tuning,
      "oob" = oob
    ),
    model_name = "CoxPH and SRF Ensemble"
  )
  output$call <- Call
  return(output)
}

##################################################################

#' Prints trained survensemble object
#'
#'@param x survensemble object
#'@param ... additional arguments to be passed
#'@return x
#'@export
print.survensemble <- function(x, ...) {
  if (!inherits(x, "survensemble")) {
    stop("Not a \"survensemble\" object")
  }
  summary.survensemble(x)
}

#' Prints summary of a trained survensemble object
#'
#'@param object survensemble object
#'@param ... additional arguments to be passed
#'@return object
#'@export
summary.survensemble <- function(object, ...) {
  if (!inherits(object, "survensemble")) {
    stop("Not a \"survensemble\" object")
  }
  cat("Survival ensemble of Cox PH and Survival Random Forest.\n")
  if (!is.null(cl <- object$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\n=> Ensemble Cox-SRF:\n")
  print(object$model)
  cat("\n=> Underlying CoxPH:\n")
  print(object$model_base)
  cat("\n=> Items available as object$item are: ")
  cat(names(object), sep = ", ")

}

##################################################################
#' Prints survensemble_cv object
#'
#'@param x survensemble_cv object
#'@param ... additional arguments to be passed
#'@return x
#'@export
print.survensemble_cv <- function(x, ...) {
  if (!inherits(x, "survensemble_cv")) {
    stop("\nNot a \"survensemble_cv\" object")
  }
  summary.survensemble_cv(x)
}


#' Prints a summary of survensemble_cv object
#'
#'@param object survensemble_cv object
#'@param ... additional arguments to be passed
#'@return object
#'@export
summary.survensemble_cv <- function(object, ...) {
  if (!inherits(object, "survensemble_cv")) {
    stop("Not a \"survensemble_cv\" object")
  }
  cat("Cross-validation results\n")

  if (!is.null(cl <- object$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nThe stats are computed from the ",
      dim(object$test)[1],
      " data splits.\n")
  print(object$summarydf)
}
