

#' Predicts event probability for a fitted survsrf_ens
#' @description
#' @param object trained survsrf_ens model
#' @param newdata test data
#' @param fixed_time  time for which probabilities are computed
#' @param ... other parameters to pass
#' @return matrix of predictions for observations in newdata by times
#' @export
survsrfens_predict <- function(object,
                                 newdata,
                                 fixed_time) {
  if (!inherits(object, "survensemble")){stop("Not a \"survensemble\" object")}
  if (!inherits(newdata, "data.frame")){stop("The data should be a data frame")}

  # use model_base with the base Cox model to find cox_predict
  newdata$cox_predict <- survcox_predict(object$model_base,newdata, fixed_time)
  predicted_event_prob <- survsrf_predict(object$model, newdata, fixed_time)
  return(predicted_event_prob)
}

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
#' @param tuningparams list of mtry, nodedepth and nodesize, to use default supply empty list()
#' @param useCoxLasso FALSE/TRUE, FALSE by default
#' @param var_importance_calc FALSE/TRUE, TRUE by default
#' @return trained object of class survsrf_ens
#' @export
survsrfens_train <- function(df_train,
                               predict.factors,
                               fixed_time = NaN,
                               inner_cv = 3,
                               randomseed = NULL,
                               tuningparams = list(),
                               useCoxLasso = FALSE,
                               max_grid_size =10,
                               var_importance_calc = FALSE) {
  # the function trains Cox model, then adds its predictions
  # into Survival Random Forest model
  # to mimic stacking procedure and reduce overfitting,
  # we train Cox model on 0.9 of the data and predict
  # on the rest 0.1 for each 1/10s fold
  # so we pass out-of-the-bag prediction to SRF

  if (length(eligible_params(predict.factors, df_train)) == 0) {
    print("No eliible params")
    return(NULL)
  }
  if (is.nan(fixed_time)| (length(fixed_time) > 1)) {
    fixed_time <-
      round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
  }
  if (is.nan(randomseed)) {
    randomseed <- round(stats::runif(1) * 1e9, 0)
  }
  set.seed(randomseed)
  #out-of-sample Cox predictions
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
    # predict for cox_oob - linear predictors, not event probabilities
    cox_predict_oob <- predict(cox_m_cv, cox_oob, type = "lp") #beta x X

    # adding Cox prediction to the df_train in the column "cox_predict"
    df_train[cv_folds == cv_iteration, "cox_predict"] <- cox_predict_oob
  }
  predict.factors.plusCox <- c(predict.factors, "cox_predict")

  srf.ens <-
    survsrf_train(
      df_train = df_train,
      predict.factors = predict.factors.plusCox,
      fixed_time = fixed_time,
      max_grid_size = max_grid_size,
      inner_cv = inner_cv,
      randomseed = randomseed,
      tuningparams = tuningparams
    )
  if (var_importance_calc) {
    v <- randomForestSRC::vimp(
      srf.ens$model,importance = "permute",seed = randomseed)
    var_importance <- sort(v$importance, decreasing = TRUE)
    vimp10 <- var_importance[1:min(length(var_importance), 20)]
  } else { vimp10 <- c(NaN)  }

  #base cox model
  cox_base_model <-
    survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

  #output
  output = list()
  output$model_name = "SRF_ensemble"
  output$model <- srf.ens$model
  output$model_base <- cox_base_model
  output$randomseed <- randomseed
  output$bestparams <- srf.ens$bestparams
  output$call <-  match.call()
  output$vimp10 <- vimp10
  class(output) <- "survensemble"
  return(output)
}

#' Cross-validates predictive performance for SRF Ensemble
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param outer_cv number of folds in outer CV, default 3
#' @param inner_cv number of folds for model tuning CV, default 3
#' @param repeat_cv number of CV repeats, if NULL, runs once
#' @param randomseed random seed
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param useCoxLasso TRUE/FALSE, default is FALSE
#' @param tuningparams list of tuning parameters for random forest: 1) NULL for using a default tuning grid, or 2) a list("mtry"=c(...), "nodedepth" = c(...), "nodesize" = c(...))
#' @examples \donttest{
#' \dontshow{rfcores_old <- options()$rf.cores; options(rf.cores=1)}
#' df <- simulate_nonlinear()
#' ens_cv <- survsrf_ens_cv(df, names(df)[1:4])
#' summary(ens_cv)
#' \dontshow{options(rf.cores=rfcores_old)}
#' }
#' @return list of outputs
#' @export
survsrfens_cv <- function(df,
                            predict.factors,
                            fixed_time = NaN,
                            outer_cv = 3,
                            inner_cv = 3,
                            repeat_cv = 2,
                            randomseed = NULL,
                            return_models = FALSE,
                            useCoxLasso = FALSE,
                            tuningparams = list(),
                            max_grid_size = 10,
                            parallel = FALSE) {
  Call <- match.call()
  inputs <- list(df , predict.factors, fixed_time,
                 outer_cv,inner_cv, repeat_cv,
                 randomseed, return_models,
                 useCoxLasso,tuningparams, parallel)
  inputclass<- list(df = "data.frame", predict.factors = "character",
                    fixed_time = "numeric",outer_cv = "numeric",
                    inner_cv = "numeric", repeat_cv = "numeric",
                    randomseed = "numeric",return_models = "logical",
                    useCoxLasso="logical", tuningparams = "list",
                   parallel = "logical")

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
    train_function = survsrfens_train,
    predict_function = survsrfens_predict,
    model_args = list(
      "useCoxLasso" = useCoxLasso,
      "tuningparams" = tuningparams,
      "fixed_time" = fixed_time,
      "max_grid_size" = max_grid_size,
      "randomseed" = randomseed
    ),
    model_name = "SRF_ensemble",
    parallel= parallel
  )
  output$call <- Call
  return(output)
}

