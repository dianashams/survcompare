
#' Computes event probabilities from a trained ensemble of CoxPH and DeepHit
#' @param trained_object a trained model, output of survdhens_train()
#' @param newdata new data for which predictions are made
#' @param fixed_time time of interest for which event probabilities are computed
#' @param predict.factors list of predictor names
#' @param extrapsurvival if probabilities are extrapolated beyond trained times (constant)
#' @return vector of predicted event probabilities
#' @export
survdhens_predict <-
  function(trained_object,
           newdata,
           fixed_time,
           predict.factors,
           extrapsurvival = TRUE) {

    if (!inherits(trained_object, "survensemble")) {
      stop("Not a \"survensemble\" object")
    }
    if (!inherits(newdata, "data.frame")) {
      stop("The data should be a data frame")
    }
    # use model_base with the base Cox model to find cox_predict
    newdata$cox_predict <- survcox_predict(trained_object$model_base,
                                           newdata = newdata,fixed_time = fixed_time)
    # use deephit_predict()
    predict_eventprob <-
      survdeephit_predict(trained_object$model,newdata, fixed_time,
                      predict.factors = c(predict.factors, "cox_predict"),
                      extrapsurvival = extrapsurvival
      )
    return(predict_eventprob)
  }


#' Trains ensemble of the CoxPH and DeepHit
#' @param df_train  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param randomseed random seed to control tuning including data splits
#' @param useCoxLasso if CoxLasso is used (TRUE) or not (FALSE, default)
#' @param tuningparams if given, list of hyperparameters, list(mod_alpha=c(), ...), otherwise a wide default grid is used
#' @param max_grid_size number of random grid searches for model tuning
#' @param verbose FALSE(default)/TRUE
#' @export
survdhens_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NaN,
           useCoxLasso = FALSE,
           tuningparams = list(),
           max_grid_size =10,
           verbose = verbose) {

    Call <- match.call()

    if (is.nan(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    if (length(eligible_params(predict.factors, df_train)) == 0) {
      print("No eliible params")
      return(NULL)
    }

    # setting fixed_time if not given
    if (is.nan(fixed_time)| (length(fixed_time) > 1)) {
      fixed_time <-
        round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
    }

    #out-of-sample Cox predictions
    set.seed(randomseed)
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

    # adding Cox predictions as a new factor to tune SRF,
    predict.factors.plusCox <- c(predict.factors, "cox_predict")

    # train the deephit model
    deephit.ens <-
      survdeephit_train(df_train = df_train,
                    predict.factors = predict.factors.plusCox,
                    fixed_time = fixed_time,
                    tuningparams = tuningparams,
                    max_grid_size = max_grid_size,
                    inner_cv = inner_cv,
                    randomseed = randomseed,
                    verbose = verbose)

    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

    #output
    output = list()
    output$model_name = "DeepHit_ensemble"
    output$model <- deephit.ens$model
    output$model_base <- cox_base_model
    output$randomseed <- randomseed
    output$bestparams <- deephit.ens$bestparams
    output$grid <- deephit.ens$grid_of_hyperparams
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }

#' Cross-validates the performance of DeepHit and CoxPH Ensemble
#' @param df  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param outer_cv number of cross-validation folds for model validation
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param repeat_cv number of CV repeats, if NaN, runs once
#' @param randomseed random seed to control tuning including data splits
#' @param return_models TRUE/FALSE, if TRUE returns all trained models
#' @param useCoxLasso if CoxLasso is used (TRUE) or not (FALSE, default)
#' @param tuningparams if given, list of hyperparameters, list(mod_alpha=c(), ...), otherwise a wide default grid is used
#' @param max_grid_size number of random grid searches for model tuning
#' @param parallel if parallel calculations are used
#' @param verbose FALSE(default)/TRUE
#' @param package_path survcompare package path if not installed as a library
#' @param python_path python path for survivalmodels
#' @export
survdhens_cv <- function(df,
                         predict.factors,
                         fixed_time = NaN,
                         outer_cv = 3,
                         inner_cv = 3,
                         repeat_cv = 2,
                         randomseed = NaN,
                         return_models = FALSE,
                         useCoxLasso = FALSE,
                         tuningparams = list(),
                         max_grid_size =10,
                         parallel = FALSE,
                         verbose = FALSE,
                         package_path = NaN,
                         python_path = NaN
) {
  Call <- match.call()

  if (sum(is.na(df[c("time", "event", predict.factors)])) > 0) {
    stop("Missing data can not be handled. Please impute first.")
  }
  if (parallel&(is.nan(package_path)|(is.nan(python_path)))){
    stop("Please supply package and python paths for parallel computations.")
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
    train_function = survdhens_train,
    predict_function = survdhens_predict,
    model_args = list("tuningparams" = tuningparams,
                      "useCoxLasso" = useCoxLasso,
                      "max_grid_size" = max_grid_size,
                      "randomseed" = randomseed,
                      "fixed_time" = fixed_time,
                      "verbose"= verbose),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "DeepHit_ensemble",
    parallel = parallel,
    package_path = package_path,
    python_path = python_path
  )
  output$call <- Call
  return(output)
}
