#' Predicts event probability by a trained stacked ensemble of Survival Random Forest and CoxPH
#'
#' @param trained_object a trained model, output of survsrfstack_train()
#' @param newdata new data for which predictions are made
#' @param fixed_time time of interest, for which event probabilities are computed
#' @param predict.factors list of predictor names
#' @param extrapsurvival if probabilities are extrapolated beyond trained times (constant)
#' @return vector of predicted event probabilities
#' @export
survsrfstack_predict <-
  function(trained_object,
           newdata,
           fixed_time,
           predict.factors,
           extrapsurvival=TRUE
  ) {
    predictdata <- newdata[predict.factors]
    l <- trained_object$lambda

    # use model_base with the base Cox model to find cox_predict
    predictdata$cox_predict <-
      survcox_predict(trained_model = trained_object$model_base_cox,
                      newdata = newdata, fixed_time = fixed_time)
    # if it is just Cox model, i.e. lambda = 0, return Cox predictions
    if ((!use_alternative_model) & (l==0)) {return (predictdata$cox_predict)}
    # otherwise compute ML predictions
    predictdata$ml_predict <-
      survsrf_predict(trained_model = trained_object$model_base_ml,
                      newdata = newdata,
                      fixed_time = fixed_time,
                      extrapsurvival = extrapsurvival)
    #weighted sum
    p <- predictdata$cox_predict +
      l * (predictdata$ml_predict- predictdata$cox_predict)
    return(p)
  }

# Trains stacked ensemble of the CoxPH and Survival Random Forest
#' @param df_train  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param randomseed random seed to control tuning including data splits
#' @param useCoxLasso if CoxLasso is used (TRUE) or not (FALSE, default)
#' @param tuningparams if given, list of hyperparameters, list(mtry=c(), nodedepth=c(),nodesize=c()), otherwise a wide default grid is used
#' @param max_grid_size number of random grid searches for model tuning
#' @param verbose FALSE(default)/TRUE
#' @return output = list(bestparams, allstats, model)
#' @examples
#' d <-simulate_nonlinear(100),
#' p<- names(d)[1:4]
#' tuningparams = list(
#'  "mtry" = c(5,10,15),
#'  "nodedepth" = c(5,10,15,20),
#'  "nodesize" =    c(20,30,50)
#' )
#' m_srf<- survsrf_train(d,p,tuningparams=tuningparams)
#' @export
survsrfstack_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NaN,
           useCoxLasso = FALSE,
           tuningparams = list(),
           max_grid_size =10,
           verbose = FALSE) {

    if (is.nan(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    # setting fixed_time if not given
    if (is.nan(fixed_time)| (length(fixed_time) > 1)) {
      fixed_time <-
        round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
    }
    # if fixed_time is outside of the event times, we return an error -
    #  extrapolation is not supproted
    if (fixed_time > max(df_train$time)) {
      print ("The calibration time given is outside of the train data.")
      return (NULL)
    }

    if(verbose) cat("\nTraining baseline learners: SRF...\n")
    #base srf
    ml_base_model <-
      survsrf_train(df_train = df_train,
                    predict.factors = predict.factors,
                    fixed_time = fixed_time,
                    max_grid_size = max_grid_size,
                    tuningparams = tuningparams,
                    inner_cv = inner_cv,
                    randomseed = randomseed,
                    verbose = verbose)
    if (verbose) cat("CoxPH...")
    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

    #first, using k-fold CV, underlying models are trained
    # (using already tuned parameters)
    #and out-of-sample predictions for Cox and ML are computed
    #then, a regression (outcome ~ cox_predictions + ml_predictions)
    #constitutes the meta-learner

    #meta-learner: compute out-of-sample Cox and ML (srf) predictions
    if(verbose) cat("\nComputing out-of-sample predictions...\t")

    k_for_oob = 5
    goodsplit <- function(max_tries = 10) {
      # if maximum event time is earlier than fixed_time, SRF won't predict
      # we try 10 different splits to get this right, otherwise show NaNs
      for (split_trial in 1:max_tries) {
        {set.seed(randomseed + split_trial-1)
          cv_folds <- caret::createFolds(df_train$event, k = k_for_oob, list = FALSE)}
        max.eventtimes =
          unlist(lapply(
            1:k_for_oob,
            FUN = function(i) max(df_train[df_train$event == 1 & cv_folds != split_trial, "time"])
          ))
        # return a good split straight away
        if (sum(max.eventtimes < fixed_time) == 0) {return(cv_folds)}
      }
      return(cv_folds) # if not, return the last one
    }
    cv_folds <- goodsplit(10)
    cindex_train <- vector(length = k_for_oob)
    cindex_test <- vector(length = k_for_oob)

    for (cv_iteration in 1:k_for_oob) {
      if(verbose) cat("\t", cv_iteration, "/", k_for_oob)
      data_cvtrain <- df_train[cv_folds != cv_iteration, ]
      data_oob <- df_train[cv_folds == cv_iteration, ]

      # train cox model on data_cvtrain
      cox_m_cv <-
        survcox_train(data_cvtrain,predict.factors,useCoxLasso = useCoxLasso)
      # predict for data_oob
      cox_predict_oob <-
        survcox_predict(cox_m_cv, data_oob, fixed_time)
      # adding Cox prediction to the df_train in the column "cox_predict"
      df_train[cv_folds == cv_iteration, "cox_predict"] <- cox_predict_oob
      srf.cv <-
          survsrf_train(df_train = data_cvtrain,
                      predict.factors = predict.factors,
                      fixed_time = fixed_time,
                      tuningparams = as.list(ml_base_model$bestparams),
                      max_grid_size = 1,
                      inner_cv = inner_cv,
                      randomseed = randomseed + cv_iteration,
                      verbose = FALSE)

      # predict for data_oob
      ml_predict_oob <-
        survsrf_predict(trained_model = srf.cv,
                        newdata = data_oob,
                        fixed_time= fixed_time)
      # adding ML prediction to the df_train in the column "ml_predict"
      df_train[cv_folds == cv_iteration, "ml_predict"] <- ml_predict_oob
      remove(cox_m_cv)
      remove(srf.cv)
    }

    if(verbose) cat("\n Calibrating meta-learner ...")
    # find lambda that gives highest c-score using oob cox and ml predictions
    c_score <- c()
    lambdas <- seq(0,1,0.01)
    for (i in 1:length(lambdas)){
      y_hat <- with(df_train, cox_predict + lambdas[i]*(ml_predict - cox_predict))
      temp <-try(survival::concordancefit(
        survival::Surv(df_train$time, df_train$event),-1 * y_hat),
        silent = TRUE)
      c_score[i] <-
        ifelse((inherits(temp, "try-error")) |
                 is.null(temp$concordance), NaN, temp$concordance)
    }
    best_i = ifelse(sum(!is.nan(c_score))==0, 1, which.max(c_score))
    worst_i = ifelse(sum(!is.nan(c_score))== 0, 1, which.min(c_score))
    bestparams_meta <-
      c("lambda" = lambdas[best_i],
        "c_score" = c_score[best_i],
        "lambda_worst" = lambdas[worst_i],
        "c_score_worst" = c_score[worst_i]
      )

    if(verbose) cat("\t Lambda = ", lambdas[best_i],", in Cox + lambda * (ML - Cox).")

    # alternative stacked model (unrestricted to lambda in 0-1)
    # 1/0 by fixed_time:
    df_train$event_t <-
      ifelse(df_train$time <= fixed_time & df_train$event == 1, 1, 0)
    logit = function(x) {
      log(pmax(pmin(x, 0.9999), 0.0001) / (1 - pmax(pmin(x, 0.9999), 0.0001)))
    }
    df_train$cox_predict_logit <- logit(df_train$cox_predict)
    df_train$ml_predict_logit <- logit(df_train$ml_predict)
    # Excluding censored observations before fixed_time, leave those with known state
    df_train_in_scope <-
      df_train[(df_train$time >= fixed_time) |
                 (df_train$time < fixed_time & df_train$event == 1),]
    model_meta_alternative <-
      glm(event_t ~ cox_predict_logit + ml_predict_logit,
          data = df_train_in_scope,
          family = "binomial")

    #output
    output = list()
    output$model_name <- "Stacked_SRF_CoxPH"
    output$oob_predictions <- df_train
    output$lambda <- lambdas[best_i]
    output$model_meta_alternative <- model_meta_alternative
    output$model_base_cox <- cox_base_model
    output$model_base_ml <- ml_base_model
    output$randomseed <- randomseed
    output$bestparams <- c( ml_base_model$bestparams, bestparams_meta)
    output$grid = ml_base_model$grid_of_hyperparams
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }

#' Cross-validates stacked ensemble of the CoxPH and Survival Random Forest models
#' @param df  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param outer_cv number of cross-validation folds for model validation
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param repeat_cv number of CV repeats, if NaN, runs once
#' @param randomseed random seed to control tuning including data splits
#' @param useCoxLasso if CoxLasso is used (TRUE) or not (FALSE, default)
#' @param tuningparams if given, list of hyperparameters, list(mtry=c(), nodedepth=c(),nodesize=c()), otherwise a wide default grid is used
#' @param max_grid_size number of random grid searches for model tuning
#' @param parallel if parallel calculations are used
#' @param verbose FALSE(default)/TRUE
#' @export
survsrfstack_cv <- function(df,
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
                         verbose = FALSE
) {
  Call <- match.call()

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
    train_function = survsrfstack_train,
    predict_function = survsrfstack_predict,
    model_args = list("tuningparams" = tuningparams,
                      "useCoxLasso" = useCoxLasso,
                      "max_grid_size" = max_grid_size,
                      "randomseed" = randomseed,
                      "fixed_time" = fixed_time,
                      "verbose" = verbose),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "Stacked_SRF_CoxPH",
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
