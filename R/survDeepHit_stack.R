#' survdhstack_predict Computes event probabilities from a trained stacked ensemble of DeepHit and CoxPH
#' @param trained_object a trained model, output of survdhstack_train()
#' @param newdata new data for which predictions are made
#' @param fixed_time time of interest for which event probabilities are computed
#' @param predict.factors list of predictor names
#' @param extrapsurvival if probabilities are extrapolated beyond trained times (constant)
#' @return vector of predicted event probabilities
#' @export
survdhstack_predict <-
  function(trained_object,
           newdata,
           fixed_time,
           predict.factors,
           extrapsurvival = FALSE
  ) {
    predictdata <- newdata[predict.factors]
    l <- trained_object$lambda

    # use model_base with the base Cox model to find cox_predict
    predictdata$cox_predict <-
      survcox_predict(trained_model = trained_object$model_base_cox,
                      newdata = newdata, fixed_time = fixed_time)
    # if it is just Cox model, i.e. lambda = 0, return Cox predictions
    if (l==0) {return (predictdata$cox_predict)}
    # otherwise compute ML predictions
    predictdata$ml_predict <-
      survdeephit_predict(trained_model = trained_object$model_base_ml,
                      newdata = newdata, predict.factors = predict.factors,
                      fixed_time = fixed_time,extrapsurvival = extrapsurvival)
    #weighted sum
    p <- predictdata$cox_predict +
      l * (predictdata$ml_predict- predictdata$cox_predict)
    # return weighted sum
    return(p)
  }

#' Trains stacked ensemble of the CoxPH and DeepHit
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
survdhstack_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NaN,
           useCoxLasso = FALSE,
           tuningparams = list(),
           max_grid_size =10,
           verbose = FALSE) {

    # setting fixed_time if not given
    if (is.nan(fixed_time)| (length(fixed_time) > 1)) {
      fixed_time <-
        round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
    }
    if (is.nan(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    if(verbose) cat("\nTraining baseline learners ...\n")
    #base DeepHit
    ml_base_model <-
      survdeephit_train(df_train = df_train,
                    predict.factors = predict.factors,
                    fixed_time = fixed_time,
                    tuningparams = tuningparams,
                    max_grid_size = max_grid_size,
                    inner_cv = inner_cv,
                    randomseed = randomseed )

    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

    #------ Stack model training:---
    # first, using k-fold CV, underlying models are trained (using already tuned parameters though)
    # and out-of-sample predictions for Cox and ML are computed
    # then, a regression (outcome ~ cox_predictions + ml_predictions)
    # constitutes the meta-learner
    # meta-learner: compute out-of-sample Cox and ML (DeepHit) predictions
    # find lambda with the highest c-score using obtained cox and ml predictions
    if(verbose) cat("\nTraining meta-learner... computing out-of-sample predictions...\t")

    #avoid variance in predictions for small data, can use higher number of folds(not activated)
    k_for_oob = ifelse(dim(df_train)[1]<=250, 5, 5)

    goodsplit <- function(max_tries = 10) {
      # if maximum event time is earlier than fixed_time, DH won't predict
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
      data_train <- df_train[cv_folds != cv_iteration, ]
      data_oob <- df_train[cv_folds == cv_iteration, ]
      # train cox model on data_train
      cox_m_cv <- survcox_train(data_train,predict.factors,useCoxLasso = useCoxLasso)
      # predict for data_oob
      cox_predict_oob <-
        survcox_predict(cox_m_cv, data_oob, fixed_time)
      # adding Cox prediction to the df_train in the column "cox_predict"
      df_train[cv_folds == cv_iteration, "cox_predict"] <- cox_predict_oob
      #train ML model on data_train
      deephit.cv <-
        survdeephit_train(df_train = data_train,
                      predict.factors = predict.factors,
                      fixed_time = fixed_time,
                      tuningparams = as.list(ml_base_model$bestparams),
                      max_grid_size = 1,
                      inner_cv = inner_cv,
                      randomseed = randomseed + cv_iteration,
                      verbose = verbose)
      # predict for data_oob
      ml_predict_oob <-
        survdeephit_predict(trained_model = deephit.cv,
                        newdata = data_oob,
                        predict.factors = predict.factors,
                        fixed_time= fixed_time)
      # adding ML prediction to the df_train in the column "ml_predict"
      df_train[cv_folds == cv_iteration, "ml_predict"] <- ml_predict_oob
    }

    if(verbose) cat("\t calibrating meta-learner ...")
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

    #output
    output = list()
    output$model_name <- "Stacked_DeepHit_CoxPH"
    output$oob_predictions <- df_train
    output$lambda <- lambdas[best_i]
    output$model_base_cox <- cox_base_model
    output$model_base_ml <- ml_base_model
    output$randomseed <- randomseed
    output$bestparams <- c(ml_base_model$bestparams, bestparams_meta)
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }


#' Cross-validates Stacked Ensemble of CoxPH and DeepHit
#' @param df  data, "time" and "event" should describe survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time time at which performance is maximized
#' @param outer_cv number of cross-validation folds for model validation
#' @param inner_cv number of cross-validation folds for hyperparameters' tuning
#' @param repeat_cv number of cross-validation repetitions
#' @param randomseed random seed to control tuning including data splits
#' @param return_models TRUE/FALSE, if TRUE returns all trained models
#' @param useCoxLasso if CoxLasso is used (TRUE) or not (FALSE, default)
#' @param tuningparams if given, list of hyperparameters, list(mtry=c(), nodedepth=c(),nodesize=c()), otherwise a wide default grid is used
#' e.g. (mod_alpha = c(0.2,0.5),dropout= c(0.1,0.8), learning_rate= c(0.001,0.01), num_nodes = list(c(64,64)), epochs = 10, sigma = c(0.1,1,10), cuts =50,batch_size=50, early_stopping = FALSE, weight_decay= 0)
#' @param max_grid_size number of random grid searches for model tuning
#' @param parallel if parallel calculations are used
#' @param verbose FALSE(default)/TRUE
#' @param package_path survcompare package path if not installed as a library
#' @param python_path python path for survivalmodels
#' @param suppresswarn TRUE/FALSE, TRUE by default
#' @export
survdhstack_cv <- function(df,
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
                           python_path = NaN,
                           suppresswarn = TRUE
) {
  Call <- match.call()
  if (sum(is.na(df[c("time", "event", predict.factors)])) > 0) {
    stop("Missing data can not be handled. Please impute first.")
  }
  if (suppresswarn){ user_warn <-options()$warn; options(warn=-1)}
  output <- surv_CV(
    df = df,
    predict.factors = predict.factors,
    fixed_time = fixed_time,
    outer_cv = outer_cv,
    inner_cv = inner_cv,
    repeat_cv = repeat_cv,
    randomseed = randomseed,
    return_models = return_models,
    train_function = survdhstack_train,
    predict_function = survdhstack_predict,
    model_args = list("tuningparams" = tuningparams,
                      "useCoxLasso" = useCoxLasso,
                      "max_grid_size" = max_grid_size,
                      "randomseed" = randomseed,
                      "fixed_time" = fixed_time,
                      "verbose"= verbose),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "Stacked_DeepHit_CoxPH",
    parallel = parallel,
    package_path = NaN,
    python_path = NaN
  )
  if (suppresswarn){ options(warn=user_warn)}
  output$call <- Call
  return(output)
}
