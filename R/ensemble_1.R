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
#' @param srf_tuning list of mtry, nodedepth and nodesize, default is NULL
#' @param fast_version  TRUE/FALSE, TRUE by default
#' @param oob FALSE/TRUE, TRUE by default
#' @param useCoxLasso FALSE/TRUE, FALSE by default
#' @param var_importance_calc FALSE/TRUE, TRUE by default
#' @examples
#' df <- simsurv_nonlinear()
#' @return trained object of class survensemble1
#' @export
survensemble1 <- function(df_train,
                          predict.factors,
                          fixed_time = NaN,
                          inner_cv = 3,
                          randomseed = NULL,
                          srf_tuning = NULL,
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
  if (is.null(randomseed)){randomseed <- round(stats::runif(1) * 1e9, 0)}
  set.seed(randomseed)
  #creating folds
  cv_folds <-
    caret::createFolds(df_train$event, k = 10, list = FALSE)
  cindex_train <- vector(length = 10)
  cindex_test <- vector(length = 10)
  for (cv_iteration in 1:10) {
    cox_train <- df_train[cv_folds != cv_iteration,]
    cox_oob <- df_train[cv_folds == cv_iteration,]
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
    survrf_train(
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
    survcox_train(df_train,predict.factors,useCoxLasso = useCoxLasso)

  #output
  ensemble1_model$vimp10 <- vimp10
  ensemble1_model$model_base <- cox_base_model
  ensemble1_model$randomseed <- randomseed
  ensemble1_model$call <- Call
  class(ensemble1_model)<- "survensemble1"
  return(ensemble1_model)
}



#' Predicts event probability for a fitted Ensemble 1 method
#' @description
#' \link[survcompare:predict.survensemble1]{predict.survensemble1}
#'
#' @param survens1 trained survensemble1 model
#' @param df_test test data
#' @param fixed_time  time for which probabilities are computed
#' @param oob TRUE/FALSE , default is FALSE
#' @param useCoxLasso TRUE/FALSE , default is FALSE
#' @examples
#' df <- simsurv_nonlinear()
#' @return returns matrix(nrow = length(newdata), ncol = length(times))
#' @export
predict.survensemble1 <- function(survens1,
                                  df_test,
                                  fixed_time,
                                  oob = FALSE,
                                  useCoxLasso = FALSE) {
  if (!inherits(survens1, "survensemble1")) {stop("Not a \"survensemble1\" object")}
  if (is.null(df_test)) {stop("The data for predictions is not supplied")}

  # use model_base with the base Cox model to find cox_predict
  df_test$cox_predict <- survcox_predict(survens1$model_base,
                                         df_test, fixed_time)
  # now use "model" which is SRF which needs additional risk factor
  # "cox_predict" which was created in the previous row
  predicted_event_prob <-
    1 - srf_survival_prob_for_time(survens1$model, df_test, fixed_time, oob = oob)
  return(predicted_event_prob)
}


#' Cross-validates predictive performance for Ensemble 1
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param cv_number k in k-fold CV, default 3
#' @param inner_cv kk in the inner look of kk-fold CV, default 3
#' @param randomseed random seed
#' @param useCoxLasso TRUE/FALSE, default is FALSE
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param srf_tuning list of mtry, nodedepth, nodesize to tune, default is NULL
#' @param repeat_cv if NULL, runs once, otherwise repeats CV
#' @examples
#' df <- simsurv_nonlinear()
#' @return output list: output$train, test, testaverage, traintaverage, time
#' @export
survensemble1cv <- function(df,
                           predict.factors,
                           fixed_time = NaN,
                           cv_number = 3,
                           inner_cv = 3,
                           randomseed = NULL,
                           useCoxLasso = FALSE,
                           return_models = FALSE,
                           srf_tuning = NULL,
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
    print("No eliible params")
    return(NULL)
  }

  Call <- match.call()

  #defining number of repeated cv
  if (is.null(repeat_cv)) {repeat_cv = 1}
  if (is.numeric(repeat_cv) & repeat_cv > 1) {
    repeat_cv = round(repeat_cv, 0)
  }else{
    repeat_cv = 1
  }

  print(paste(
    "Cross-validating Survival Random Forest - Cox model ensemble ( ",
    repeat_cv," repeat(s), ", cv_number," outer, ",
    inner_cv," inner loops)",sep="")
    )

  modelstats_train <- list()
  modelstats_test <- list()
  models_for_each_cv <- list()
  #progress bar
  pb <- utils::txtProgressBar(0, cv_number*repeat_cv, style = 3)
  utils::setTxtProgressBar(pb, cv_number*repeat_cv / 50)

  for (rep_cv in 1:repeat_cv) {
    set.seed(randomseed + rep_cv)
    if (rep_cv!=1) {df <- df[sample(1:nrow(df)), ]}
    cv_folds <-
      caret::createFolds(df$event, k = cv_number, list = FALSE)
    for (cv_iteration in 1:cv_number) {
      #print(paste("CV ", cv_iteration,"/",cv_number, sep=""))
      df_train_cv <- df[cv_folds != cv_iteration,]
      df_test_cv <- df[cv_folds == cv_iteration,]
      model.tuned <- survensemble1(
        df_train = df_train_cv,
        predict.factors = predict.factors,
        fixed_time = fixed_time,
        inner_cv = inner_cv,
        randomseed = randomseed,
        fast_version = TRUE,
        oob = TRUE,
        useCoxLasso = useCoxLasso,
        srf_tuning = srf_tuning
      )
      #  calculating tuned model predictions
      y_predict_test <-
        predict.survensemble1(model.tuned, df_test_cv, fixed_time, oob = FALSE)
      y_predict_train <-
        predict.survensemble1(model.tuned, df_train_cv, fixed_time, oob = FALSE)
      modelstats_test[[cv_iteration + (rep_cv-1)*cv_number]] <-
        survval(y_predict_test, fixed_time, df_train_cv,
                        df_test_cv, weighted = 1)
      modelstats_train[[cv_iteration + (rep_cv-1)*cv_number]] <-
        survval(y_predict_train,fixed_time,df_train_cv,
                        df_train_cv,weighted = 1)
      if (return_models) {
        models_for_each_cv[[cv_iteration + (rep_cv-1)*cv_number]] <-
          model.tuned
      }
      utils::setTxtProgressBar(pb, cv_iteration + (rep_cv-1)*cv_number)
    }
  }
  df_modelstats_test<- data.frame(modelstats_test[[1]])
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

  output <- list()
  output$test <- df_modelstats_test
  output$train <- df_modelstats_train
  output$testaverage <- sapply(df_modelstats_test, mean, na.rm = 1)
  output$trainaverage <-  sapply(df_modelstats_train, mean, na.rm = 1)
  output$tuned_cv_models <- models_for_each_cv
  output$randomseed <- randomseed
  output$call <- Call
  class(output) <- "survensemble1cv"
  time_1 <- Sys.time()
  output$time <- time_1 - time_0
  print(time_1 - time_0)
  return(output)
}

#' prints survensemble1 object
#' @param obj1 survensemble1 object
#' @export
print.survensemble1<- function(obj1){
  if (!inherits(obj1, "survensemble1")) {stop("Not a \"survensemble1\" object")}
  summary.survensemble1(obj1)
}

#' summary for a trained Cox-SRF ensemble, should be a survensemble1 object
#' @param obj1 survensemble1 object
#' @export
summary.survensemble1<- function(obj1){
  if (!inherits(obj1, "survensemble1")) {stop("Not a \"survensemble1\" object")}
  cat("Survival ensemble of Cox PH and Survival Random Forest.\n")
  if(!is.null(cl<- obj1$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\n=> Ensemble Cox-SRF:\n")
  print(obj1$model)
  cat("\n=> Underlying CoxPH:\n")
  print(obj1$model_base)
  cat("\n=> Items available as object$item are: ")
  cat(names(obj1), sep=", ")
}
