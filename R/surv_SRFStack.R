
################  stack_srf_train ####################
# Create out-of-bag Cox predictions, then train srf
#' @export
stack_srf_train <-
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
    if (!is.nan(randomseed)) {set.seed(randomseed)}
    if(verbose) cat("\nTraining baseline learners ...\n")
    #base srf
    ml_base_model <-
      survsrf_train(df_train = df_train,
                    predict.factors = predict.factors,
                    fixed_time = fixed_time,
                    tuningparams = tuningparams,
                    #max_grid_size = max_grid_size,
                    inner_cv = inner_cv,
                    randomseed = randomseed )

    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)

    #------ Computationally demanding fit of the stack model:---
    #first, using k-fold CV, underlying models are trained (using already tuned parameters though)
    #and out-of-sample predictions for Cox and ML are computed
    #then, a regression (outcome ~ cox_predictions + ml_predictions)
    #constitutes the meta-learner
    #meta-learner: compute out-of-sample Cox and ML (srf) predictions
    if(verbose) cat("\nTraining meta-learner... computing out-of-sample predictions...\t")

    tuningparams_tuned = list(mtry = ml_base_model$bestparams$mtry,
                              nodesize = ml_base_model$bestparams$nodesize,
                              nodedepth = ml_base_model$bestparams$nodedepth)
    bestparams_base <- ml_base_model$bestparams

    #avoid variance in predictions for small data, use higher number of folds(not activated)
    k_for_oob = ifelse(dim(df_train)[1]<=250, 5, 5)

    cv_folds <-
      caret::createFolds(df_train$event, k = k_for_oob, list = FALSE)
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
      srf.cv <-
        survsrf_train(df_train = data_train,
                      predict.factors = predict.factors,
                      fixed_time = fixed_time,
                      tuningparams = tuningparams_tuned,
                      #max_grid_size = 1,
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
    output$bestparams <- c(bestparams_base, bestparams_meta)
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }

################  stack_srf_predict ####################

#same as srf_predict
#' @export
stack_srf_predict <-
  function(trained_object,
           newdata,
           fixed_time,
           predict.factors,
           use_alternative_model = FALSE
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
                      fixed_time = fixed_time)
    #weighted sum
    p <- predictdata$cox_predict +
      l * (predictdata$ml_predict- predictdata$cox_predict)

    # alternative_model = logistic regression of Cox and ML predictions,
    if(use_alternative_model) {
      m <- trained_object$model_meta_alternative
      logit = function(x) {
        log(pmax(pmin(x, 0.9999), 0.0001) / (1 - pmax(pmin(x, 0.9999), 0.0001)))
      }
      predictdata$cox_predict_logit <- logit(predictdata$cox_predict)
      predictdata$ml_predict_logit <- logit(predictdata$ml_predict)
      return (predict(m, newdata = predictdata))
    }
    # return weighted sum
    return(p)
  }

############### stack_srf_cv #############
#' @export
stack_srf_cv <- function(df,
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
                         parallel = FALSE
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
    train_function = stack_srf_train,
    predict_function = stack_srf_predict,
    model_args = list("tuningparams" = tuningparams,
                      "useCoxLasso" = useCoxLasso,
                      "max_grid_size" = max_grid_size,
                      "randomseed" = randomseed,
                      "fixed_time" = fixed_time),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "Stacked_SRF_CoxPH",
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
