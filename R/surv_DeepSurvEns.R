
################  ens_deepsurv_train ####################
# Create out-of-bag Cox predictions, then train deepsurv
#' @export
ens_deepsurv_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NaN,
           useCoxLasso = FALSE,
           tuningparams = list(),
           max_grid_size = 10) {

    if (is.nan(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }

    # setting fixed_time if not given
    if (is.nan(fixed_time)| (length(fixed_time) > 1)) {
      fixed_time <-
        round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
    }

    #creating folds
    set.seed(randomseed)
    cv_folds <-
      caret::createFolds(df_train$event, k = 5, list = FALSE)
    cindex_train <- vector(length = 5)
    cindex_test <- vector(length = 5)
    for (cv_iteration in 1:5) {
      cox_train <- df_train[cv_folds != cv_iteration, ]
      cox_oob <- df_train[cv_folds == cv_iteration, ]
      # train cox model on cox_train
      cox_m_cv <-
        survcox_train(cox_train,
                      eligible_params(predict.factors, cox_train),
                      useCoxLasso = useCoxLasso)
      # predict for cox_oob
      cox_predict_oob <-
        predict(cox_m_cv, cox_oob, type = "lp") #(beta x X)
        #survcox_predict(cox_m_cv, cox_oob, fixed_time)
      # adding Cox prediction to the df_train in the column "cox_predict"
      df_train[cv_folds == cv_iteration, "cox_predict"] <- cox_predict_oob
    }

    # adding Cox predictions as a new factor to tune SRF,
    predict.factors.plusCox <- c(predict.factors, "cox_predict")

    # train the deepsurv model
    deepsurv.ens <-
      deepsurv_train(df_train = df_train,
                    predict.factors = predict.factors.plusCox,
                    fixed_time = fixed_time,
                    tuningparams = tuningparams,
                    max_grid_size = max_grid_size,
                    inner_cv = inner_cv,
                    randomseed = randomseed )
    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)
    #output
    output = list()
    output$model_name = "DeepSurv_ensemble"
    output$model <- deepsurv.ens$model
    output$model_base <- cox_base_model
    output$randomseed <- randomseed
    output$bestparams <- deepsurv.ens$bestparams
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }

################  ens_deepsurv_predict ####################

#same as deepsurv_predict
#' @export
ens_deepsurv_predict <-
  function(trained_object,
           newdata,
           fixed_time,
           predict.factors
  ) {
    trained_model <- trained_object$model
    predictdata <- newdata[predict.factors]
    # use model_base with the base Cox model to find cox_predict
    predictdata$cox_predict <- survcox_predict(trained_object$model_base,
                                               newdata, fixed_time)
    s1 <- predict(trained_model, newdata = predictdata, type = "survival")
    ## if this failed, take Cox predictions
    if (is.na(s1[1,1])) {
      predict_eventprob<- predictdata$cox_predict
    }else{
      f <- function(i) {
        approxfun(as.double(colnames(s1)), s1[i,], method = "linear")(fixed_time)
      }
      predict_eventprob <- 1 - unlist(lapply(1:dim(s1)[1], f))
    }
    return(predict_eventprob)
  }

############### ens_deepsurv_CV #############
#' @export
ens_deepsurv_cv <- function(df,
                           predict.factors,
                           fixed_time = NaN,
                           outer_cv = 3,
                           inner_cv = 3,
                           repeat_cv = 2,
                           randomseed = NaN,
                           return_models = FALSE,
                           useCoxLasso = FALSE,
                           tuningparams = list(),
                           max_grid_size = 10,
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
    train_function = ens_deepsurv_train,
    predict_function = ens_deepsurv_predict,
    model_args = list("tuningparams" = tuningparams,
                      "useCoxLasso" = useCoxLasso,
                      "max_grid_size" = max_grid_size,
                      "randomseed" = randomseed,
                      "fixed_time" = fixed_time),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "DeepSurv_ensemble",
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
