
################  ens_deepsurv_train ####################
# Create out-of-bag Cox predictions, then train DeepSurv
#' @export
ens_deepsurv_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NULL,
           useCoxLasso = FALSE,
           deepsurvparams = list()) {
    
    #setting random seed
    if (is.null(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    set.seed(randomseed)
    
    # setting fixed_time if not given
    if (sum(is.nan(fixed_time)) > 0 | (length(fixed_time) > 1)) {
      fixed_time <-
        round(quantile(df_train[df_train$event == 1, "time"], 0.9), 1)
    }
    
    #creating folds
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
        survcox_predict(cox_m_cv, cox_oob, fixed_time)
      # adding Cox prediction to the df_train in the column "cox_predict"
      df_train[cv_folds == cv_iteration, "cox_predict"] <- cox_predict_oob
    }
    
    # cox_m<- survcox_train(df_train, predict.factors,
    #                                  useCoxLasso = useCoxLasso)
    #df_train$cox_predict <- survcox_predict(cox_m,df_train,fixed_time)
    
    # adding Cox predictions as a new factor to tune SRF,
    predict.factors.plusCox <- c(predict.factors[-1], "cox_predict")
    
    # train the DeepSurv model
    deepsurv.ens <-
      deepsurv_train(df_train,predict.factors.plusCox,deepsurvparams)
    
    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)
    
    #output
    output = list()
    output$model <- deepsurv.ens
    output$model_base <- cox_base_model
    output$randomseed <- randomseed
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
    trained_model<- trained_object$model
    predictdata <- newdata[predict.factors[-1]]
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
ens_deepsurv_CV <- function(df,
                            predict.factors,
                            fixed_time = NaN,
                            outer_cv = 3,
                            inner_cv = 3,
                            repeat_cv = 2,
                            randomseed = NULL,
                            return_models = FALSE,
                            useCoxLasso = FALSE,
                            deepsurvparams = list()
) {
  Call <- match.call()
  # inputs <- list(df , predict.factors, fixed_time,
  #                outer_cv,inner_cv, repeat_cv,
  #                randomseed, return_models,
  #                useCoxLasso,deepsurvparams)
  # inputclass<- list(df = "data.frame", predict.factors = "character", fixed_time = "numeric",
  #                   outer_cv = "numeric",inner_cv = "numeric", repeat_cv = "numeric",
  #                   randomseed = "numeric",return_models = "logical",
  #                   useCoxLasso="logical", deepsurvparams = "list")
  # cp<- check_call(inputs, inputclass, Call)
  # if (cp$anyerror) stop (paste(cp$msg[cp$msg!=""], sep=""))
  
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
    model_args = list("deepsurvparams" = deepsurvparams),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "CoxPH and DeepSurv Ensemble"
  )
  output$call <- Call
  return(output)
}