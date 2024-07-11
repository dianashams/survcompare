

################## deepsurv_predict ##################
# The function to predict event probability by a trained deepsurvmodel
#' @export
deepsurv_predict <-
  function(trained_model,
           newdata,
           fixed_time,
           predict.factors
  ){
    s1 <-
      predict(trained_model, newdata = newdata[predict.factors], type = "survival")
    #times<- as.double(colnames(s1))
    f <- function(i) {
      approxfun(as.double(colnames(s1)), s1[i,], method = "linear")(fixed_time)
    }
    predict_eventprob <- 1 - unlist(lapply(1:dim(s1)[1], f))
    return(predict_eventprob)
  }

################## deepsurv_train ##################
# The function trains DeepSurv model with given params
#' @export
deepsurv_train <-
  function(df_train,
           predict.factors,
           deepsurvparams = list()
  ) {
    if (length(deepsurvparams) ==0) {
      deepsurvm <- survivalmodels::deepsurv(
        data = df_train,
        x = df_train[predict.factors],
        y = Surv(df_train$time, df_train$event))
    }else{
      if(is.null(deepsurvparams$dropout)) {deepsurvparams$dropout= 0.2}
      if(is.null(deepsurvparams$learning_rate)) {deepsurvparams$learning_rate = 0.01}
      if(is.null(deepsurvparams$num_nodes)) {deepsurvparams$num_nodes =  c(16, 16)}
      if(is.null(deepsurvparams$batch_size)) {deepsurvparams$batch_size =
                                              round(max(100, dim(df_train)[1] / 4), 0) }
      if(is.null(deepsurvparams$epochs)) {deepsurvparams$epochs = 50}
      deepsurvm <- survivalmodels::deepsurv(
        data = df_train,
        x = df_train[predict.factors],
        y = Surv(df_train$time, df_train$event),
        dropout = deepsurvparams$dropout,
        learning_rate = deepsurvparams$learning_rate,
        num_nodes = deepsurvparams$num_nodes,
        batch_size = deepsurvparams$batch_size,
        epochs = deepsurvparams$epochs
      )
    }
    return(deepsurvm)
  }

################## deepsurv_tune ##################
#' This function repeats 3-fold CV repeat-tune times
#' and finds grid with highest c-index
#' @export
deepsurv_tune <-
  function(repeat_tune,
           df_tune,
           predict.factors,
           fixed_time,
           inner_cv = 3) {

    #fixed_time
    if (is.nan(fixed_time) | length(fixed_time) > 1) {
      # not implemented for multiple time
      fixed_time <- round(quantile(df_tune[df_tune$event == 1, "time"], 0.9), 2)
    }

    means = c()
    for (i in (1:repeat_tune)) {
      print (repeat_tune)
      ds_tune_fus <-
        deepsurv_tune_single(df_tune, predict.factors,fixed_time,inner_cv)
      means <- cbind(means, ds_tune_fus$cindex_mean)
    }
    grid = ds_tune_fus$grid
    allmeans <- apply(means, 1, mean)
    # grid and average cindex ordered by cindex (highest first)
    output = list()
    output$cindex_ordered <-
      cbind(grid, allmeans)[order(allmeans, decreasing = TRUE),]
    output$bestparams <- output$cindex_ordered[1,]
    return(output)
  }

#' The function trains DeepSurv model with given params
#' @export
deepsurv_tune_single <-
  function(df_tune,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3) {

    #fixed_time
    if (is.nan(fixed_time) | length(fixed_time) > 1) {
      # not implemented for multiple time
      fixed_time <- round(quantile(df_tune[df_tune$event == 1, "time"], 0.9), 2)
    }

    # defining the tuning grid
    grid_of_values <- expand.grid(
      "dropout" = c(0.1, 0.3),
      "learning_rate" = c(0.001, 0.01),
      "num_nodes" =  list(c(16, 16), c(32, 32), c(16, 16, 16), c(32, 32, 32))
    ) #2*2*4=16
    grid_size <- dim(grid_of_values)[1]
    #placeholder for c-index
    cind = matrix(NA, nrow = grid_size, ncol = inner_cv)
    # tuning cross-validation loop

    #progress bar
    pb <- utils::txtProgressBar(0, inner_cv * grid_size, style = 3)
    utils::setTxtProgressBar(pb, inner_cv * grid_size / 50)

    for (cv_iteration in 1:inner_cv) {
      # print(paste("DeepSurv tuning CV step", cv_iteration, "of", inner_cv))
      cv_folds <-
        caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
      df_train_cv <- df_tune[cv_folds != cv_iteration, ]
      df_test_cv <- df_tune[cv_folds == cv_iteration, ]
      #Grid search
      for (i in 1:grid_size) {
        deepsurvm <- deepsurv(
          data = df_train,
          x = df_train_cv[predict.factors],
          y = Surv(df_train_cv$time, df_train_cv$event),
          shuffle = TRUE,
          epochs = 50,
          #early_stopping = TRUE,
          frac = 0.2,
          batch_size = round(max(100, dim(df_train_cv)[1] / 4), 0),
          dropout = grid_of_values[i, "dropout"],
          learning_rate = grid_of_values[i, "learning_rate"],
          num_nodes = grid_of_values[i, "num_nodes"][[1]]
        )
        #check test performance
        pp <-
          deepsurv_predict(trained_model = deepsurvm,newdata =  df_test_cv, fixed_time =fixed_time, predict.factors = predict.factors)
        cind[i, cv_iteration] =
          surv_validate(pp, fixed_time, df_train_cv, df_test_cv)[1, "C_score"]
        #cind[i, cv_iteration] =
        #  concordance(Surv(df_test_cv$time, df_test_cv$event) ~ pp)$concordance
        utils::setTxtProgressBar(pb, grid_size + (i - 1) * cv_iteration)
      }
      # here we already have cindex for allgrid for cv_iteration
    }
    utils::setTxtProgressBar(pb, inner_cv * grid_size)
    close(pb)
    remove(deepsurvm)
    output = list()
    output$grid = grid_of_values
    output$cindex = cind
    output$cindex_mean = apply(cind, 1, mean)
    return(output)
  }

################## deepsurv_CV ########################

# The function trains DeepSurv model with given params
#' @export
deepsurv_cv <- function(df,
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
    train_function = deepsurv_train,
    predict_function = deepsurv_predict,
    model_args = list("deepsurvparams" = deepsurvparams),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "DeepSurv"
  )
  output$call <- Call
  return(output)
}
