

################## deephit_predict ##################
# The function to predict event probability by a trained deephitmodel
#' @export
deephit_predict <-
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

################## deephit_train ##################
# The function trains deephit model with given params
#' @export
deephit_train <-
  function(df_train,
           predict.factors,
           deephitparams = list()
  ) {

    if (length(deephitparams) ==0) {
      deephitm <- deephit(
        data = df_train,
        x = df_train[predict.factors],
        y = Surv(df_train$time, df_train$event))
    }else{
      if(is.null(deephitparams$dropout)) {deephitparams$dropout= 0.2}
      if(is.null(deephitparams$learning_rate)) {deephitparams$learning_rate = 0.01}
      if(is.null(deephitparams$num_nodes)) {deephitparams$num_nodes =  c(16, 16,16,2)}
      if(is.null(deephitparams$batch_size)) {deephitparams$batch_size =  round(max(100, dim(df_train)[1] / 4), 0) }
      if(is.null(deephitparams$epochs)) {deephitparams$epochs =  30}
      deephitm <- deephit(
        data = df_train,
        x = df_train[predict.factors],
        y = Surv(df_train$time, df_train$event),
        dropout = deephitparams$dropout,
        learning_rate = deephitparams$learning_rate,
        num_nodes = deephitparams$num_nodes,
        batch_size = deephitparams$batch_size,
        epochs = deephitparams$epochs
      )
    }
    return(deephitm)
  }

################## deephit_tune ##################
# This function repeats 3-fold CV repeat-tune times
# and finds grid with highest c-index
#' @export
deephit_tune <-
  function(repeat_tune,
           df_tune,
           predict.factors,
           fixed_time,
           inner_cv = 3) {
    means = c()
    for (i in (1:repeat_tune)) {
      print (repeat_tune)
      ds_tune_fus <-
        deephit_tune_single(df_tune, predict.factors,fixed_time,inner_cv)
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

################## deepsurv_tune ##################
#' 3-fold CV tuning
#' @export
deephit_tune_single <-
  function(df_tune,
           predict.factors,
           fixed_time,
           inner_cv = 3) {
    # defining the tuning grid
    grid_of_values <- expand.grid(
      "dropout" = c(0.1,0.3),
      "epochs" = c(5,50),
      "learning_rate" =c(0.001),
      "frac" = c(0.1, 0.3),
      "num_nodes" =    list(c(8,8), c(16,16,16,2), c(32,32,32,4))
    ) #20
    grid_size <- dim(grid_of_values)[1]
    #placeholder for c-index
    cind = matrix(NA, nrow = grid_size, ncol = inner_cv)

    #progress bar
    pb <- utils::txtProgressBar(0, inner_cv * grid_size, style = 3)
    utils::setTxtProgressBar(pb, inner_cv * grid_size / 50)

    # tuning cross-validation loop
    for (cv_iteration in 1:inner_cv) {
      # print(paste("deephit tuning CV step", cv_iteration, "of", inner_cv))
      cv_folds <-
        caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
      df_train_cv <- df_tune[cv_folds != cv_iteration, ]
      df_test_cv <- df_tune[cv_folds == cv_iteration, ]
      #Grid search
      for (i in 1:grid_size) {
        deephitm <- survivalmodels::deephit(
          data = df_train,
          x = df_train_cv[predict.factors],
          y = Surv(df_train_cv$time, df_train_cv$event),
          shuffle = TRUE,
          early_stopping = TRUE,
          dropout = grid_of_values[i, "dropout"],
          epochs = grid_of_values[i, "epochs"],
          frac = grid_of_values[i, "frac"],
          learning_rate = grid_of_values[i, "learning_rate"],
          num_nodes = grid_of_values[i, "num_nodes"][[1]],
          batch_size = round(max(100, dim(df_train_cv)[1] / 4), 0)
        )
        #check test performance
        pp <-
          deephit_predict(deephitm, df_test_cv, predict.factors, fixed_time)
        cind[i, cv_iteration] =
          surv_validate(pp, fixed_time, df_train_cv, df_test_cv)[1, "C_score"]
        #cind[i, cv_iteration] =
        #  concordance(Surv(df_test_cv$time, df_test_cv$event) ~ pp)$concordance
        utils::setTxtProgressBar(pb, grid_size + (i - 1) * cv_iteration)
      } # here we already have cindex for allgrid for cv_iteration
    }
    utils::setTxtProgressBar(pb, inner_cv * grid_size)
    close(pb)
    remove(deephitm)
    output = list()
    output$grid = grid_of_values
    output$cindex = cind
    output$cindex_mean = apply(cind, 1, mean)
    return(output)
  }


################## deephit_cv ########################

# The function trains deephit model with given params
#' @export
deephit_cv <- function(df,
                        predict.factors,
                        fixed_time = NaN,
                        outer_cv = 3,
                        inner_cv = 3,
                        repeat_cv = 2,
                        randomseed = NULL,
                        return_models = FALSE,
                        useCoxLasso = FALSE,
                        deephitparams = list()
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
    train_function = deephit_train,
    predict_function = deephit_predict,
    model_args = list("deephitparams" = deephitparams),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "DeepHit"
  )
  output$call <- Call
  return(output)
}
