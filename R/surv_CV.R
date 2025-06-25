surv_CV <-
  function(df,
           predict.factors,
           fixed_time = NaN,
           outer_cv = 5,
           inner_cv = 3,
           repeat_cv = 2,
           randomseed = NaN,
           return_models = FALSE,
           train_function,
           predict_function,
           model_args = list(),
           predict_args = list(),
           model_name = "my model",
           impute = 0,
           impute_method = "missForest") {
    time_0 <- Sys.time()

    if (is.nan(fixed_time)) {fixed_time <-quantile(df[df$event == 1, "time"], 0.9, na.rm = TRUE)}
    if (is.nan(randomseed)) {randomseed <- round(stats::runif(1) * 1e9, 0)}
    if (!is.null(model_args$max_grid_size)) { 
      max_grid_size = model_args$max_grid_size } else{ max_grid_size = NaN }
    
    if (!any(is.data.frame(df),
      predict.factors %in% colnames(df),
      c("time", "event") %in% colnames(df)
    )) {
      stop("Data should be a data frame, predictors
           should correspond to the columns.")
    }
    predict.factors <- eligible_params(predict.factors, df)
    if (length(predict.factors) == 0) {
      print("No eligible params")
      return(NULL)
    }
    Call <- match.call()
    
    # Missing data handling
    if (sum(is.na(df[c("time", "event", predict.factors)])) > 0){ 
      df = missing_data_handle(df, predict.factors, impute)
      }
    
    #defining number of repeated cv
    if (is.null(repeat_cv)) {      repeat_cv = 1    }
    if (is.numeric(repeat_cv) & repeat_cv > 1) {
      repeat_cv = round(repeat_cv, 0)    } else{      repeat_cv = 1    }

    print(paste("Cross-validating ",model_name," using ",repeat_cv,
        " repeat(s), ",outer_cv," outer, ",inner_cv," inner loops).",
        sep = ""))

    modelstats_train <- list()
    modelstats_test <- list()
    models_for_each_cv <- list()
    params_for_each_cv <- list()

    # repeat_cv loop
    for (rep_cv in 1:repeat_cv) {
      print(paste("Repeated CV", rep_cv, "/", repeat_cv))
      set.seed(randomseed + rep_cv)
      if (rep_cv != 1) {  df <- df[sample(1:nrow(df)), ]}
      cv_folds <- caret::createFolds(df$event, k = outer_cv, list = FALSE)

      # cross-validation loop:
      pb <- utils::txtProgressBar(0, outer_cv, style = 3) # progress bar
      for (cv_iteration in 1:outer_cv) {
        utils::setTxtProgressBar(pb, cv_iteration) #progress bar update
        df_train_cv <- df[cv_folds != cv_iteration,]
        df_test_cv <- df[cv_folds == cv_iteration,]
        
        # impute train and test if impute == 1 and there are missing values 
        if (impute == 1){
          temp = impute1(df_train_cv, df_test_cv, predict.factors)
          df_train_cv = temp$train
          df_test_cv = temp$test
          remove(temp)
        }
        predict.factors.cv <- eligible_params(predict.factors, df_train_cv)
        
        # tune the model using train_function
        trained_model <-
          do.call(train_function,
                  append(list(df_train_cv, predict.factors.cv), model_args))

        # compute event probability predictions for given times
        y_predict_test <-
          do.call(predict_function,
                  append(
                    list(trained_model, df_test_cv, fixed_time),
                    predict_args
                  ))
        y_predict_train <-
          do.call(predict_function,
                  append(
                    list(trained_model, df_train_cv, fixed_time),
                    predict_args
                  ))
        modelstats_test[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
          surv_validate(y_predict_test,
                        fixed_time,
                        df_train_cv,
                        df_test_cv,
                        weighted = 1)
        modelstats_train[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
          surv_validate(y_predict_train,
                        fixed_time,
                        df_train_cv,
                        df_train_cv,
                        weighted = 1)

        modelstats_train[[
          cv_iteration + (rep_cv - 1) * outer_cv]][, "repeat_cv"] = rep_cv
        modelstats_train[[
          cv_iteration + (rep_cv - 1) * outer_cv]][, "outer_cv"] = cv_iteration
        modelstats_test[[
          cv_iteration + (rep_cv - 1) * outer_cv]][, "repeat_cv"] = rep_cv
        modelstats_test[[
          cv_iteration + (rep_cv - 1) * outer_cv]][, "outer_cv"] = cv_iteration

        if (return_models) {
          models_for_each_cv[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
            trained_model
          if(!is.null(trained_model$bestparams)){
          params_for_each_cv[[cv_iteration + (rep_cv - 1) * outer_cv]]<-
            trained_model$bestparams}
        }else{
          # save tuned parameters
          if(!is.null(trained_model$bestparams))
            {params_for_each_cv[[cv_iteration + (rep_cv - 1) * outer_cv]]<-
              trained_model$bestparams}#end if
          }#end else
      } #end of cv loop
      close(pb) #close progress bar to start new one
    }#end of repeat loop

    # glue the list into a data frame
    df_modelstats_test = as.data.frame(do.call(rbind, modelstats_test))
    df_modelstats_train = as.data.frame(do.call(rbind, modelstats_train))
    row.names(df_modelstats_train) <- 1:(outer_cv * repeat_cv)
    row.names(df_modelstats_test) <- 1:(outer_cv * repeat_cv)
    df_modelstats_test$test <- 1
    df_modelstats_train$test <- 0

    bestparams = as.data.frame(do.call(rbind, params_for_each_cv))
    if (dim(bestparams)[1]==outer_cv*repeat_cv) {
      bestparams$C_score_outer = df_modelstats_test$C_score
      bestparams = bestparams[order(-bestparams$C_score_outer), ]
      }

    #summary for printing and summary(obj)
    stats_summary <- function(x) {
      #remove the last column ("test" 1 or 0)
      x = x[, 1:(dim(x)[2] - 1)]
      as.data.frame(round(
        cbind(
          "mean" = apply(x, 2, mean, na.rm = TRUE),
          "sd" = apply(x, 2, sd, na.rm = TRUE),
          "95CILow" = apply(x, 2, quantile, 0.025, na.rm = TRUE),
          "95CIHigh" = apply(x, 2, quantile, 0.975, na.rm = TRUE),
          "median" = apply(x, 2, median, na.rm = TRUE)
        ),4))
    }
    #if CV was repeated many times, get pooled CV stats too
    pooled_test <- function(x) {
      temp2<- data.frame()
      for (i in 1:repeat_cv) {
        shift = (i-1)*outer_cv #get a chunk of x for i-th repetition
        xi= apply(x[(shift+1):(shift+outer_cv), ], 2, mean,na.rm=TRUE)
        temp2<- rbind(temp2, xi)
      }
      names(temp2) = names(x)
      return(temp2)
    }
    #mean, sd and confidence intervals for all test and train datasets
    summarydf <- as.data.frame(cbind(
      "test" = stats_summary(df_modelstats_test),
      "train" = stats_summary(df_modelstats_train)
    ))
    #Same for the pooled results over CV repetitions (if >1)
    if (repeat_cv>1){
      summarydf_pooled <- as.data.frame(cbind(
        "test" = stats_summary(pooled_test(df_modelstats_test)),
        "train" = stats_summary(pooled_test(df_modelstats_train))
      ))
    }else{summarydf_pooled = summarydf}

    output <- list()
    output$test <- df_modelstats_test
    output$train <- df_modelstats_train
    output$test_pooled <- pooled_test(df_modelstats_test)
    output$testaverage <- sapply(df_modelstats_test, mean, na.rm = TRUE)
    output$testmedian <- sapply(df_modelstats_test, median, na.rm = TRUE)
    output$trainaverage <- sapply(df_modelstats_train, mean, na.rm = TRUE)
    output$tuned_cv_models <- models_for_each_cv
    output$randomseed <- randomseed
    output$bestparams<- bestparams
    output$call <- Call
    output$cv <- c("outer_cv" = outer_cv, "inner_cv"=inner_cv,
                   "repeat_cv"=repeat_cv, "max_grid_size"=max_grid_size)
    output$summarydf <- summarydf
    output$summarydf_pooled <- summarydf_pooled
    output$model_name <- model_name
    class(output) <- "survensemble_cv"
    time_1 <- Sys.time()
    output$time <- difftime(time_1, time_0, units = "secs")
    print(time_1 - time_0)
    return(output)
  }

impute1 <- 
  function(df_train, df_test, predict.factors) {
  miss_train = sum(is.na(df_train[predict.factors]))
  miss_test = sum(is.na(df_test[ predict.factors]))
  if (miss_train > 0 | miss_test > 0) {
    mf = missForestPredict::missForest(
      df_train[predict.factors], save_models = TRUE,
      num.trees = 100, maxiter = 5, verbose = FALSE)
    df_train[predict.factors] = mf$ximp
    if (miss_test > 0) {
      df_test[predict.factors] = 
        missForestPredict::missForestPredict(
          mf, newdata = df_test[predict.factors])
    }
  }
  output = list()
  output$train = df_train
  output$test = df_test
  return (output)
}

missing_data_handle <- function(df, predict.factors, impute){
  msg0 <- "Set 'impute' to 0 for no imputation, 1 for proper imputation by missForest, 2 for fast imputation by missForest, 3 for complete cases."
  # If impute = 0, then there is no imputation, we stop with a message.
  if (impute == 0) {
    stop(paste("Some data are missing. By default ('impute' = 0), no imputation is performed.", msg0))
  }
  # If impute = 3, then analyse complete cases.
  if (impute == 3) {
    msg3 =  "Some data are missing, 'impute' is set to 3 for complete cases analysis."
    print(paste(msg3, msg0))
    print(paste("Share of complete cases = ", round(mean(complete.cases(df)),2)))
    df <- df[complete.cases(df), ]
    return(df)
  }
  # If impute = 2, then we impute first, then perform all validations
  if (impute == 2) {
    msg2 <- "Some data are missing, 'impute' is set to 2 for fast imputation (impute all the data, then cross-validate)."
    print(paste(msg2,msg0))
    temp <-
      missForestPredict::missForest(
        df[c("time", "event", predict.factors)],
        save_models = FALSE,
        maxiter = 5,
        num.trees = 100,
        verbose = FALSE
      )
    df = temp$ximp
    remove(temp)
    return(df)
  }
  # If impute = 2, then we impute properly, train, then test. Only information message here
  if (impute == 1) {
    msg1 <- "Some data are missing, 'impute' is set to 1 for proper imputation (impute train, then test for each data split)."
    print(paste(msg1,msg0))
    return(df)
  }
  return(df)
}