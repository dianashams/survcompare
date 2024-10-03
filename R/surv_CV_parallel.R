#' @importFrom devtools load_all
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom foreach  %dopar%
#' @importFrom foreach  foreach
#' @importFrom reticulate use_python
#' @importFrom survivalmodels deephit
surv_CV_parallel =
  function(df,
           predict.factors,
           package_path=NaN,
           python_path=NaN,
           fixed_time = NaN,
           outer_cv = 3,
           inner_cv = 3,
           repeat_cv = 2,
           randomseed = NaN,
           return_models = FALSE,
           train_function,
           predict_function,
           model_args = list(),
           predict_args = list(),
           model_name = "my model") {

    if (is.nan(package_path)){
      stop("Please supply package path for parallel computations.")
    }

    time_0 <- Sys.time()
    print(package_path)

    # Make parallel cluster, use all but 2 cores for parallel CV
    num_cores <- parallel::detectCores()
    `%dopar%` <- foreach::`%dopar%`

    if (is.nan(randomseed)) {randomseed <- round(stats::runif(1) * 1e9, 0)}

    if (!any(
      is.data.frame(df),
      predict.factors %in% colnames(df),
      c("time", "event") %in% colnames(df)
    )) {
      stop("Data should be a data frame, predictors should correspond to the columns.")
    }

    if (sum(is.na(df[c("time", "event", predict.factors)])) > 0) {
      stop("Missing data can not be handled. Please impute first.")
    }

    if (sum(is.nan(fixed_time)) > 0) {
      fixed_time <- round(quantile(df[df$event == 1, "time"], 0.9), 1)
    }
    predict.factors <- eligible_params(predict.factors, df)

    if (length(predict.factors) == 0) {
      print("No eligible params")
      return(NULL)
    }
    Call <- match.call()

    #defining number of repeated cv
    if (is.null(repeat_cv)) {
      repeat_cv = 1
    }
    if (is.numeric(repeat_cv) & repeat_cv > 1) {
      repeat_cv = round(repeat_cv, 0)
    } else{
      repeat_cv = 1
    }

    note_srf = ifelse(model_name == "Survival Random Forest",
      "For SRF inner CV is not used if oob = TRUE (default)","")
    print(paste("Cross-validating with doParallel,",
                model_name," using ",
                repeat_cv," repeat(s), ",
                outer_cv," outer, ",
                inner_cv," inner loops).",
                note_srf,sep = ""))

    modelstats_train <- list()
    modelstats_test <- list()
    models_for_each_cv <- list()
    params_for_each_cv <- list()

    # repeat_cv loop
    for (rep_cv in 1:repeat_cv) {
      print(paste("Repeated CV", rep_cv, "/", repeat_cv))
      set.seed(randomseed + rep_cv)
      if (rep_cv != 1) {
        df <- df[sample(1:nrow(df)), ]
      }
      cv_folds <-
        caret::createFolds(df$event, k = outer_cv, list = FALSE)

      ############################################################
      # parallel cluster-register
      cl= parallel::makeCluster(num_cores-2)
      parallel::clusterExport(cl, varlist = c("package_path", "python_path"))
      doParallel::registerDoParallel(cl)
      # cross-validation loop:
      results <- foreach::foreach(
        cv_iteration = 1:outer_cv
        ) %dopar% {
          devtools::load_all(package_path)
          if(!is.nan(python_path)) {reticulate::use_python(python_path)}

        df_train_cv <- df[cv_folds != cv_iteration,]
        df_test_cv <- df[cv_folds == cv_iteration,]
        predict.factors.cv <- eligible_params(predict.factors, df_train_cv)

        # tune the model using train_function
        trained_model <-
          do.call(
            train_function,
            append(list(df_train_cv, predict.factors.cv), model_args))

        # compute event probability predictions for given times
        y_predict_test <-
          do.call(
            predict_function,
            append(list(trained_model, df_test_cv, fixed_time),predict_args))
        y_predict_train <-
          do.call(
            predict_function,
            append(list(trained_model, df_train_cv, fixed_time),predict_args))

        modelstats_test <-
          surv_validate(y_predict_test,fixed_time,df_train_cv,df_test_cv)
        modelstats_train <-
          surv_validate(y_predict_train,fixed_time,df_train_cv,df_train_cv)

        # save models and tuned parameters
        if (return_models) {  models_for_each_cv <- trained_model}
        if(!is.null(trained_model$bestparams)){
          params_for_each_cv<-trained_model$bestparams
          }else{params_for_each_cv<-NaN}

        # return from each parallel run of CV fold
        list(modelstats_train=modelstats_train,
             modelstats_test=modelstats_test,
             models_for_each_cv=models_for_each_cv,
             params_for_each_cv=params_for_each_cv)

      } #end of cv loop

      for (i in 1:length(results)){
        modelstats_train[[i + (rep_cv - 1) * outer_cv]] = results[[i]]$modelstats_train
        modelstats_train[[i + (rep_cv - 1) * outer_cv]][, "repeat_cv"] = rep_cv
        modelstats_train[[i + (rep_cv - 1) * outer_cv]][, "outer_cv"] = i
        modelstats_test[[i + (rep_cv - 1) * outer_cv]] = results[[i]]$modelstats_test
        modelstats_test[[i + (rep_cv - 1) * outer_cv]][, "repeat_cv"] = rep_cv
        modelstats_test[[i + (rep_cv - 1) * outer_cv]][, "outer_cv"] = i
        models_for_each_cv[[i + (rep_cv - 1) * outer_cv]] = results[[i]]$models_for_each_cv
        params_for_each_cv[[i + (rep_cv - 1) * outer_cv]] = results[[i]]$params_for_each_cv
      }
      parallel::stopCluster(cl)

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
      bestparams$C_score_outer = df_modelstats_test$C_score}

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
    output$testaverage <-
      sapply(df_modelstats_test, mean, na.rm = TRUE)
    output$testmedian <-
      sapply(df_modelstats_test, median, na.rm = TRUE)
    output$trainaverage <-
      sapply(df_modelstats_train, mean, na.rm = TRUE)
    output$tuned_cv_models <- models_for_each_cv
    output$randomseed <- randomseed
    output$bestparams<- bestparams
    output$call <- Call
    output$cv <- c(outer_cv, inner_cv, repeat_cv)
    output$summarydf <- summarydf
    output$summarydf_pooled <- summarydf_pooled
    output$model_name <- model_name
    class(output) <- "survensemble_cv"
    time_1 <- Sys.time()
    output$time <- difftime(time_1, time_0, units = "secs")
    print(time_1 - time_0)
    return(output)
  }
