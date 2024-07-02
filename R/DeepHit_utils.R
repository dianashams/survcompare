# 
# ################## Pre-tuned hyperparameters ################## 
# # Pre-tuned hyperparameters for the deephit
# datasets <- c("FUS", "ADNI", "HF", "W500", "GBSG2", "HNSCC", "ELSA")
# params <- c("dropout", "learning_rate", "num_nodes")
# dftune<- data.frame("datasets" = datasets)
# dftune$num_nodes <- 
#   list(c(16,16), c(32,32), c(32,32), c(32,32), c(16,16),c(16,16),c(16,16))
# dftune$dropout<- c(0.1, 0.3, .1, .3, .1, .1, .3)
# dftune$learning_rate<- c(0.1, 0.1, .3, .1, .1, .3, .1)
# dftune$epochs<- rep(10, length(datasets))
# dftune$batch_size <- c(100, 100, 100, 100, 100, 100, 250)
# row.names(dftune) = dftune$datasets
# deephitparams_tuned <- list()
# for (dt in datasets) {
#   deephitparams_tuned[[dt]]<- as.list(dftune[dt,])
#   deephitparams_tuned[[dt]]$num_nodes = dftune[dt,"num_nodes"][[1]]}
# 

################## deephit_predict ################## 
# The function to predict event probability by a trained deephitmodel
deephit_predict <-
  function(trained_model,
           newdata,
           predict.factors,
           fixed_time) {
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
deephit_train <-
  function(df_train,
           predict.factors,
           deephitparams = list(
             "dropout" = 0.2,
             "learning_rate" = 0.01,
             "num_nodes" = c(32, 32),
             "batch_size" = 100,
             "epochs" = 50,
             "frac"=0.1)
  ) {
    deephitm <- deephit(
      data = df_train,
      x = df_train[predict.factors],
      y = Surv(df_train$time, df_train$event),
      shuffle = TRUE,
      dropout = deephitparams$dropout,
      learning_rate = deephitparams$learning_rate,
      num_nodes = deephitparams$num_nodes,
      batch_size = deephitparams$batch_size,
      epochs = deephitparams$epochs,
      frac = deephitparams$frac
    )
    return(deephitm)
  }

################## deephit_tune ################## 
# This function repeats 3-fold CV repeat-tune times
# and finds grid with highest c-index
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
      "frac" = c(0.1, 0.3), #frac= c(0, 0.3)
      "num_nodes" =    list(c(8,8), c(16,16,16,2), c(32,32,32,4))
    ) #20
    grid_size <- dim(grid_of_values)[1]
    #placeholder for c-index
    cind = matrix(NA, nrow = grid_size, ncol = inner_cv)
    # tuning cross-validation loop
    
    #progress bar
    pb <- utils::txtProgressBar(0, inner_cv * grid_size, style = 3)
    utils::setTxtProgressBar(pb, inner_cv * grid_size / 50)
    
    for (cv_iteration in 1:inner_cv) {
      # print(paste("deephit tuning CV step", cv_iteration, "of", inner_cv))
      cv_folds <-
        caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
      df_train_cv <- df_tune[cv_folds != cv_iteration, ]
      df_test_cv <- df_tune[cv_folds == cv_iteration, ]
      #Grid search
      for (i in 1:grid_size) {
        deephitm <- deephit(
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
      }
      # here we already have cindex for allgrid for cv_iteration
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



################  ens_deephit_train ####################  
# Create out-of-bag Cox predictions, then train deephit
ens_deephit_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           inner_cv = 3,
           randomseed = NULL,
           useCoxLasso = FALSE,
           deephitparams = list()) {
    
    #setting random seed
    if (is.null(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    set.seed(randomseed)
    #creating folds
    cv_folds <-
      caret::createFolds(df_train$event, k = 10, list = FALSE)
    cindex_train <- vector(length = 10)
    cindex_test <- vector(length = 10)
    for (cv_iteration in 1:10) {
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
    
    # adding Cox predictions as a new factor to tune SRF,
    # taking out the first predictor to avoid collinearity 
    predict.factors.plusCox <- c(predict.factors[-1], "cox_predict")
    
    # train the deephit model
    deephit.ens <- 
      deephit_train(df_train,predict.factors.plusCox,deephitparams)
    
    #base cox model
    cox_base_model <-
      survcox_train(df_train, predict.factors, useCoxLasso = useCoxLasso)
    
    #output
    output = list()
    output$model <- deephit.ens
    output$model_base <- cox_base_model
    output$randomseed <- randomseed
    output$call <-  match.call()
    class(output) <- "survensemble"
    return(output)
  }

################  ens_deephit_predict ####################  

#same as deephit_predict
ens_deephit_predict <-
  function(trained_object,
           newdata,
           predict.factors,
           fixed_time) {
    
    trained_model<- trained_object$model
    predictdata <- newdata[predict.factors]
    # use model_base with the base Cox model to find cox_predict
    predictdata$cox_predict <- survcox_predict(trained_object$model_base,
                                               newdata, fixed_time)
    # now use "model" which is SRF which needs additional risk factor
    # "cox_predict" which was created in the previous row
    
    s1 <-
      predict(trained_model, newdata = predictdata, type = "survival")
    f <- function(i) {
      approxfun(as.double(colnames(s1)), s1[i,], method = "linear")(fixed_time)
    }
    predict_eventprob <- 1 - unlist(lapply(1:dim(s1)[1], f))
    
    return(predict_eventprob)
  }

############### ens_deephit_CV #############

ens_deepsurv_CV <- function(df,
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
    train_function = ens_deephit_train,
    predict_function = ens_deephit_predict,
    model_args = list("useCoxLasso" = useCoxLasso,
                      "deephitparams" = deephitparams),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "CoxPH and DeepHit Ensemble"
  )
  output$call <- Call
  return(output)
}

############### deephit_CV #############

deephit_CV <- function(df,
                        predict.factors,
                        fixed_time = NaN,
                        outer_cv = 3,
                        inner_cv = 3,
                        repeat_cv = 2,
                        randomseed = NULL,
                        return_models = FALSE,
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


# survcompare <- function(df_train,
#                         predict_factors,
#                         predict_time = NULL,
#                         randomseed = NULL,
#                         useCoxLasso = FALSE,
#                         outer_cv = 3,
#                         inner_cv = 3,
#                         srf_tuning = list(),
#                         return_models = FALSE,
#                         repeat_cv = 2,
#                         train_srf = FALSE)
#=>  cox_cv <- survcox_cv(
#   df = df_train,
#   predict.factors = predict_factors,
#   fixed_time = predict_time ,
#   outer_cv = outer_cv,
#   randomseed = randomseed,
#   useCoxLasso = useCoxLasso,
#   return_models = return_models,
#   repeat_cv = repeat_cv
# )
# 
# ens1_cv <- survensemble_cv(
#   df = df_train,
#   predict.factors = predict_factors,
#   fixed_time = predict_time,
#   outer_cv = outer_cv,
#   inner_cv = inner_cv,
#   randomseed = randomseed,
#   return_models = return_models,
#   srf_tuning = srf_tuning,
#   useCoxLasso = useCoxLasso,
#   repeat_cv = repeat_cv
# )
# Survensemble calls =>  output <- surv_CV(
#   df = df,
#   predict.factors = predict.factors,
#   fixed_time = fixed_time,
#   outer_cv = outer_cv,
#   inner_cv = inner_cv,
#   repeat_cv = repeat_cv,
#   randomseed = randomseed,
#   return_models = return_models,
#   train_function = survensemble_train,
#   predict_function = predict.survensemble,
#   model_args = list(
#     "useCoxLasso" = useCoxLasso,
#     "srf_tuning" = srf_tuning,
#     "oob" = oob
#   ),
#   model_name = "CoxPH and SRF Ensemble"
# )

eligible_params <- function(params, df) {
  # This function checks eligible predictors from params list for split
  # It deletes those which are
  # 1) not in df and
  # 2) taking only 1 value (constants)
  # TODOLater may delete collinear factors
  if (length(params) == 0) {
    return(NULL)
  }
  # take only columns which are in df
  z <- params %in% names(df)
  if (sum(!z) == length(params)) {
    return(NULL) # no eligible params
  } else {
    params <- params[z] # there are some potentially eligible
  }
  params_eligible <- params
  for (i in 1:length(params)) {
    if (length(unique(df[, params[i]])) < 2) {
      params_eligible <- params_eligible[params_eligible != params[i]]
    }
  }
  return(params_eligible)
}

surv_CV <-
  function(df,
           predict.factors,
           fixed_time = NaN,
           outer_cv = 3,
           inner_cv = 3,
           repeat_cv = 2,
           randomseed = NULL,
           return_models = FALSE,
           train_function,
           predict_function,
           model_args = list(),
           predict_args = list(),
           model_name = "a model") {
    time_0 <- Sys.time()
    if (is.null(randomseed)) {
      randomseed <- round(stats::runif(1) * 1e9, 0)
    }
    
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
    
    note_srf = ifelse(
      model_name == "Survival Random Forest",
      "For SRF inner CV is not used if oob = TRUE (default).",
      ""
    )
    
    print(
      paste(
        "Cross-validating ",
        model_name,
        " using ",
        repeat_cv,
        " repeat(s), ",
        outer_cv,
        " outer, ",
        inner_cv,
        " inner loops).",
        note_srf,
        sep = ""
      )
    )
    
    modelstats_train <- list()
    modelstats_test <- list()
    models_for_each_cv <- list()
    #progress bar
    pb <- utils::txtProgressBar(0, outer_cv * repeat_cv, style = 3)
    utils::setTxtProgressBar(pb, outer_cv * repeat_cv / 50)
    
    # cross-validation loop
    for (rep_cv in 1:repeat_cv) {
      set.seed(randomseed + rep_cv)
      if (rep_cv != 1) {
        df <- df[sample(1:nrow(df)), ]
      }
      cv_folds <-
        caret::createFolds(df$event, k = outer_cv, list = FALSE)
      
      for (cv_iteration in 1:outer_cv) {
        df_train_cv <- df[cv_folds != cv_iteration,]
        df_test_cv <- df[cv_folds == cv_iteration,]
        # tune the model using train_function
        trained_model <-
          do.call(train_function,
                  append(list(df_train_cv, predict.factors), model_args))
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
        if (return_models) {
          models_for_each_cv[[cv_iteration + (rep_cv - 1) * outer_cv]] <-
            trained_model
        }
        utils::setTxtProgressBar(pb, cv_iteration + (rep_cv - 1) * outer_cv)
      }
    }
    
    df_modelstats_test <- data.frame(modelstats_test[[1]])
    df_modelstats_train <- data.frame(modelstats_train[[1]])
    
    for (i in 2:(outer_cv * repeat_cv)) {
      df_modelstats_test[i, ] <- modelstats_test[[i]]
      df_modelstats_train[i, ] <- modelstats_train[[i]]
    }
    row.names(df_modelstats_train) <- 1:(outer_cv * repeat_cv)
    row.names(df_modelstats_test) <- 1:(outer_cv * repeat_cv)
    
    utils::setTxtProgressBar(pb, outer_cv * repeat_cv)
    close(pb)
    df_modelstats_test$test <- 1
    df_modelstats_train$test <- 0
    
    #summary for printing and summary(obj)
    stats_summary <- function(x) {
      #remove the last column ("test" 1 or 0)
      x = x[, 1:(dim(x)[2] - 1)]
      as.data.frame(round(
        cbind(
          "mean" = apply(x, 2, mean, na.rm = 1),
          "sd" = apply(x, 2, sd, na.rm = 1),
          "95CILow" = apply(x, 2, quantile, 0.025, na.rm = 1),
          "95CIHigh" = apply(x, 2, quantile, 0.975, na.rm = 1)
        ),
        4
      ))
    }
    #mean, sd and confidence intervals for all test and train datasets
    summarydf <- as.data.frame(cbind(
      "test" = stats_summary(df_modelstats_test),
      "train" = stats_summary(df_modelstats_train)
    ))
    
    output <- list()
    output$test <- df_modelstats_test
    output$train <- df_modelstats_train
    output$testaverage <-
      sapply(df_modelstats_test, mean, na.rm = 1)
    output$trainaverage <-
      sapply(df_modelstats_train, mean, na.rm = 1)
    output$tuned_cv_models <- models_for_each_cv
    output$randomseed <- randomseed
    output$call <- Call
    output$cv <- c(outer_cv, inner_cv, repeat_cv)
    output$summarydf <- summarydf
    class(output) <- "survensemble_cv"
    time_1 <- Sys.time()
    output$time <- time_1 - time_0
    print(time_1 - time_0)
    return(output)
  }
