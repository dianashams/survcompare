

################## deepsurv_predict ##################
# The function to predict event probability by a trained deepsurvmodel
#' @export
deepsurv_predict <-
  function(trained_model,
           newdata,
           fixed_time,
           predict.factors)
  {
    if (inherits(trained_model, "list")) {trained_model <- trained_model$model}
    if (class(try(
      predict(trained_model, newdata[predict.factors], type = "survival"),
      silent= TRUE))[1]=="try-error"){
      return(        rep(NaN, dim(newdata)[1])              )
    }
    s1 <- predict(trained_model, newdata[predict.factors], type = "survival")
    f <- function(i) {
      if (class(try(
        approxfun(as.double(colnames(s1)), s1[i, ], method = "linear")(fixed_time),
        silent = TRUE   )      )  ==  "try-error") {return(NaN)
        } else{
          return(
            approxfun(as.double(colnames(s1)), s1[i, ], method = "linear")(fixed_time)
            )
      }
    }
    predict_eventprob <- 1 - unlist(lapply(1:dim(s1)[1], f))
    return(predict_eventprob)
  }

################## deepsurv_train ##################
# The function trains deepsurv model with given params
#' @export
deepsurv_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           tuningparams = list(),
           max_grid_size = 10,
           inner_cv = 3,
           randomseed = NaN,
           verbose = TRUE) {

    grid_of_hyperparams <-
      ml_hyperparams_deepsurv(
        mlparams = tuningparams,
        dftune_size = dim(df_train)[1],
        max_grid_size = max_grid_size
      )

    if (dim(grid_of_hyperparams)[1] == 1) {
      #if only one option is given, we use it to train
      bestparams = grid_of_hyperparams[1,]
    } else{
      #if there are options, we tune the grid of hyperparams
      tuning_m <- deepsurv_tune(
        df_train,
        predict.factors,
        repeat_tune = 1,
        fixed_time = fixed_time,
        tuningparams = tuningparams,
        max_grid_size = max_grid_size,
        inner_cv = inner_cv,
        randomseed = randomseed
      )
      bestparams = tuning_m$bestparams
    }
    if (verbose) {
      cat("\n") # !!!!!!!!!!!!!!!!!!!!!!!!!! #
      print(bestparams) # !!!!!!!!!!!!!!!!!!!!!!!!!! #
    }
    deepsurvm <- deepsurv(
      data = df_train,
      x = df_train[predict.factors],
      y = Surv(df_train$time, df_train$event),
      early_stopping = bestparams$early_stopping,
      dropout = bestparams$dropout,
      learning_rate = bestparams$learning_rate,
      num_nodes = bestparams$num_nodes[[1]],
      batch_size = bestparams$batch_size,
      epochs = bestparams$epochs,
      weight_decay = bestparams$weight_decay,
      batch_norm = bestparams$batch_norm
    )
    output = list()
    output$model_name = "DeepSurv"
    output$model = deepsurvm
    output$bestparams = bestparams
    return(output)
  }

################## deepsurv_tune ##################
#'A repeated 3-fold CV over a hyperparameters grid
#'
#' @param repeat_tune number of repeats
#' @param df_tune data
#' @param predict.factors predictor names
#' @param fixed_time  not used here, but for some models the time for which performance is optimized
#' @param tuningparams list with the range of hyperparameters
#' @param inner_cv number of folds for each CV
#' @param max_grid_size 25 by default. If grid size > max_grid_size, a random search is performed for max_grid_size iterations. Set this to a small number for random search
#' @param randomseed to choose random subgroup of hyperparams
#' @return  output=list(cindex_ordered, bestparams)
#' @examples
#' d <-simulate_nonlinear(100),
#' p<- names(d)[1:4]
#' tuningparams = list(
#'  dropout = c(0.1, 0.3),
#'  learning_rate = c(0.001),
#'  num_nodes =    list(c(8, 8), c(16, 16, 16, 4), c(32, 32, 32, 4)),
#'  batch_size = min(max(64, round(dim(df_train_cv)[1]/8, 0)),256),
#'  epochs = c(5, 50),
#'  weight_decay = 0.1
#' )
#' deepsurv_tune(d, p,tuningparams = tuningparams, max_grid_size = 10)
#' @export
deepsurv_tune <-
  function(df_tune,
           predict.factors,
           repeat_tune=1,
           fixed_time = NaN,
           tuningparams = list(),
           max_grid_size = 10,
           inner_cv = 3,
           randomseed = NaN) {

    deepsurvgrid <-
      ml_hyperparams_deepsurv(
        mlparams = tuningparams,
        dftune_size = dim(df_tune)[1],
        max_grid_size = max_grid_size,
        randomseed = randomseed
      )
    means = c()
    for (i in (1:repeat_tune)) {
      ds_tune_fus <-
        deepsurv_tune_single(
          df_tune = df_tune,
          predict.factors =  predict.factors,
          fixed_time = fixed_time,
          grid_hyperparams = deepsurvgrid,
          inner_cv = inner_cv,
          randomseed=randomseed
        )
      means <- cbind(means, ds_tune_fus$cindex_mean)
      remove(ds_tune_fus)
    }

    allmeans <- apply(means, 1, mean)
    # grid and average cindex ordered by cindex (highest first)
    output = list()
    output$cindex_ordered <-
      cbind(deepsurvgrid, allmeans)[order(allmeans, decreasing = TRUE),]
    output$bestparams <- output$cindex_ordered[1,]
    return(output)
  }

#' Internal function for getting grid of hyperparameters
#' for random or grid search of size = max_grid_size
#' @export
ml_hyperparams_deepsurv <- function(mlparams = list(),
                           max_grid_size = 10,
                           dftune_size = 1000,
                           randomseed = NaN) {
  default_grid <- list(
      "dropout" = c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),
      "learning_rate" = c(0.001,0.01,0.1),
      "num_nodes" = list(
        c(16, 16),  c(32, 32),  c(64, 64), c(128,128),
        c(8, 8, 8),  c(16, 16, 16), c(32,32,32),c(64,64,64),
        c(16, 16, 16, 16), c(32, 32, 32, 32), c(64, 64, 64, 64)),
      "batch_size" =
        seq(min(64, round(dftune_size/ 8)), 256, 64),
      "epochs" = c(10, 50, 100, 200),
      "early_stopping" = c(TRUE, FALSE),
      "weight_decay" = c(0,0.1,0.2,0.3,0.4,0.5),
      "batch_norm" = TRUE
    )

  if (length(mlparams) == 0) {
    grid_of_hyperparams <- expand.grid(default_grid)
  } else{
    if (is.null(mlparams$dropout)) {mlparams$dropout = default_grid$dropout    }
    if (is.null(mlparams$learning_rate)) {mlparams$learning_rate = default_grid$learning_rate}
    if (is.null(mlparams$num_nodes)) {mlparams$num_nodes = default_grid$num_nodes}
    if (is.null(mlparams$batch_size)) {mlparams$batch_size = default_grid$batch_size}
    if (is.null(mlparams$epochs)) {mlparams$epochs = default_grid$epochs}
    if (is.null(mlparams$early_stopping)) {mlparams$early_stopping = default_grid$early_stopping}
    if (is.null(mlparams$weight_decay)) {mlparams$weight_decay = default_grid$weight_decay}
    if (is.null(mlparams$batch_norm)) {mlparams$batch_norm = default_grid$batch_norm}

    grid_of_hyperparams <- expand.grid(
      "dropout" = mlparams$dropout,
      "learning_rate" = mlparams$learning_rate,
      "num_nodes" = mlparams$num_nodes,
      "batch_size" = mlparams$batch_size,
      "epochs" = mlparams$epochs,
      "early_stopping"= mlparams$early_stopping,
      "weight_decay" = mlparams$weight_decay,
      "batch_norm" = mlparams$batch_norm
    )
  }

  grid_size <- dim(grid_of_hyperparams)[1]

  final_grid <- grid_of_hyperparams
  if (grid_size > max_grid_size) {
    if (!is.nan(randomseed)) {set.seed(randomseed)}
    final_grid <-
      grid_of_hyperparams[sample(1:grid_size, max_grid_size, replace = FALSE), ]
    grid_size <- max_grid_size
  }
  remove(grid_of_hyperparams); remove(default_grid)

  return(final_grid)
}


################## deepsurv_tune_single ##################
#' Internal function for deepsurv_tune(), performs 1 CV
#' @param df_tune data
#' @param predict.factors predictor names
#' @param fixed_time predictions for which time are computed for c-index
#' @param grid_hyperparams  hyperparameters grid (or a default will be used )
#' @param inner_cv number of folds for each CV
#' @return  output=list(grid, cindex, cindex_mean)
#' @export
deepsurv_tune_single <-
  function(df_tune,
           predict.factors,
           fixed_time = NaN,
           grid_hyperparams=c(),
           inner_cv = 3,
           randomseed = NaN) {

    #fixed_time
    if (is.nan(fixed_time) | length(fixed_time) > 1) {
      # not implemented for multiple time
      fixed_time <-
        round(quantile(df_tune[df_tune$event == 1, "time"], 0.9), 2)
    }
    if (length(grid_hyperparams) == 0) {
      grid_hyperparams <-
        ml_hyperparams_deepsurv(mlparams = list(),
                       dftune_size = dim(df_tune)[1])
    }
    grid_size <- dim(grid_hyperparams)[1]
    #placeholder for c-index
    cind = matrix(NA, nrow = grid_size, ncol = inner_cv)
    #progress bar
    pb0 <- utils::txtProgressBar(0, inner_cv * grid_size, style = 3)
    utils::setTxtProgressBar(pb0, inner_cv * grid_size / 50)

    # tuning cross-validation loop
    for (cv_iteration in 1:inner_cv) {
      if (!is.nan(randomseed)) {
        set.seed(randomseed + 123 + cv_iteration)
      }
      # print(paste("deepsurv tuning CV step", cv_iteration, "of", inner_cv))
      cv_folds <-
        caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
      df_train_cv <- df_tune[cv_folds != cv_iteration,]
      df_test_cv <- df_tune[cv_folds == cv_iteration,]
      #Grid search
      for (i in 1:grid_size) {
        deepsurvm <- survivalmodels::deepsurv(
          data = df_train,
          x = df_train_cv[predict.factors],
          y = Surv(df_train_cv$time, df_train_cv$event),
          shuffle = TRUE,
          early_stopping = grid_hyperparams[i, "early_stopping"],
          dropout = grid_hyperparams[i, "dropout"],
          learning_rate = grid_hyperparams[i, "learning_rate"],
          num_nodes = grid_hyperparams[i, "num_nodes"][[1]],
          batch_size = grid_hyperparams[i, "batch_size"],
          epochs = grid_hyperparams[i, "epochs"],
          weight_decay = grid_hyperparams[i, "weight_decay"],
          batch_norm = grid_hyperparams[i, "batch_norm"]
        )
        #check test performance
        pp <-
          deepsurv_predict(
            trained_model = deepsurvm,
            newdata = df_test_cv,
            predict.factors = predict.factors,
            fixed_time = fixed_time
          )
        cind[i, cv_iteration] =
          surv_validate(pp, fixed_time, df_train_cv, df_test_cv)[1, "C_score"]
        #cind[i, cv_iteration] =
        #  concordance(Surv(df_test_cv$time, df_test_cv$event) ~ pp)$concordance
        utils::setTxtProgressBar(pb0, grid_size + (i - 1) * cv_iteration)
      } # here we already have cindex for allgrid for cv_iteration
    }
    utils::setTxtProgressBar(pb0, inner_cv * grid_size)
    close(pb0)
    remove(deepsurvm)
    cindex_mean <- apply(cind, 1, mean)
    output = list()
    output$cindex = cind
    output$cindex_mean = cindex_mean
    output$cindex_ordered <-
      cbind(grid_hyperparams, cindex_mean)[order(cindex_mean, decreasing = TRUE), ]
    output$bestparams <- output$cindex_ordered[1, ]
    return(output)
  }


################## deepsurv_cv ########################

# The function trains deepsurv model with given params
#' @export
deepsurv_cv <- function(df,
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
    train_function = deepsurv_train,
    predict_function = deepsurv_predict,
    model_args = list("tuningparams" = tuningparams,
                      "max_grid_size"= max_grid_size),
    predict_args = list("predict.factors" = predict.factors),
    model_name = "DeepSurv",
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
