
################## survsrf_predict ##################
# The function to predict event probability by a trained srfmodel
#' @export
survsrf_predict <-
  function(trained_model,
           newdata,
           fixed_time
  )
  {
    if (inherits(trained_model, "list")) {trained_model <- trained_model$model}
    predictions <- predict(trained_model, newdata = newdata)
    train_times <- predictions$time.interest
    s1 <- predictions$survival
    if(fixed_time > max(train_times)) {return(rep(NaN, dim(newdata)[1]))}
    f <- function(i) {
      approxfun(train_times, s1[i, ], method = "constant")(fixed_time)}
    predict_eventprob <- 1 - unlist(lapply(1:dim(s1)[1], f))
    return(predict_eventprob)
  }

################## survsrf_train ##################
# The function trains srf model with given params
#' @export
survsrf_train <-
  function(df_train,
           predict.factors,
           fixed_time = NaN,
           tuningparams = list(),
           max_grid_size = 10,
           inner_cv = 3,
           randomseed = NaN,
           verbose = TRUE) {
    
    grid_of_hyperparams <-
      ml_hyperparams_srf(
        mlparams = tuningparams,
        dftune_size = dim(df_train)[1],
        max_grid_size = max_grid_size,
        randomseed = randomseed,
        p = length(predict.factors)
      )
    if (dim(grid_of_hyperparams)[1] == 1) {
      #if only one option is given, we use it to train
      bestparams = grid_of_hyperparams[1,]
    } else{
      #if there are options, we tune the grid of hyperparams
      tuning_m <- survsrf_tune(
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
    if (verbose) {print(bestparams)}
    
    srfm <- randomForestSRC::rfsrc(
      as.formula(
        paste("Surv(time, event) ~",
              paste(predict.factors, collapse = "+"))),
      data = df_train,
      nodesize = as.integer(bestparams$nodesize),
      #  AVERAGE size, so we want this to be quite high
      ntree = 500,
      mtry = as.integer(bestparams$mtry),
      nodedepth = as.integer(bestparams$nodedepth),
      nsplit = 50,
      splitrule = "logrank",
      statistics = FALSE,
      membership = TRUE,
      importance = "none",# speed up, switching off VIMP calculations
      seed = randomseed
    )
    output = list()
    output$model_name = "Survival Random Forest"
    output$model = srfm
    output$bestparams = bestparams
    return(output)
  }

################## survsrf_tune ##################
#'A repeated 3-fold CV over a hyperparameters grid
#'
#' @param repeat_tune number of repeats
#' @param df_tune data
#' @param predict.factors predictor names
#' @param fixed_time  not used here, but for some models the time for which performance is optimized
#' @param tuningparams list with the range of hyperparameters
#' @param inner_cv number of folds for each CV
#' @param max_grid_size 10 by default. If grid size > max_grid_size, a random search is performed for max_grid_size iterations. Set this to a small number for random search
#' @param randomseed to choose random subgroup of hyperparams
#' @return  output=list(cindex_ordered, bestparams)
#' @examples
#' d <-simulate_nonlinear(100),
#' p<- names(d)[1:4]
#' tuningparams = list(
#'  "dropout" = c(0.1, 0.3),
#'  "learning_rate" = c(0.001),
#'  "num_nodes" =    list(c(8, 8), c(16, 16, 16, 4), c(32, 32, 32, 4)),
#'  "batch_size" = min(max(64, round(dim(df_train_cv)[1]/8, 0)),256),
#'  "epochs" = c(5, 50),
#'  "mod_alpha" = 0.2,
#'  "sigma" = 0.1,
#'  "cuts" = 10
#' )
#' srf_tune(d, p,tuningparams = tuningparams, max_grid_size = 10)
#' @export
survsrf_tune <-
  function(df_tune,
           predict.factors,
           repeat_tune=1,
           fixed_time = NaN,
           tuningparams = list(),
           max_grid_size = 10,
           inner_cv = 3,
           randomseed = NaN) {
    
    srfgrid <-
      ml_hyperparams_srf(
        mlparams = tuningparams,
        dftune_size = dim(df_tune)[1],
        max_grid_size = max_grid_size,
        randomseed = randomseed,
        p = length(predict.factors)
      )
    means = c()
    for (i in (1:repeat_tune)) {
      ds_tune_fus <-
        survsrf_tune_single(
          df_tune = df_tune,
          predict.factors =  predict.factors,
          fixed_time = fixed_time,
          grid_hyperparams = srfgrid,
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
      cbind(srfgrid, allmeans)[order(allmeans, decreasing = TRUE),]
    output$bestparams <- output$cindex_ordered[1,]
    return(output)
  }


################## survsrf_tune_single ##################
#' Internal function for survsrf_tune(), performs 1 CV
#' @param df_tune data
#' @param predict.factors predictor names
#' @param fixed_time predictions for which time are computed for c-index
#' @param grid_hyperparams  hyperparameters grid (or a default will be used )
#' @param inner_cv number of folds for each CV
#' @return  output=list(grid, cindex, cindex_mean)
#' @export
survsrf_tune_single <-
  function(df_tune,
           predict.factors,
           fixed_time = NaN,
           grid_hyperparams=c(),
           inner_cv = 3,
           randomseed = NaN,
           progressbar = FALSE) {
    #fixed_time
    if (is.nan(fixed_time) | length(fixed_time) > 1) {
      # not implemented for multiple time
      fixed_time <-
        round(quantile(df_tune[df_tune$event == 1, "time"], 0.9), 2)
    }
    if (length(grid_hyperparams) == 0) {
      grid_hyperparams <-
        ml_hyperparams_srf(mlparams = list(),
                           dftune_size = dim(df_tune)[1],
                           p = length(predict.factors),
                           randomseed = randomseed)
    }
    grid_size <- dim(grid_hyperparams)[1]
    #placeholder for c-index
    cind = matrix(NA, nrow = grid_size, ncol = inner_cv)
    #progress bar
    if(progressbar) {
      pb <- utils::txtProgressBar(0, inner_cv * grid_size, style = 3)
      utils::setTxtProgressBar(pb, inner_cv * grid_size / 50)
    }
    # tuning cross-validation loop
    for (cv_iteration in 1:inner_cv) {
      if (!is.nan(randomseed)) {
        set.seed(randomseed + 123 + cv_iteration)
      }

      cv_folds <-
        caret::createFolds(df_tune$event, k = inner_cv, list = FALSE)
      df_train_cv <- df_tune[cv_folds != cv_iteration,]
      df_test_cv <- df_tune[cv_folds == cv_iteration,]
      #Grid search
      for (i in 1:grid_size) {
        srfm <- 
          randomForestSRC::rfsrc(
            as.formula(paste(
              "Surv(time, event) ~",
              paste(predict.factors, collapse = "+")
            )),
            data = df_train_cv,
            nodesize =  grid_hyperparams[i, "nodesize"],
            # this is AVERAGE size, so we want this to be quite high
            ntree = 300,
            mtry = grid_hyperparams[i, "mtry"],
            nodedepth = grid_hyperparams[i, "nodedepth"],
            nsplit = 50,
            splitrule = "logrank",
            statistics = FALSE,
            membership = TRUE,
            importance = "none"
          )
        #check test performance
        pp <-
          survsrf_predict(
            trained_model = srfm,
            newdata = df_test_cv,
            fixed_time = fixed_time
          )
        cind[i, cv_iteration] =
          surv_validate(pp, fixed_time, df_train_cv, df_test_cv)[1, "C_score"]
        if (progressbar) {
          utils::setTxtProgressBar(pb, grid_size + (i - 1) * cv_iteration)
        }
      } # here we already have cindex for allgrid for cv_iteration
    }
    if(progressbar) {
      utils::setTxtProgressBar(pb, inner_cv * grid_size)
      close(pb)
    }
    remove(srfm)
    cindex_mean <- apply(cind, 1, mean)
    output = list()
    output$cindex = cind
    output$cindex_mean = cindex_mean
    output$cindex_ordered <-
      cbind(grid_hyperparams, cindex_mean)[order(cindex_mean, decreasing = TRUE), ]
    output$bestparams <- output$cindex_ordered[1, ]
    return(output)
  }

######################## ml_hyperparams_srf ################################

#' Internal function for getting grid of hyperparameters
#' for random or grid search of size = max_grid_size
#' @export
ml_hyperparams_srf <- function(mlparams = list(),
                              p = 10,
                              max_grid_size = 10,
                              dftune_size = 1000,
                              randomseed = NaN) {
  mtry = ifelse(
      p <= 10, c(2, 3, 4, 5),
      ifelse(p <= 25, c(3, 5, 7, 10, 15),
             round(c(p/10, p/5, p/3, p/2, 
                     max(1, round(sqrt(p), 0))), 0)))
  nodesize <- seq(5, 50, 10)
  nodedepth <- c(5, 25)
  default_grid <- list(mtry = mtry, nodesize = nodesize, nodedepth = nodedepth)
  
  if (length(mlparams) == 0) {
    grid_of_hyperparams <- expand.grid(default_grid)
  } else{
    if (is.null(mlparams$mtry)) {mlparams$mtry = default_grid$mtry}
    if (is.null(mlparams$nodesize)) {mlparams$nodesize = default_grid$nodesize}
    if (is.null(mlparams$nodedepth)) {mlparams$nodedepth = default_grid$nodedepth}
    
    mlparams$mtry <- mlparams$mtry[mlparams$mtry <= p]
    if (length(mlparams$mtry)==0 ) {mlparams$mtry <- c(max(1,round(sqrt(p),0)))}

    grid_of_hyperparams <- expand.grid(
      "mtry" = mlparams$dropout,
      "nodesize" = mlparams$learning_rate,
      "nodedepth" = mlparams$num_nodes
      )
  }
  grid_size <- dim(grid_of_hyperparams)[1]
  final_grid <- grid_of_hyperparams
  if (grid_size > max_grid_size) {
    if (!is.nan(randomseed)) {set.seed(randomseed)}
    final_grid <-
      grid_of_hyperparams[sample(1:grid_size, max_grid_size, replace = FALSE),]
    grid_size <- max_grid_size
  }
  remove(grid_of_hyperparams); remove(default_grid)
  return(final_grid)
}

################## survsrf_cv ########################

#' Cross-validates SRF model
#'
#' @param df data frame with the data, "time" and "event" for survival outcome
#' @param predict.factors list of predictor names
#' @param fixed_time  at which performance metrics are computed
#' @param outer_cv k in k-fold CV, default 3
#' @param repeat_cv if NULL, runs once, otherwise repeats CV
#' @param randomseed random seed
#' @param return_models TRUE/FALSE, if TRUE returns all CV objects
#' @param inner_cv k in the inner loop of k-fold CV for SRF hyperparameters tuning, default is 3
#' @param tuningparams list of tuning parameters for random forest: 1) NULL for using a default tuning grid, or 2) a list("mtry"=c(...), "nodedepth" = c(...), "nodesize" = c(...))
#' @param parallel TRUE/FALSE if use parallel computations
#' @examples \donttest{
#' \dontshow{rfcores_old<- options()$rf.cores; options(rf.cores = 1)}
#' df <- simulate_nonlinear()
#' srf_cv <- survsrf_cv(df, names(df)[1:4])
#' summary(srf_cv)
#' \dontshow{options(rf.cores=rfcores_old)}
#' }
#' @return list of outputs
#' @export
survsrf_cv <- function(df,
                       predict.factors,
                       fixed_time = NaN,
                       outer_cv = 3,
                       inner_cv = 3,
                       repeat_cv = 2,
                       randomseed = NaN,
                       return_models = FALSE,
                       tuningparams = list(),
                       max_grid_size = 10,
                       parallel = FALSE) {
  
  Call <- match.call()
  inputs <- list(df,
                 predict.factors,
                 fixed_time,
                 outer_cv,
                 inner_cv,
                 repeat_cv,
                 randomseed,
                 return_models,
                 tuningparams,
                 max_grid_size,
                 parallel
                 )
  inputclass<- list(df = "data.frame",
                    predict.factors = "character",
                    fixed_time = "numeric",
                    outer_cv = "numeric",
                    inner_cv = "numeric", 
                    repeat_cv = "numeric",
                    randomseed = "numeric",
                    return_models = "logical",
                    tuningparams = "list",
                    max_grid_size = "list",
                    parallel = "logical")
  
  cp<- check_call(inputs, inputclass, Call)
  if (cp$anyerror) stop (paste(cp$msg[cp$msg!=""], sep=""))
  
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
    train_function = survsrf_train,
    predict_function = survsrf_predict,
    model_args = list(tuningparams = tuningparams,
                      max_grid_size= max_grid_size,
                      fixed_time = fixed_time,
                      randomseed = randomseed),
    model_name = "Survival Random Forest",
    parallel = parallel
  )
  output$call <- Call
  return(output)
}
