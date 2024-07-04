

################## deephit_predict ##################
# The function to predict event probability by a trained deephitmodel
#' @export
deephit_predict <-
  function(trained_model,
           newdata,
           fixed_time
  ){
    s1 <-
      predict(trained_model, newdata, type = "survival")
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
      if(is.null(deephitparams$num_nodes)) {deephitparams$num_nodes =  c(16, 16,16,16)}
      if(is.null(deephitparams$batch_size)) {
        deephitparams$batch_size =  min(max(64, round(dim(df_train)[1]/8, 0)),256)}
      if(is.null(deephitparams$epochs)) {deephitparams$epochs =  30}
      if(is.null(deephitparams$mod_alpha)) {deephitparams$mod_alpha =  0.2}
      if(is.null(deephitparams$sigma)) {deephitparams$sigma =  0.1}
      if(is.null(deephitparams$cuts)) {deephitparams$cuts =  10}

      grid_of_hyperparams <- expand.grid(
        dropout = deephitparams$dropout,
        learning_rate = deephitparams$learning_rate,
        num_nodes = deephitparams$num_nodes,
        batch_size = deephitparams$batch_size,
        epochs = deephitparams$epochs,
        mod_alpha = deephitparams$mod_alpha,
        sigma = deephitparams$sigma,
        cuts = deephitparams$cuts
      )

      #if only one option is given, we use it to train
      if (dim(grid_of_hyperparams)[1]==1){
      deephitm <- deephit(
        data = df_train,
        x = df_train[predict.factors],
        y = Surv(df_train$time, df_train$event),
        dropout = deephitparams$dropout,
        learning_rate = deephitparams$learning_rate,
        num_nodes = deephitparams$num_nodes,
        batch_size = deephitparams$batch_size,
        epochs = deephitparams$epochs,
        mod_alpha = deephitparams$mod_alpha,
        sigma = deephitparams$sigma,
        cuts = deephitparams$cuts
      )
      #if there are options, we tune
      }else{
        tuning_results <-
          deephit_tune(1,df_train,predict.factors,)


      }
    }
    return(deephitm)
  }

################## deephit_tune ##################
#'A repeated 3-fold CV over a hyperparameters grid
#'
#' @param repeat_tune number of repeats
#' @param df_tune data
#' @param predict.factors predictor names
#' @param fixed_time  not used here, but for some models the time for which performance is optimized
#' @param deephitparams list with the range of hyperparameters
#' @param inner_cv number of folds for each CV
#' @max_grid_size 50 by default. If grid size > max_grid_size, a random search is performed for max_grid_size iterations. Set this to a small number for random search
#' @return  output=list(cindex_ordered, bestparams)
#' @export
deephit_tune <-
  function(repeat_tune,
           df_tune,
           predict.factors,
           fixed_time = NaN,
           deephitparams = list(),
           inner_cv = 3,
           max_grid_size = 50
           ) {
    means = c()
    for (i in (1:repeat_tune)) {
      print (repeat_tune)
      ds_tune_fus <-
        deephit_tune_single(df_tune,
                            predict.factors,
                            fixed_time,
                            deephitparams,
                            inner_cv,
                            max_grid_size)
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
#' Internal function for deephit_tune(), performs 1 CV
#' @param df_tune data
#' @param predict.factors predictor names
#' @param fixed_time  not used here, but for some models the time for which performance is optimized
#' @param deephitparams list with the range of hyperparameters
#' @param inner_cv number of folds for each CV
#' @max_grid_size 1000 by default. If grid size > max_grid_size, a random search is performed for max_grid_size iterations. Set this to a small number for random search
#' #' @examples
#' \dontshow{rfcores_old <- options()$rf.cores; options(rf.cores=1)}
#' d <-simulate_nonlinear(100),
#' p<- names(d)[1:4]
#' deephitparams = list(
#'  "dropout" = c(0.1, 0.3),
#'  "learning_rate" = c(0.001),
#'  "num_nodes" =    list(c(8, 8), c(16, 16, 16, 4), c(32, 32, 32, 4)),
#'  "batch_size" = min(max(64, round(dim(df_train_cv)[1]/8, 0)),256),
#'  "epochs" = c(5, 50),
#'  "mod_alpha" = 0.2,
#'  "sigma" = 0.1,
#'  "cuts" = 10
#' )
#' deephit_tune_single(d, p,deephitparams = deephitparams, max_grid_size = 10)
#' \dontshow{options(rf.cores=rfcores_old)}
#' @return  output=list(grid, cindex, cindex_mean)
#' @export
deephit_tune_single <-
  function(df_tune,
           predict.factors,
           fixed_time = NaN,
           deephitparams=list(),
           inner_cv = 3,
           max_grid_size = 50) {

    default_grid <-
      list(
      "dropout" = c(0.1, 0.3),
      "learning_rate" = c(0.001),
      "num_nodes" =    list(c(8, 8), c(16, 16, 16, 4), c(32, 32, 32, 4)),
      "batch_size" = min(max(64, round(dim(df_tune)[1]/8, 0)),256),
      "epochs" = c(5, 50),
      "mod_alpha" = 0.2,
      "sigma" = 0.1,
      "cuts" = 10
      ) #size = 2x3x2= 12

    if (length(deephitparams) == 0) {
      grid_of_hyperparams <- default_grid
    }else{
      if (is.null(deephitparams$dropout)) {deephitparams$dropout = default_grid$dropout}
      if (is.null(deephitparams$learning_rate)) {deephitparams$learning_rate = default_grid$learning_rate}
      if (is.null(deephitparams$num_nodes)) {deephitparams$num_nodes = default_grid$num_nodes}
      if (is.null(deephitparams$batch_size)) {deephitparams$batch_size= default_grid$batch_size}
      if (is.null(deephitparams$epochs)) {deephitparams$epochs = default_grid$epochs}
      if (is.null(deephitparams$mod_alpha)) {deephitparams$mod_alpha = default_grid$mod_alpha}
      if (is.null(deephitparams$sigma)) {deephitparams$sigma = default_grid$sigma}
      if (is.null(deephitparams$cuts)) {deephitparams$cuts = default_grid$cuts}
    }

    grid_of_hyperparams <- expand.grid(
        "dropout" = deephitparams$dropout,
        "learning_rate" = deephitparams$learning_rate,
        "num_nodes" = deephitparams$num_nodes,
        "batch_size" = deephitparams$batch_size,
        "epochs" = deephitparams$epochs,
        "mod_alpha" = deephitparams$mod_alpha,
        "sigma" = deephitparams$sigma,
        "cuts" = deephitparams$cuts
    )

    print(grid_of_hyperparams)
    grid_size <- dim(grid_of_hyperparams)[1]

    # choose random sub-grid if grid_size too large
    if (grid_size > max_grid_size) {
      grid_of_values<- grid_of_hyperparams[sample(1:grid_size, max_grid_size,replace = FALSE), ]
      grid_size <- max_grid_size
    }else{
      grid_of_values <- grid_of_hyperparams
      }

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
          learning_rate = grid_of_values[i, "learning_rate"],
          num_nodes = grid_of_values[i, "num_nodes"][[1]],
          batch_size = grid_of_values[i, "batch_size"],
          epochs = grid_of_values[i, "epochs"],
          mod_alpha = grid_of_values[i, "mod_alpha"],
          sigma = grid_of_values[i, "sigma"],
          cuts = grid_of_values[i, "cuts"]
        )
        #check test performance
        pp <-
          deephit_predict(deephitm, df_test_cv, fixed_time)
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
    cindex_mean <- apply(cind, 1, mean)
    output = list()
    output$grid = grid_of_values
    output$cindex = cind
    output$cindex_mean = cindex_mean
    output$cindex_ordered <-
      cbind(grid_of_values, cindex_mean)[order(cindex_mean, decreasing = TRUE),]
    output$bestparams <- output$cindex_ordered[1,]
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
