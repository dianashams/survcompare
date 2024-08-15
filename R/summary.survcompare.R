#' Summary of survcompare results
#'
#' @param object output object of the survcompare function
#' @param ... additional arguments to be passed
#' @return object
#' @export
summary.survcompare <-
  function(object, ...) {
    # prints summary statement for the output of the 'survcompare' function

    #check
    if (!inherits(object, "survcompare")) {stop("Not a \"survcompare\" object")}

    # Cox model Lasso or PH
    coxend <- object$model_name_base

    # e.g. SRF DeepHit DeepSurv
    mlmodel = object$model_name

    # Summary title to print
    x <- object
    cat(
      "\nInternally validated test performance of",
      coxend,       "and" ,       mlmodel,
      "over", object$cv[3], "repeated",
      object$cv[1], "fold cross-validations (inner k =",
      object$cv[2], ").",
      "Mean performance:\n")

    # PRINT MEAN TEST RESULTS
    printing_columns= c("T", "C_score", "AUCROC","Calib_slope","sec")
    xx = x$results_mean[printing_columns]
    print(round(x$results_mean[printing_columns], 4))

    cat("\nMedian performance:\n")
    # MEDIAN RESULTS
    printing_columns= c("T", "C_score", "AUCROC","Calib_slope","sec")
    print(round(x$results_median[printing_columns], 4))

    cat("\n")
    pv <- x$difftest["pvalue", "C_score"]
    m <- x$difftest["m", "C_score"]
    msd<- x$difftest["std", "C_score"]
    pvstars <-
      ifelse(pv < 0.001, "***", ifelse(pv <= 0.01, "**", ifelse(pv <= 0.05, "*", "")))

    # compile the output message:
    if (x$difftest["pvalue", "C_score"] < 0.05) {
      t1 <- paste(
        mlmodel, " has outperformed ",
        coxend, "by ", round(m, 4)," in C-index.\n",sep = ""    )
      t2 <-
        "The difference is statistically significant with the p-value "
      t3 <-
        paste(".\nThe supplied data may contain non-linear or cross-term dependencies, \nbetter captured by ",
              mlmodel, ".\n", sep="")
    } else{
      t1 <- paste(
        mlmodel, " has NOT outperformed ",
        coxend, " with the mean c-index difference of ",
        round(m, 4), ".\n",sep = "")
      t2 <-
        "The difference is not statistically significant with the p-value = "
      t3 <-
        paste(". \nThe data may NOT contain considerable non-linear or cross-term dependencies\nthat could be captured by ",
              mlmodel, ".\n", sep="")
    }
    # print the output message
    cat(paste(t1, t2, ifelse(pv < 0.001, round(pv, 8), round(pv, 4)), pvstars, t3, sep = ""))

    #print the main stats
    ms <- x$main_stats_pooled
    f <- function(i) {
      paste(round(mean(ms[i, "mean"]), 4),
            "(95CI=",round(mean(ms[i, "95CILow"]), 4),"-",
            round(mean(ms[i, "95CIHigh"]), 4),";SD=",
            round(mean(ms[i, "sd"]), 4),")",sep = "")}
    #' to check which row is for AUCROC for Cox and ML,
    #' if baseline ML was also reported, these are rows 4 and 5
    #' if not, then rows 3 and 4, so we add 1 if needed to 3 and 4
    ii = ifelse(dim(ms)[1]==4, 0, 1)
    cat(
      paste( "Mean C-score: \n  ",coxend,"  ",
             f(1),"\n  ", mlmodel, " ", f(2),
             "\nMean AUCROC:\n  ",coxend,"  ",f(3+ii),
             "\n",mlmodel, " ", f(4+ii),sep = ""))
    invisible(object)
  }

#' Print survcompare object
#' @param x output object of the survcompare function
#' @param ... additional arguments to be passed
#' @return x
#' @export
print.survcompare <- function(x, ...) {
  #check
  if (!inherits(x, "survcompare")) {
    stop("Not a \"survcompare\" object")
  }
  summary.survcompare(x)
  cat("\n", "See other items as x$item. Items available:\n")
  cat(names(x), sep = ", ")
  invisible(x)
}
