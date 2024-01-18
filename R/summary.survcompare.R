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
    if (!inherits(object, "survcompare")) {
      stop("Not a \"survcompare\" object")
    }

    # Cox model Lasso or PH
    coxend <- ifelse(object$useCoxLasso, "CoxLasso  ", "CoxPH    ")
    x <- object
    cat(
      "\nInternally validated test performance of",
      coxend,
      "and Survival Random Forest ensemble:\n"
    )
    # printing the table with mean test results, diff to CoxPH
    print(round(x$results_mean, 4))
    cat("\n")
    pv <- x$difftest["pvalue", "C_score"]
    m <- x$difftest["m", "C_score"]
    msd<- x$difftest["std", "C_score"]
    pvstars <-
      ifelse(pv < 0.001, "***", ifelse(pv <= 0.01, "**", ifelse(pv <= 0.05, "*", "")))

    # compile the output message:
    if (x$difftest["pvalue", "C_score"] < 0.05) {
      t1 <- paste(
        "Survival Random Forest ensemble has outperformed ",
        coxend,
        "by ",
        round(m, 4),
        " in C-index.\n",
        sep = ""
      )
      t2 <-
        "The difference is statistically significant with the p-value "
      t3 <-
        ".\nThe supplied data may contain non-linear or cross-term dependencies, \nbetter captured by the Survival Random Forest.\n"
    } else{
      t1 <- paste(
        "Survival Random Forest ensemble has NOT outperformed ",
        coxend,
        "with mean c-index difference of",
        round(m, 4),
        ".\n",
        sep = ""
      )
      t2 <-
        "The difference is not statistically significant with the p-value = "
      t3 <-
        ". \nThe data may NOT contain considerable non-linear or cross-term dependencies, \nthat could be captured by the Survival Random Forest.\n"
    }

    # print the output message
    cat(paste(t1, t2, ifelse(pv < 0.001, round(pv, 8), round(pv, 4)), pvstars, t3, sep = ""))

    #print the main stats
    ms <- x$main_stats
    f <- function(i) {
      paste(
        round(mean(ms[i, "mean"]), 4),
        "(95CI=",
        round(mean(ms[i, "95CILow"]), 4),
        "-",
        round(mean(ms[i, "95CIHigh"]), 4),
        ";SD=",
        round(mean(ms[i, "sd"]), 4),
        ")",
        sep = ""
      )
    }
    cat(
      paste(
        "C-score: \n  ",
        coxend,
        "  ",
        f(1),
        "\n  SRF_Ensemble ",
        f(2),
        "\nAUCROC:\n  ",
        coxend,
        "  ",
        f(3),
        "\n  SRF_Ensemble ",
        f(4),
        sep = ""
      )
    )
    invisible(object)
  }

#' Print survcompare object
#'
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
