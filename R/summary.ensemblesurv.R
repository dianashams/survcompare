
summary.ensemblesurv <-
  function(output_object) {
    if (!inherits(output_object, "list")) {stop("Not a legitimate object")}
    x <- output_object
    cat(
      "\nInternally validated test performance of Cox-PH, Cox-Lasso",
      "and Ensemble 1 (Cox-Survival Random Forest):\n"
    )
    # printing the table with mean test results, diff to Cox-PH
    print(round(x$results_mean, 6))
    cat("\n")
    pv <- x$difftest["pvalue", "C_score"]
    m <- x$difftest["m", "C_score"]
    if (x$difftest["pvalue", "C_score"] < 0.05) {
      cat(
        paste(
          "Survival Random Forest ensemble has outperformed Cox-PH by ",
          round(m, 4),
          " in C-index.\n",
          "The difference is statistically significant with p-value = ",
          ifelse(pv < 0.001, round(pv, 6), round(pv, 4)),
          "*", ifelse(pv < 0.01, "*", ""), ifelse(pv < 0.001, "*", ""),".\n",
          "The supplied data may contain non-linear or cross-term dependencies",
          " better captured by survival random forest.\n",
          sep = ""
        )
      )
    } else {
      cat(
        paste(
          "Survival Random Forest ensemble has NOT outperformed ",
          "Cox-PH model.\nThe difference in validated C-index is ",
          round(m, 4),
          " , which is not statistically significant with p-value = ",
          round(pv, 4),".\n",sep=""
        )
      )
    }

    ms<- x$main_stats
    f<- function(i){
      paste(round(mean(ms[i, "mean"]),4),
            "(95CI=", round(mean(ms[i,"95CILow"]),4),"-", round(mean(ms[i,"95CIHigh"]),4),
            ";SD=", round(mean(ms[i,"sd"]),4), ")", sep="")
    }
    cat(
      paste(
        "C-score: \n      CoxPH  ", f(1),"\n  Ensemble1 ", f(2),
        "\nAUCROC:\n      CoxPH  ", f(3),"\n  Ensemble1 ", f(4), sep=""
      )
    )
  }
