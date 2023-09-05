
summary.survcompare <-
  function(output_object, useCoxLasso) {
    if (!inherits(output_object, "list")) {stop("Not a legitimate object")}
    coxend <- ifelse(useCoxLasso, "CoxLasso  ", "CoxPH    ")
    x <- output_object
    cat(
      "\nInternally validated test performance of",coxend,"and Survival Random Forest ensemble:\n"
    )
    # printing the table with mean test results, diff to CoxPH
    print(round(x$results_mean, 4))
    cat("\n")
    pv <- x$difftest["pvalue", "C_score"]
    m <- x$difftest["m", "C_score"]
    pvstars <-  ifelse(pv < 0.001, "***", ifelse(pv <= 0.01, "**", ifelse(pv <= 0.05, "*", "")))
    
    #output message:
    if (x$difftest["pvalue", "C_score"] < 0.05) {
      t1<- paste("Survival Random Forest ensemble has outperformed ",
                 coxend, "by ", round(m, 4),pvstars, " in C-index.\n",sep="")
      t2<- "The difference is statistically significant, p-value = "
      t3<- ".\nThe supplied data may contain non-linear or cross-term dependencies, \nbetter captured by the Survival Random Forest.\n"
    }else{
      t1<- paste("Survival Random Forest ensemble has NOT outperformed ",
                 coxend, "with mean c-index difference of", round(m, 4),pvstars, ".\n",sep="")
      t2<- "The difference is not statistically significant, p-value = "
      t3<- ". \nThe data may NOT contain considerable non-linear or cross-term dependencies, \nthat could be captured by the Survival Random Forest.\n"
    }
    # print the resulting message
    cat(paste(t1, t2, ifelse(pv < 0.001, round(pv, 6), round(pv, 4)), t3, sep = ""))
    
    #print main stats
    ms<- x$main_stats
    f<- function(i){
      paste(round(mean(ms[i, "mean"]),4),
            "(95CI=", round(mean(ms[i,"95CILow"]),4),"-", round(mean(ms[i,"95CIHigh"]),4),
            ";SD=", round(mean(ms[i,"sd"]),4), ")", sep="")
    }
    cat(
      paste(
        "C-score: \n  ",coxend,"  ", f(1),"\n  SRF_Ensemble ", f(2),
        "\nAUCROC:\n  ",coxend,"  ", f(3),"\n  SRF_Ensemble ", f(4), sep=""
      )
    )
  }
