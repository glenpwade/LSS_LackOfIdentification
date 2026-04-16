calcStats <- function(resultSet,TVpars,Garchpars) {
    
    # Debug:
    if(FALSE){
        resultSet = results2S
        resultSet = resultsIter
        TVpars = TVpars
        Garchpars = GARCHpars
    }
    #
    
    resultSet <- resultSet[,c(2:8)]  #Extract the parameters

    biasSet <- colMeans(resultSet,na.rm = TRUE) - c(TVpars,Garchpars)
    sdSet <- c(sd(resultSet[,1],na.rm = TRUE),sd(resultSet[,2],na.rm = TRUE),sd(resultSet[,3],na.rm = TRUE),sd(resultSet[,4],na.rm = TRUE)
               ,sd(resultSet[,5],na.rm = TRUE),sd(resultSet[,6],na.rm = TRUE),sd(resultSet[,7],na.rm = TRUE))
    
    tblResultsGt <- matrix( c(biasSet[1],sdSet[1], biasSet[2],sdSet[2], biasSet[3],sdSet[3], biasSet[4],sdSet[4]), nrow = 1, ncol = 8 )
    colnames(tblResultsGt) <- c("d0","d0_se","d1","d1_se","spd","spd_se","loc","loc_se")
    rownames(tblResultsGt) <- c("meanBias, se: ")
    
    tblResultsHt <- matrix( c(biasSet[5],sdSet[5], biasSet[6],sdSet[6], biasSet[7],sdSet[7] ), nrow = 1, ncol = 6 )
    colnames(tblResultsHt) <- c("omega","omega_se","alpha","alpha_se","beta","beta_se")
    rownames(tblResultsHt) <- c("meanBias, se: ")
    
    resTableG <- kable(round(tblResultsGt,getOption("digits")),caption="g(t)")
    resTableH <- kable(round(tblResultsHt,getOption("digits")),caption="h(t)")
    
    print(resTableG)
    print(resTableH)
    
}

