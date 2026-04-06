calcStats <- function(resultSet,TVpars,Garchpars) {
    
    # Debug:
    if(FALSE){
        resultSet = results_2S
        resultSet = results_Iter
        TVpars = TVparScale
        Garchpars = GARCHparScale
    }
    #
    
    resultSet <- resultSet[,c(2:8)]  #Extract the parameters
    
    biasSet <- colMeans(resultSet) - c(TVpars,Garchpars)
    sdSet <- c(sd(resultSet[,1]),sd(resultSet[,2]),sd(resultSet[,3]),sd(resultSet[,4]),sd(resultSet[,5]),sd(resultSet[,6]),sd(resultSet[,7]))
    
    tblResultsGt <- matrix( c(biasSet[1],sdSet[1], biasSet[2],sdSet[2], biasSet[3],sdSet[3], biasSet[4],sdSet[4]), nrow = 1, ncol = 8 )
    colnames(tblResultsGt) <- c("d0","d0_se","d1","d1_se","spd","spd_se","loc","loc_se")
    rownames(tblResultsGt) <- c("meanBias, se: ")
    
    tblResultsHt <- matrix( c(biasSet[5],sdSet[5], biasSet[6],sdSet[6], biasSet[7],sdSet[7] ), nrow = 1, ncol = 6 )
    colnames(tblResultsHt) <- c("omega","omega_se","alpha","alpha_se","beta","beta_se")
    rownames(tblResultsHt) <- c("meanBias, se: ")
    
    resTableG <- kable(round(tblResultsGt,4),caption="g(t)")
    resTableH <- kable(round(tblResultsHt,4),caption="h(t)")
    
    print(resTableG)
    print(resTableH)
    
}