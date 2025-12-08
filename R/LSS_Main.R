

library(tvgarch)
library(doParallel)
library(knitr)

#setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")
setwd("C:\\Repos\\LSS_LackOfIdentification")

# T2000:  ####

T5000_Data <- readRDS("T5000_Data.RDS")
T2000_Data <- Ref_Data[1:2000,]


# Test 1 Series: ####

if(FALSE){
    e = T2000_Data[,55]
    plot(e,type='l')
    

    nr.trans <- 1
    initVals <- list()
    initVals$intercept.g = 0.5  # delta0
    initVals$size = 1.5         # delta1
    #initVals$speed = 10.0      # throws an error for some reason
    initVals$location = 0.5
    #
    initVals$intercept.h = 0.1  # omega
    initVals$arch = 0.1         # alpha
    initVals$garch = 0.8         # beta
    
    tvg_mod <- tvgarch::tvgarch(e,nr.trans,initial.values = initVals, trace = TRUE)
    
    summary(tvg_mod)
    plot(tvg_mod)
    
   
}



# Run Simulation:  ####

# We want to store the parameters and their standard errors for both TV & GARCH
# We'll just use a matrix with TV$pars, GARCH$pars
Nr.Series <- 12000

# Setup the parallel backend
numCores <- 10
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

nr.trans <- 1

results_defStartPars = foreach(i=1:Nr.Series,.combine = rbind,.inorder = FALSE,.packages = "tvgarch")%dopar%{

    # Using package default starting pars
    tvg_mod <- tvgarch::tvgarch(T2000_Data[,i],nr.trans,turbo = TRUE)
    # Return:
    c(tvg_mod$iter, coef(tvg_mod, "tvgarch"))
    
}


# Stop the parallel cluster
stopCluster(cl)

saveRDS(results_setStartPars,"T2000_Results_lss_defStartPars.RDS")

# Set start Pars
initVals <- list()
initVals$intercept.g = 0.4  # delta0
initVals$size = 1.0         # delta1
initVals$speed = 2.0        # throws an error for values > 7 (suggests eta)
initVals$location = 0.4
#
initVals$intercept.h = 0.1  # omega
initVals$arch = 0.1         # alpha
initVals$garch = 0.8         # beta

results_setStartPars = foreach(i=1:Nr.Series,.combine = rbind,.inorder = FALSE,.packages = "tvgarch")%dopar%{
    
    # Do first iteration of estimation:
    tvg_mod <- tvgarch::tvgarch(T2000_Data[,i],nr.trans,initial.values = initVals,turbo = TRUE)
    # Return:
    c(tvg_mod$iter, coef(tvg_mod, "tvgarch"))
}

# Stop the parallel cluster
stopCluster(cl)

saveRDS(results_setStartPars,"T2000_Results_lss_setStartPars.RDS")

# Analyse the results:  ####

calcStats <- function(resultSet){
    
    resultSet <- resultSet[,c(2:8)]
    
    biasSet <- colMeans(resultSet) - c(0.5,1.5,logb(10),0.5,0.1,0.1,0.8)
    sdSet <- c(sd(resultSet[,1]),sd(resultSet[,2]),sd(resultSet[,3]),sd(resultSet[,4]),sd(resultSet[,5]),sd(resultSet[,6]),sd(resultSet[,7]))
    
    tblResultsGt <- matrix( c(biasSet[1],sdSet[1], biasSet[2],sdSet[2], biasSet[3],sdSet[3], biasSet[4],sdSet[4]), nrow = 1, ncol = 8 )
    colnames(tblResultsGt) <- c("d0","d0_se","d1","d1_se","spd","spd_se","loc","loc_se")
    rownames(tblResultsGt) <- c("meanBias, se: ")
    
    tblResultsHt <- matrix( c(biasSet[5],sdSet[5], biasSet[6],sdSet[6], biasSet[7],sdSet[7] ), nrow = 1, ncol = 6 )
    colnames(tblResultsHt) <- c("omega","omega_se","alpha","alpha_se","beta","beta_se")
    rownames(tblResultsHt) <- c("meanBias, se: ")
    
    print(round(tblResultsGt,4))
    
    print(round(tblResultsHt,4))
    
    #resTable <- table(tblResultsGt,dnn = dimnames(tblResultsGt))
    
}


# Analyse resuls
results <- results_setStartPars
results <- results_defStartPars

# Diverged:
divergentRes <- NROW(results[results[,1]==1000, ])
good_res <- results[results[,1]<1000, ]

##


results_InclBad <- results[1:12000, ]          # First 1000, include failed models
results_Valid <- results[results[,1]<1000, ]  # Filter out all models that diverged
results_Valid <- results_Valid[1:12000,]

sillyReport <- calcStats(results_InclBad)
stdReport <- calcStats(results_Valid)

# What are the default start pars?
tvg_mod <- tvgarch::tvgarch(T2000_Data[,1],nr.trans,turbo = TRUE,trace = TRUE)

# Initial parameters for the h component
 
#        intercept.h arch1 garch1
# Value:         0.1   0.1    0.7

# Initial parameters for the g component:
    
#     intercept.g size1 speed1 location1
# Value:        1   0.1 2.3026    0.5002


