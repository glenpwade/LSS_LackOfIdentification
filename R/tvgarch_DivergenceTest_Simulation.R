
# The CRAN package "tvgarch" has been found to fail to converge when iteratively estimating a tvgarch model.
# This issue is especially worrisome when the package is used in Simulation studies, as the only indication 
# of a failure to estimate the model is the iteration count having reached its maximum.
# The onus is therefore on the user of the package to understand and check the iteration count property
# before choosing to include or exclude a series in a simulation.
# Better packages would set the models parameters to NA (or NaN) in failure situations to reduce the risk 
# for simulation studies.


# The simulation dataset is IID (Normal Noise), infused with the specific tvgarch process.
# 
# 1. Unsophisticated Users: May rely on the default starting parameters and optimiser controls provided by the package.
# 2. General Users: Will set starting parameters 'near' to the actual parameters, to support the optimizer not getting stuck on a local maximum.
# 3. Advanced Users: Are the optimizer controls exposed to users? parameter-scaling, finite-difference step size, relative tolerance?
#
# This simulation tests 2 common scenarios for simulations, assuming we have already generated a simulation dataset.
# We look at a) The unsophisticated case (no starting params provided) and b) General user (starting params provided)
    

library(tvgarch)
library(doParallel)

setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

# T2000:  ####

T5000_Data <- readRDS("T5000_Data.RDS")


# Test 1 Series: ####

if(FALSE){
    e = T2000_Data[,55]
    plot(e,type='l')
    
    # Default Start Pars
    tvg_mod_1 <- tvgarch::tvgarch(e,nr.trans)
    summary(tvg_mod_1)
    
    # Set Start Pars
    nr.trans <- 1
    initVals <- list()
    initVals$intercept.g = 0.4  # delta0
    initVals$size = 1.0         # delta1
    initVals$speed = 2.0        # throws an error for some reason:  Reason found: The optimser used ONLY supports eta.  (They convert gamma to eta in/out of the optimizer)
    initVals$location = 0.4
    #
    initVals$intercept.h = 0.13   # omega
    initVals$arch = 0.12          # alpha
    initVals$garch = 0.75         # beta
    
    tvg_mod <- tvgarch::tvgarch(e,nr.trans,initial.values = initVals, trace = TRUE)
    
    summary(tvg_mod)
    plot(tvg_mod)
    
}



# Run Simulation:  ####

# We want to store the parameters and their standard errors for both TV & GARCH
# We'll just use a matrix with TV$pars, GARCH$pars


# Setup the parallel backend
numCores <- 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

nr.trans <- 1

results_defStartPars = foreach(i=1:1100,.combine = rbind,.inorder = FALSE,.packages = "tvgarch")%dopar%{
    
    # Do first iteration of estimation:
    tvg_mod <- tvgarch::tvgarch(T2000_Data[,i],nr.trans,turbo = TRUE)
    # Return:
    c(tvg_mod$iter, coef(tvg_mod, "tvgarch"))
    
}

# Set start Pars
initVals <- list()
initVals$intercept.g = 0.4  # delta0
initVals$size = 1.0         # delta1
#initVals$speed = 5.0      # throws an error for some reason
initVals$location = 0.4
#
initVals$intercept.h = 0.1  # omega
initVals$arch = 0.1         # alpha
initVals$garch = 0.8         # beta

results_setStartPars = foreach(i=1:1100,.combine = rbind,.inorder = FALSE,.packages = "tvgarch")%dopar%{
    
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
    
    biasSet <- c(0.5,1.5,10,0.5,0.1,0.1,0.5) - colMeans(resultSet)
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

results <- results_setStartPars

# Diverged:
results[results[,1]==1000, ]
good_res <- results[results[,1]<1000, ]

##


results_InclBad <- results[1:1000, ]          # First 1000, include failed models
results_Valid <- results[results[,1]<1000, ]  # Filter out all models that diverged
results_Valid <- results_Valid[1:1000,]

sillyReport <- calcStats(results_InclBad)
stdReport <- calcStats(results_Valid)





