# tvgarch package:
# Initial parameters for the h component
# 
#        intercept.h arch1 garch1
# Value:         0.1   0.1    0.7
# 
# Initial parameters for the g component:
#     
#        intercept.g size1 speed1 location1
# Value:           1   0.1 2.3026    0.5002
#
# Note: speed default indicates that SpeedOption = eta


# Install required Packages - once
if(FALSE){
    install.packages("tvgarch")
    install.packages("doParallel")
    install.packages("knitr")  
}

# Unload the MTVGARCH package to avoid masking issues
try(detach("package:MTVGARCH", unload = TRUE, character.only = TRUE))
# Load required packages
library(tvgarch)
library(doParallel)
library(knitr)

# Set a working directory
setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")


# T2000:  ####
T2000_Data <- readRDS("T2000_Data_Seed1984.RDS")


# Test 1 Series: ####

if(FALSE){
    e = T2000_Data[,1]
    plot(e,type='l')

    tvg_mod <- tvgarch::tvgarch(e,nr.trans = 1, trace = TRUE)
    
    summary(tvg_mod)
    plot(tvg_mod)
    
}


# Run Simulation:  ####

Nr.Series <- 1100    # Simulate 1100 estimations, to ensure we have 1000 that complete without error

# Setup the parallel backend
numCores <- 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

nr.trans <- 1

results_defStartPars = foreach(i=1:Nr.Series,.combine = rbind,.inorder = TRUE,.packages = "tvgarch")%dopar%{

    # Using package default starting pars
    tvg_mod <- tvgarch::tvgarch(T2000_Data[,i],nr.trans,turbo = TRUE)
    # Return:
    c(tvg_mod$iter, coef(tvg_mod, "tvgarch"))
    
}


# Stop the parallel cluster
stopCluster(cl)

saveRDS(results_defStartPars,"T2000_Results_lss_defStartPars.RDS")

# Experiment with setting start pars close to simulated process:
if(FALSE){
    # Set start Pars
    initVals <- list()
    initVals$intercept.g = 0.4  # delta0
    initVals$size = 1.0         # delta1
    initVals$speed = 2.0        # throws an error for values > 7 (suggests eta)
    initVals$location = 0.4
    #
    initVals$intercept.h = 0.1  # omega
    initVals$arch = 0.1         # alpha
    initVals$garch = 0.8        # beta
    
    results_setStartPars = foreach(i=1:Nr.Series,.combine = rbind,.inorder = FALSE,.packages = "tvgarch")%dopar%{
        
        # Do first iteration of estimation:
        tvg_mod <- tvgarch::tvgarch(T2000_Data[,i],nr.trans,initial.values = initVals,turbo = TRUE)
        # Return:
        c(tvg_mod$iter, coef(tvg_mod, "tvgarch"))
    }
    saveRDS(results_setStartPars,"T2000_Results_lss_setStartPars.RDS")
    
#  Note:  This experiment did not produce significantly different results from the default parameter simulation
#         While we know that starting parameters can have a significant affect on estimation, we conclude in this case
#         that the default params are sufficiently close to the simulated process to produce reasonable results
#
}



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
    
    resTableG <- kable(round(tblResultsGt,4),caption="g(t)")
    resTableH <- kable(round(tblResultsHt,4),caption="h(t)")
    
    print(resTableG)
    print(resTableH)
    
}


# Analyse results

results <- readRDS("Results/T2000_Results_lss_defStartPars.RDS")

summary(results[,"speed1"])
sd(results[,"speed1"])

# Diverged Summary:
divergentRes <- results[results[,1]==1000, ]
summary(divergentRes)

# Diverged Count:
divergentNr <- NROW(results[results[,1]==1000, ])

# Identify the series index that diverged
divergentIDs <- results[results[,1]==1000, ]
tmp <- row.names(divergentIDs[])
tmp <- substr(a,8,length(tmp))
divergentIDs <- as.integer(tmp)
saveRDS(divergentIDs,"T2000_Results_lss_divergentIDs.RDS")

results_InclBad <- results[1:1000, ]          # First 1000, include failed models
results_Valid <- results[results[,1]<1000, ]  # Filter out all models that diverged
results_Valid <- results_Valid[1:1000,]       # Take First 1000

sillyReport <- calcStats(results_InclBad)
stdReport <- calcStats(results_Valid)



