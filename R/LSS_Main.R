

library(tvgarch)

setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

# T2000:  ####

T2000_Data <- readRDS("T2000_Data.RDS")


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


# Setup the parallel backend
numCores <- 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

nr.trans <- 1
initVals <- list()
initVals$intercept.g = 0.5  # delta0

results = foreach(i=1:1100,.combine = rbind,.inorder = FALSE,.packages = "tvgarch")%dopar%{
    

    # Do first iteration of estimation:
    tvg_mod <- tvgarch::tvgarch(T2000_Data[,i],nr.trans,initial.values = initVals,turbo = TRUE)

    # Return:
    c(tvg_mod$iter, coef(tvg_mod, "tvgarch"))
    
}

# Stop the parallel cluster
stopCluster(cl)

saveRDS(results,"T2000_Results_lss.RDS")

# Analyse the results:  ####

# Diverged:
results[results[,1]==1000, ]

good_res <- results[results[,1]<1000, ]
mean(good_res[,3])

sd(results[,3])

sd(good_res[,3])





