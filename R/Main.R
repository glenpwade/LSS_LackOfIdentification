#
# Journal of Time Series Analysis:
# Robust Estimation and Inference for Time-Varying Unconditional Volatility
#
# By, Adam Lee | Rickard Sandberg | Genaro Sucarrat
#
#
#  This paper discredits the MTVGARCH model in a simulation comparison 
#  However the comparison is flawed as the MTVGARCH model has not be correctly implemented.
# 
#  The code below will correctly implement the MTVGARCH model in a similar simulation to provide a fair comparison.


# The proposed simulation is as follows:
#
# 1. Create the TVGARCH Model Spec with the provided parameters
# 2. Generate Simulated data sets injected with this process, 1000 Replications, T = 2000
# 3. Use MTVGARCH modelling processes to identify & estimate the simulated series
# 4. Report the estimation results

# Install required Packages - once
if(FALSE){
    install.packages("MTVGARCH")  # Install using R-Studio, using the install from file option - MTVGARCH_0.9.5.7.tar.gz provided
    install.packages("knitr")  
}

# Unload the tvgarch package to avoid masking issues
try(detach("package:tvgarch", unload = TRUE, character.only = TRUE))
# Load required packages
library(MTVGARCH)   # Ver. 0.9.5.7
library(knitr)

# Set a working directory
setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")
#setwd("C:\\Repos\\LSS_LackOfIdentification")

# Gen Data ####

if(FALSE){
    
    Tobs = 5000
    st = (1:Tobs)/Tobs
    shape = tvshape$single
    # Create the tv object with default parameters
    TV <- MTVGARCH::tv(st,shape)
    
    # Now override parameters with desired model specification
    TV$speedopt <- speedopt$gamma
    TV$delta0 = 0.5
    TV$pars["deltaN",1] = 1.5
    TV$pars["speedN",1] = 10
    TV$pars["locN1",1] = 0.5
    
    # Create the garch object with default parameters
    GARCH = MTVGARCH::garch(garchtype$general)
    # Now override parameters with desired model specification
    GARCH$pars["omega",1] = 0.1
    GARCH$pars["alpha",1] = 0.1
    GARCH$pars["beta",1] = 0.8
    
    
    
    # Generate 1000 series of simulated data, injected with this process:
    ## noiseDist is a named-list describing the error-distribution and parameters
    noiseDist <- list()
    noiseDist$name = 'Normal'     
    noiseDist$mean = 0  
    noiseDist$sd = 1
    # set.seed(1984)
    # Ref_Data <- MTVGARCH::generateRefData(nr.series = 11000,nr.obs = Tobs,tvObj = TV,garchObj = GARCH)
    # saveRDS(Ref_Data,"T5000_Data.RDS")
    # saveRDS(Ref_Data[1:2000, ],"T2000_Data.RDS")
    
    
    
}

# T2000:  
T2000_Data <- readRDS("T2000_Data.RDS")

# Test 1 Series: ####

if(FALSE){
    
    e = T2000_Data[,1]
    plot(e,type='l')
    
    Tobs = 2000
    st = (1:Tobs)/Tobs
    shape = tvshape$single
    # Create the TV Specification and set starting params
    TV <- MTVGARCH::tv(st,shape)
    TV$delta0 = 1.0
    TV$pars["deltaN",1] = 0.1
    TV$pars["speedN",1] = logb(10)
    TV$pars["locN1",1] = 0.5002
    TV$optimcontrol$parscale <- c(1,3,5,1)
    
    TV <- estimateTV(e,TV)
    summary(TV)
    plot(TV)
    
    # Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    
    estCtrl <- list(calcSE = TRUE, verbose = TRUE)  # Since we are only doing one - let's see what's going on & calc parameter se's.
    # Run the iterative estimation
    mod <- MTVGARCH::estimateTVGARCH(e,mod,estCtrl,autoConverge = TRUE)

    # Check convergence & iterations:
    mod$Estimated$converged
    mod@iterations
    
    #View the full estimated object:
    View(mod)
    
}


# Run Simulation:  ####

# Setup the parallel backend
numCores <- 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

# Create the TV Specification:
Tobs = 2000
st = (1:Tobs)/Tobs
shape = tvshape$single

# Set the estimation controls to suppress console output
estCtrl <- list(calcSE = FALSE, verbose = FALSE)

# Save the results for reporting:  (Col:1 'EstimationMethod': Two-Step=1, Iterative=2), (Col10: 'EstimationError': 0(FALSE) / 1(TRUE))
pars <- matrix(NA,1,10)
colnames(pars) <- c("EstMethod","NrIterations","d0","d1","spd","loc","omega","alpha","beta","EstError")

Nr.Series <- 11    # Simulate 1100 estimations, to ensure we have 1000 that complete without error

# 2-STEP: ####
cat("\nTWO-STEP:\n")
timestamp()
results_2S = foreach(i=1:Nr.Series,.combine = rbind,.inorder = TRUE,.packages = "MTVGARCH")%dopar%{

    e = T2000_Data[,i]

    # Do Two-Step estimation:
    # Create the TV Specification and set starting params = tvgarch defaults (Ensure a fair comparison)
    TV <- MTVGARCH::tv(st,shape)
    TV$delta0 = 1.0
    TV$pars["deltaN",1] = 0.1
    TV$pars["speedN",1] = logb(10)
    TV$pars["locN1",1] = 0.5002
    TV$optimcontrol$parscale <- c(1,2,5,1)
    
    # Step1: Estimate TV
    TV <- estimateTV(e,TV,estCtrl)
    # Step2: Estimate Garch on filtered data
    GARCH <- garch(garchtype$general)
    GARCH <- estimateGARCH(e,GARCH,estCtrl,TV)    # Will estimate the h(t) after filtering g(t) from e
    pars[1,] <- c(1,2,TV$Estimated$delta0,TV$Estimated$pars[1:3],GARCH$Estimated$pars,as.numeric(TV$Estimated$error) )
    # Note: Any failed estimations will be identified in the last column, so the unestimated parameters can be excluded

    # Return:
    pars
    
}
timestamp()

# ITERATIVE: ####
cat("\nITERATIVE:\n")
timestamp()
results_Iter = foreach(i=1:Nr.Series,.combine = rbind,.inorder = TRUE,.packages = "MTVGARCH")%dopar%{
    
    e = T2000_Data[,i]

    # Create the TV Specification and set starting params
    TV <- MTVGARCH::tv(st,shape)
    TV$delta0 = 1.0
    TV$pars["deltaN",1] = 0.1
    TV$pars["speedN",1] = logb(10)
    TV$pars["locN1",1] = 0.5002
    TV$optimcontrol$parscale <- c(1,2,5,1)
    
    # Do initial Estimation the g(t) assuming h(t) = 1
    TV <- estimateTV(e,TV,estCtrl)
    
    # Specify a multiplicitive tvgarch object using this TV obj and a general GARCH specification
    TVG <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    # Do Iterative estimation
    TVG <- estimateTVGARCH(e,TVG,estCtrl,autoConverge = TRUE)
    
    # Extract the estimated parameters:
    tvpars <- TVG$Estimated$tv
    garchpars <- TVG$Estimated$garch
    pars[1,] <- c(2,TVG@iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,as.numeric(!(TVG$Estimated$converged)) )
    
    # Return:
    pars
    
}
timestamp()


# Stop the parallel cluster
stopCluster(cl)

saveRDS(rbind(results_2S,results_Iter),"T2000_Results_STW.RDS")  #STW: Silvennoinen, Terasvirta, Wade


# Analyse the Results:  ####

calcStats <- function(resultSet){
    
    resultSet <- resultSet[,c(3:9)]  #Extract the parameters
    
    biasSet <- c(0.5,1.5,logb(10),0.5,0.1,0.1,0.8) - colMeans(resultSet)
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


results <- readRDS("T2000_Results_STW.RDS")

# Any Failures?
NrFails <- NROW(results[results[,10]==0, ])

# Any high iteration/slow converging?
SlowConverges <- results[results[,2] > 10, ]
#
# Conclusion: Slow but Accurate in the end!

results_2S <- results[results[,1]==1, ]  # Col#1 = 1, 2-Step
results_Iter <- results[results[,1]==2, ]  # Col#1 = 2, Iterative

summary(results_Iter[1:1000,(3:9)])  # V1-7: d0,d1,spd,loc, omega,alpha,beta
summary(results_Iter[1:1000,2])  # Col2: Iteration Count

calcStats(results_2S[1:1000,])
calcStats(results_Iter[1:1000,])


