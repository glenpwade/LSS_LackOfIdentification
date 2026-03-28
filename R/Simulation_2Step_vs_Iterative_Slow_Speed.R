if(FALSE)
{
#
# Journal of Time Series Analysis:
# Robust Estimation and Inference for Time-Varying Unconditional Volatility
#
# By, Adam Lee | Rickard Sandberg | Genaro Sucarrat
#
#
#  This paper discredits the MTVGARCH model in a simulation comparison and erroneously claims that 2-step estimation out-performs iterative.
#  However the comparison is flawed as the MTVGARCH model has not be correctly implemented and the variance level-change used is very low.
# 
#  The code below will correctly implement the MTVGARCH model in a similar simulation to provide a fair comparison between 2-Step & Iterative.



# The method is as follows:
#
# 1. Create a TVGARCH Model Spec with the parameters used in the paper. Note: The variance level shift is low (delta0 = 0.5, delta1=1.5) and the Persistence is low (GarchPars: 0.1,0.1,0.8)
# 2. Create a TVGARCH Model Spec with the parameters similar to those used in the paper, but with a higher variance level shift.  (delta0 = 0.5, delta1=4.0)
# 3. Compare 2-Step & Iterative performance when the transition location is close to a boundary. (loc = 0.8)
# 4. Compare 2-Step & Iterative performance when the transition speed is very slow. (eta = 1.5)
# 5. Compare 2-Step & Iterative performance when the transition speed is very fast. (eta = 5.5)

# Simulated Dataset Details:
#
#  Garch_Pars: 0.05, 0.05, 0.90
#
# 1. T2000_Low_VarShift        TV_Pars: d0=0.5, d1=1.5, spd=2.3, loc=0.5  
# 2. T2000_Med_VarShift     TV_Pars: d0=0.5, d1=4.0, spd=2.3, loc=0.5  
# 3. T2000_Near_Boundary    TV_Pars: d0=0.5, d1=4.0, spd=2.3, loc=0.8  
# 4. T2000_Slow_Speed       TV_Pars: d0=0.5, d1=4.0, spd=1.5, loc=0.5  
# 5. T2000_Fast_Speed       TV_Pars: d0=0.5, d1=4.0, spd=5.5, loc=0.5  
    

# Preparation:
#
# 1. Generate & save simulated data sets, with 3000 Replications and T = 2000
# 2. Each Simulated data set will be injected with the tv & garch process. 
# 3. Use MTVGARCH modelling processes to identify & estimate the simulated series
# 4. Report the estimation results
# 5. Note: 
#           The starting parameters will be set = process parameters.
#           The default optim() control parameters will be used, with the exception of parscale which is set according to the documented recommendation.
#           We will use 'eta' as the speed option.
}

# Initialise ####
library(MTVGARCH)   # Ver. 0.9.6.1
library(foreach)
library(doParallel)

# Set a working directory
setwd("C:\\Repos\\LSS_LackOfIdentification")

# Set constants
Reps <- 3000
Tobs <- 2000
parScale <- c(0.5,4.0,1.5,0.5)

# START SIMS HERE: ####

# Setup the parallel backend ####
#numCores <- parallel::detectCores() - 4
numCores <- 2
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

# Run the Simulation:  ####

# Set the iteration count 
#Reps <- 50    # Simulate 3000 estimations - set lower for debugging any parallel issues
# Set the estimation controls to suppress console output
estCtrl <- list(calcSE = FALSE, verbose = FALSE)

# Get the data:
fileName = "T2000_Slow_Speed"
filePath = paste0(".\\SimSourceData\\",fileName,".RDS")
simData <- readRDS(filePath)

# Create the TV Specification - to be estimated later:
Tobs = NROW(simData)
st = (1:Tobs)/Tobs
shape = tvshape$single
# Create the TV Specification and set starting params to match the loaded Dataset
TVspec <- MTVGARCH::tv(st,shape)
TVspec$delta0 = 0.5
TVspec$pars["deltaN",1] = 4.0
TVspec$pars["speedN",1] = 1.5
TVspec$pars["locN1",1] = 0.5
TVspec$optimcontrol$parscale <- parScale

# Save the results for reporting:  (Col:1 'EstimationMethod': Two-Step=1, Iterative=2), (Col10: 'EstimationError': 0(FALSE) / 1(TRUE))
pars <- matrix(NA,1,10)
colnames(pars) <- c("EstMethod","NrIterations","d0","d1","spd","loc","omega","alpha","beta","EstError")

# 2-STEP: ####
cat("\nTWO-STEP:\n")
timestamp()
results_2S = foreach(i=1:Reps, .combine = rbind, .inorder = TRUE, .packages = "MTVGARCH")%dopar%{

    e = simData[,i]

    # 1. Do initial Estimation of g(t) assuming h(t) = 1
    TV <- estimateTV(e,TVspec,estCtrl)
    TV$optimcontrol$parscale <- parScale
    
    # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    # 2.1 We need to set the Garch starting pars, before estimating the model:
    mod$garchpars["omega",1] = 0.05 
    mod$garchpars["alpha",1] = 0.05 
    mod$garchpars["beta",1] = 0.90
    mod$garchOptimcontrol$parscale <- c(0.05,0.05,0.90)
    
    # 3. Run the 2-Step estimation
    mod$tvOptimcontrol$reltol <- 1e-5    #Hack: This value is used as Threshold for the Iteration Convergence
    mod_2s <- MTVGARCH::estimateTVGARCH_2Step(e,mod,estCtrl)
    
    # 4. Extract the estimated parameters:
    tvpars <- mod_2s$Estimated$tv
    garchpars <- mod_2s$Estimated$garch
    pars[1,] <- c(2,mod_2s@iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,as.numeric(!(mod_2s$Estimated$converged)) )
    
    # Note: Any failed estimations will be identified in the last column, so the unestimated parameters can be excluded

    # Return:
    pars
    
}
timestamp()  # ~1 min

# ITERATIVE: ####
cat("\nITERATIVE:\n")
timestamp()
results_Iter = foreach(i=1:Reps, .combine = rbind, .inorder = TRUE, .packages = "MTVGARCH")%dopar%{
    
    e = simData[,i]

    # Do initial Estimation of g(t) assuming h(t) = 1
    TV <- estimateTV(e,TVspec,estCtrl)
    TV$optimcontrol$parscale <- parScale
    
    # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    # 2.1 We need to set the Garch starting pars, before estimating the model:
    mod$garchpars["omega",1] = 0.05 
    mod$garchpars["alpha",1] = 0.05 
    mod$garchpars["beta",1] = 0.90
    mod$garchOptimcontrol$parscale <- c(0.05,0.05,0.90)
    
    # 3. Run the Iterative estimation
    mod$tvOptimcontrol$reltol <- 1e-5    #Hack: This value is used as Threshold for the Iteration Convergence
    mod_iter <- estimateTVGARCH(e,mod,estCtrl,autoConverge = TRUE)
    
    # 4. Extract the estimated parameters:
    tvpars <- mod_iter$Estimated$tv
    garchpars <- mod_iter$Estimated$garch
    pars[1,] <- c(2,mod_iter@iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,as.numeric(!(mod_iter$Estimated$converged)) )
    
    # Return:
    pars
    
}
timestamp()  # ~? mins

# Stop the parallel cluster & remove 'cl'  
# Note: Best to not run this when executing in Parallel Mode / Background Job.  Wait until all tasks/jobs complete, then tidy Up
# stopCluster(cl)
# rm(cl)

# Save Results ####

resPath = paste0(".\\SimResults\\result_", fileName, ".RDS")
saveRDS(rbind(results_2S,results_Iter),resPath)             #STW: Silvennoinen, Terasvirta, Wade

# results <- readRDS(resPath)
# stdReport <- calcStats(results)




