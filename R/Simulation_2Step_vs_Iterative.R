if(TRUE)
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
# Install required Packages - once
if(FALSE){
    install.packages("MTVGARCH")  # Install using R-Studio, using the install from file option - MTVGARCH_0.9.5.7.tar.gz provided
    install.packages("knitr")  
}
library(MTVGARCH)   # Ver. 0.9.6.1
library(knitr)
library(foreach)
library(doParallel)

# Set a working directory
setwd("C:\\Repos\\LSS_LackOfIdentification")

# Set constants
Reps <- 3000
Tobs <- 2000
GARCHparScale <- c(0.05,0.05,0.90)
iterRelTol <- 1e-5  #Convergence tolerance for Iterative estimator
estCtrl <- list(calcSE = TRUE, verbose = TRUE) 

# Gen Data ####

# Repeat the block below for each required Dataset
if(FALSE){

    st = (1:Tobs)/Tobs
    shape = tvshape$single
    # Create the tv object with default parameters
    TV <- MTVGARCH::tv(st,shape)
    
    # Now override parameters with desired model specification
    TV$speedopt <- speedopt$eta
    TV$delta0 = 0.5
    TV$pars["deltaN",1] = 1.5
    TV$pars["speedN",1] = log(10)
    TV$pars["locN1",1] = 0.5
    
    # Create the garch object with default parameters
    GARCH = MTVGARCH::garch(garchtype$general)
    # Now override parameters with desired model specification
    GARCH$pars["omega",1] = 0.1
    GARCH$pars["alpha",1] = 0.1
    GARCH$pars["beta",1] = 0.8
    
    ## noiseDist is a named-list describing the error-distribution and parameters
    noiseDist <- list()
    noiseDist$name = 'Normal'     
    noiseDist$mean = 0  
    noiseDist$sd = 1
    set.seed(1984)
    Ref_Data <- MTVGARCH::generateRefData(nr.series = Reps, nr.obs = Tobs, tvObj = TV, garchObj = GARCH)
    
    fileName = ".\\SimSourceData\\T2000_Low_VarShift.RDS"
    saveRDS(Ref_Data,fileName)
    #
    rm(Ref_Data)
}


# Check on 1 Series ####
fileName = ".\\SimSourceData\\T2000_Low_VarShift.RDS"
simData <- readRDS(fileName)

# Check a few plots:
plot(simData[,37],type='l')
plot(simData[,77],type='l')
plot(simData[,49],type='l')
plot(simData[,55],type='l')
plot(simData[,61],type='l')

# Look at series nr. [37,77] hi delta1, [49] lo delta1, [61,55] lo,hi spd

# Quick check that TVGARCH estimation looks ok for 1 series, before starting the full simulation run:
if(FALSE){
    
    e <- simData[,37]

    # Set a seed if you want to control the starting param 'shake':
    set.seed(42)
    
    # Create the TV Specification - to be estimated later:
    st = (1:Tobs)/Tobs
    shape = tvshape$single
    # Create the TV Specification and set starting params to match the loaded Dataset
    TVspec <- MTVGARCH::tv(st,shape)
    TVspec$delta0 = runif(1,0.4,0.6)           # 0.5
    TVspec$pars["deltaN",1] = runif(1,3,5)     # 4.0
    TVspec$pars["speedN",1] = runif(1,2,3)     # 2.3
    TVspec$pars["locN1",1] = runif(1,0.4,0.6)  # 0.5
    
    TVparScale <- c(0.5,4.0,log(10),0.5)
    TVspec$optimcontrol$parscale <- TVparScale
    TVspec$optimcontrol$ndeps <- c(1e-3,1e-5,1e-5,1e-3)
    
    # 1. Do initial Estimation of g(t) assuming h(t) = 1
    TV <- estimateTV(e,TVspec,estCtrl)
    summary(TV)
    plot(TV)

    # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    # 2.1 We need to set the Garch starting pars, before estimating the model:
    mod$garchpars["omega",1] = 0.1  
    mod$garchpars["alpha",1] = 0.1  
    mod$garchpars["beta",1] = 0.8  
    mod$garchOptimcontrol$parscale <- c(0.1,0.1,0.8)
    
    # 3. Since we are only doing one - let's see what's going on & calc parameter se's.
    estCtrl <- list(calcSE = TRUE, verbose = TRUE)  
    
    # 4. Run the 2-Step estimation
    mod_2s <- MTVGARCH::estimateTVGARCH_2Step(e,mod,estCtrl)
    
    # 5. Run the iterative estimation
    # 5.1 But We don't need to calculate statistics for each iteration
    estCtrl <- list(calcSE = FALSE, verbose = TRUE) 
    mod_iter <- MTVGARCH::estimateTVGARCH(e,mod,estCtrl)

    # Check convergence & iterations:
    mod_2s$Estimated$converged
    mod_2s@iterations
    mod_iter$Estimated$converged
    mod_iter@iterations
    
    #View the full estimated object:
    View(mod_2s)
    View(mod_iter)
    
    # View Estimated Params:
    mod <- mod_2s
    cat("2-Step Model Pars, loglik: ", mod$Estimated$value )
    #
    mod <- mod_iter
    cat("Iterative Model Pars, loglik: ", mod$Estimated$value )
    #
    cat("delta0 ", round(mod$Estimated$tv$delta0,3) )
    print(round(mod$Estimated$tv$pars,3) )
    print(round(mod$Estimated$garch$pars,3) )
}

# START SIMS HERE: ####

# Setup the parallel backend ####
#numCores <- parallel::detectCores() - 4
numCores <- 2
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

# Run the Simulation:  ####

# Set the iteration count 
Reps <- 3000    # Simulate 3000 estimations - set lower for debugging any parallel issues
# Set the estimation controls to suppress console output
estCtrl <- list(calcSE = FALSE, verbose = FALSE)

# Get the data:
fileName = "T2000_Low_VarShift"
filePath = paste0(".\\SimSourceData\\",fileName,".RDS")
    
simData <- readRDS(filePath)

# Create the TV Specification - to be estimated later:
Tobs = NROW(simData)
st = (1:Tobs)/Tobs
shape = tvshape$single
# Create the TV Specification and set starting params to match the loaded Dataset
TVspec <- MTVGARCH::tv(st,shape)
TVspec$delta0 = 0.5
TVspec$pars["deltaN",1] = 4
TVspec$pars["speedN",1] = log(10)
TVspec$pars["locN1",1] = 0.5
TVspec$optimcontrol$parscale <- c(0.5,4,log(10),0.5)

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
    
    # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    # 2.1 We need to set the Garch starting pars, before estimating the model:
    mod$garchpars["omega",1] = 0.1  
    mod$garchpars["alpha",1] = 0.1  
    mod$garchpars["beta",1] = 0.8  
    mod$garchOptimcontrol$parscale <- c(0.1,0.1,0.8)
    
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
    
    # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    # 2.1 We need to set the Garch starting pars, before estimating the model:
    mod$garchpars["omega",1] = 0.1  
    mod$garchpars["alpha",1] = 0.1  
    mod$garchpars["beta",1] = 0.8  
    mod$garchOptimcontrol$parscale <- c(0.1,0.1,0.8)
    
    # 3. Run the Iterative estimation
    mod$tvOptimcontrol$reltol <- 1e-5    #Hack: This value is used as Threshold for the Iteration Convergence
    mod_iter <- estimateTVGARCH_Iterate(e,mod,estCtrl)
    
    # 4. Extract the estimated parameters:
    tvpars <- mod_iter$Estimated$tv
    garchpars <- mod_iter$Estimated$garch
    pars[1,] <- c(2,mod_iter@iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,as.numeric(!(mod_iter$Estimated$converged)) )
    
    # Return:
    pars
    
}
timestamp()  # ~20 mins

# Stop the parallel cluster & remove 'cl'
stopCluster(cl)
rm(cl)

resPath = paste0(".\\SimResults\\res_",fileName,".RDS")
saveRDS(rbind(results_2S,results_Iter),resPath)  #STW: Silvennoinen, Terasvirta, Wade
#results_Sim <- rbind(results_2S,results_Iter)


# Analyse the Results:  ####

library(MTVGARCH)   # Ver. 0.9.6.1
library(knitr)

calcStats <- function(resultSet,TVpars,Garchpars) {
    
    # Debug:
    if(FALSE){
        resultSet = results_2S
        resultSet = results_Iter
        TVpars = TVparScale
        Garchpars = GARCHparScale
    }
    #
    
    resultSet <- resultSet[,c(3:9)]  #Extract the parameters
    
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

resfileName = paste0("result_",fileName)
resPath = paste0(".\\SimResults\\", resfileName, ".RDS")
results <- readRDS(resPath)

# Remove any simulation runs that threw an error:
sum(results[,10])

stdReport <- calcStats(results,TVparScale,GARCHparScale)

# Look at series nr. [3,45] lo delta1, [61] hi delta1, lo spd, [87] hi spd




# Any Failures?
NrFails <- NROW(results[results[,10]==0, ])

# Any high iteration/slow converging?
SlowConverges <- results[results[,2] > 10, ]
#
# Conclusion: Slow but Accurate in the end!
