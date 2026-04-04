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
library(MTVGARCH)   # Ver. 0.9.8.27
library(knitr)
library(foreach)
library(doParallel)

# Set a working directory
#setwd("C:\\Repos\\LSS_LackOfIdentification")

setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

# Set constants
Reps <- 100
Tobs <- 2000

# START SIMS HERE: ####

# Setup the parallel backend ####
#numCores <- parallel::detectCores() - 2
numCores <- 5
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

# Run the Simulation:  ####

# Set the estimation controls to suppress console output
estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startparAdjust=10)

# Get the data:
fileName = "T2000_Med_VarShift"
filePath = paste0(".\\SimSourceData\\",fileName,".RDS")
simData <- readRDS(filePath)

# Create the TV Specification - to be estimated later:
Tobs = NROW(simData)
st = (1:Tobs)/Tobs
shape = tvshape$single
# Create the TV Specification and set starting params to match the loaded Dataset
TVspec <- tv(st,shape)
TVspec$delta0 = 0.5
TVspec$pars["deltaN",1] = 4.0
TVspec$pars["speedN",1] = log(10)
TVspec$pars["locN1",1] = 0.5
TVspec$optimcontrol$parscale <- c(0.5,4.0,2.3,0.5)


GARCHspec <- garch(garchtype$general)
GARCHspec$pars["omega",1] = 0.05           
GARCHspec$pars["alpha",1] = 0.05           
GARCHspec$pars["beta",1]  = 0.90           
GARCHspec$optimcontrol$parscale <- c(0.05,0.05,0.9)
GARCHspec$optimcontrol$ndeps <- c(1e-3,1e-3,1e-3)

#myG <- calculate_g(TVspec)
#plot(myG,type='l')

#myH <- calculate_h(e,GARCHspec)
#plot(myH,type='l')

#processLoglik <- loglik.tvgarch.univar(e,myG,myH)

#rm(mod_2s)

# 2-STEP: ####
cat("\nTWO-STEP:\n")
timestamp()
results = foreach(i=1:Reps, .combine = rbind, .inorder = TRUE, .packages = "MTVGARCH")%dopar%{

#results = foreach(i=91:Reps, .combine = rbind, .inorder = TRUE, .packages = "MTVGARCH", .verbose=TRUE)%do%{

    e = simData[,i]
    
    #cat("\nStarting to do e(n) = ", i, "\n")

    # 1. Set desired Iterations & calc the "true" process LogLik Value
    estCtrl$maxIter <- 1      # Equal to LSS 2-Step
    myG <- calculate_g(TVspec)
    myH <- calculate_h(e,GARCHspec)
    processLoglik <- unname(loglik.tvgarch.univar(e,myG,myH))

    # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification
    mod <- tvgarch(TVspec,GARCHspec)

    # 3. Run the 2-Step estimation
    mod$iterationReltol <- 1e-5    #Hack: This value is used as Threshold for the Iteration Convergence
    mod_2s <- estimateTVGARCH(e,mod,estCtrl)

    # 4. Extract the estimated parameters:
    if(!(is.null(mod_2s))){
        
        tvpars <- mod_2s$Estimated$tv
        garchpars <- mod_2s$Estimated$garch
        # Return
        nr.iterations <- mod_2s@iterations  # Need to avoid the @ in the line below, as the foreach() wrapup (rbind) cannot access the model after the fact
        c(nr.iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,as.numeric(!(mod_2s$Estimated$converged)),(processLoglik - unname(mod_2s$Estimated$value)) )
        
        # Note: Any failed estimations will be identified in the last column, so the unestimated parameters can be excluded
    }else{
        # Just add the starting pars and set failed column = 1, loglik.value difference = Inf
        c(2,mod$tvdelta0,mod$tvpars[1:3],mod$garchpars,1,Inf )
    }
    

}
timestamp()  # ~1 min

# Stop the parallel cluster & remove 'cl'  
# Note: Best to not run this when executing in Parallel Mode / Background Job.  Wait until all tasks/jobs complete, then tidy Up
stopCluster(cl)
rm(cl)

# Save Results ####

resPath = paste0(".\\SimResults\\result_", fileName, "_2S.RDS")
saveRDS(results,resPath)             #STW: Silvennoinen, Terasvirta, Wade

# Identify 2-Step by the IterationCount column1.  Should=1 or 2, depending on 2-Step definition. Or 2+ for Iterative

# results <- readRDS(resPath)
# stdReport <- calcStats(results)
# 



