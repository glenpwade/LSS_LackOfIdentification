{#Introduction  ####
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
library(MTVGARCH)   # Ver. 0.9.8.53
library(knitr)
library(foreach)
library(doParallel)

# Set a working directory
setwd("C:\\Repos\\LSS_LackOfIdentification")

#setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

# Set constants
Reps <- 3000
Tobs <- 2000

# START SIMS HERE: ####

# Setup the parallel backend ####
#numCores <- parallel::detectCores() - 2
numCores <- 10
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)

# Run the Simulation:  ####

# Get the data:
fileName = "T2000_Fast_Speed"
filePath = paste0(".\\SimSourceData\\",fileName,".RDS")
simData <- readRDS(filePath)

# Create the TV Specification - to be estimated later:
Tobs = NROW(simData)
st = (1:Tobs)/Tobs
shape = tvshape$single
# Create the TV Specification and set starting params to be near (but not too close) the process values
TVspec <- tv(st,shape)
TVspec$delta0 = 1.0
TVspec$pars["deltaN",1] = 2.0
TVspec$pars["speedN",1] = 3.0
TVspec$pars["locN1",1] = 0.66
TVspec$optimcontrol$parscale <- c(0.5,4.0,5.5,0.5)
TVspec$optimcontrol$ndeps <- c(1e-5,1e-5,1e-5,1e-5)

# Create the GARCH Specification and set starting params to be near (but not too close) the process values
GARCHspec <- garch(garchtype$general)
GARCHspec$pars["omega",1] = 0.15           
GARCHspec$pars["alpha",1] = 0.15           
GARCHspec$pars["beta",1]  = 0.70           
GARCHspec$optimcontrol$parscale <- c(0.05,0.05,0.90)    # Keep parScale = process params.
GARCHspec$optimcontrol$ndeps <- c(1e-5,1e-5,1e-5)

# ITERATIVE: ####
cat("\nITERATIVE:\n")
timestamp()
results_Iter = foreach(i=1:Reps, .combine = rbind, .inorder = TRUE, .packages = "MTVGARCH")%dopar%{    
    
    # Set the estimation controls to suppress console output
    estCtrl <- list(calcSE=FALSE, verbose=FALSE, maxIter=100, fixStartPars=FALSE, startparAdjust=10)
    
    # Attempt the estimation
    mod <- tryCatch({
        # 1. Set desired Iterations & calc the "true" process LogLik Value
        estCtrl$maxIter <- 100      
        myG <- calculate_g(TVspec)
        myH <- calculate_h(simData[,i],GARCHspec)
        processLoglik <- unname(loglik.tvgarch.univar(simData[,i],myG,myH))
        
        # 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification
        modTVG <- tvgarch(TVspec,GARCHspec)
        
        # 3. Run the 2-Step estimation
        modTVG$iterationReltol <- 1e-5    #Hack: This value is used as Threshold for the Iteration Convergence
        mod <- estimateTVGARCH(simData[,i],modTVG,estCtrl)
        
    }, error = function(e) {
        # If a hard error occurs, return a placeholder with the flag set
        message(paste("Error in series", i, ":", e$message))
        return(NULL) 
    })
    
    # Check if estimation succeeded and converged
    if (is.null(mod)) {
        # FAILED
        # Return NA's to keep rbind happy & set ConvergeError (col 10) =1, Error (col 11) =1
        return(c(0,rep(NA,8),1,1 ))
        
    } else {
        # Estimation succeeded and converged
        tvpars <- mod$Estimated$tv
        garchpars <- mod$Estimated$garch
        nr.iterations <- mod@iterations  # May need to avoid the @ in the line below, as the foreach() wrapup (rbind) cannot access the model after the fact
        # Return
        return(c(nr.iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,(processLoglik - unname(mod$Estimated$value)),as.numeric(!(mod$Estimated$converged)),as.numeric(mod$Estimated$error) ))
    }
    
}
timestamp()

# Stop the parallel cluster & remove 'cl'  
# Note: Best to not run this when executing in Parallel Mode / Background Job.  Wait until all tasks/jobs complete, then tidy Up
stopCluster(cl)
rm(cl)


# Save the results for reporting:  (Col:1 IterationCount), (Col10: 'EstimationError': 0(FALSE) / 1(TRUE))
results <- results_Iter
colnames(results) <- c("NrIterations","d0","d1","spd","loc","omega","alpha","beta","loglikDeviation","NotConverged","EstError")

resPath = paste0(".\\SimResults\\result_", fileName, "_Iter.RDS")
saveRDS(results,resPath)             #STW: Silvennoinen, Terasvirta, Wade


# Analyse the results:  ----

#  Took 5:40 hr:min - 6 cores

# source("R\\calcStats.R")
# fileName = "T2000_Low_VarShift"
# 
# resPath = paste0(".\\SimResults\\result_", fileName, ".RDS")
# results <- readRDS(resPath)
# # Index the rownames function, not the results matrix
# rownames(results)[1:3000] <- paste0("two_Steps.", 1:3000)
# rownames(results)[3001:6000] <- paste0("Iterative.", 1:3000)
# ## Note: the names above must be the same length!!  Will rely on this later.
# 
# # Modify the start pars:
# TVpars <- c(0.5,1.5,log(10),0.5)
# GARCHpars <- c(0.10,0.10,0.80)
# 
# ## -- Only analyse the data that successfully estimated both 2-Step & Iterative:  ----
# 
# # Grab every row that starts with "Iterative"
# resultsIter <- results[grep("^Iterative\\.", rownames(results)), ]
# # Remove all the failed estimations:
# resultsIter <- resultsIter[resultsIter[,11]==0,]  # Keep rows where Error = 0
# # Remove all that failed to converge
# resultsIter <- resultsIter[resultsIter[,10]==0,]  # Keep rows where NotConverged = 0
# # Remove all that hit MaxIter
# resultsIter <- resultsIter[resultsIter[,1] < 100,]
# 
# 
# # Grab every row that starts with "two_Steps"
# results2S <- results[grep("^two_Steps\\.", rownames(results)), ]
# # Remove all the failed estimations:
# results2S <- results2S[results2S[,11]==0,]  # Keep rows where Error = 0
# # Remove all that failed to converge
# results2S <- results2S[results2S[,10]==0,]  # Keep rows where NotConverged = 0
# 
# # Now we need to exclude the 2-Step estimations that failed as Iterative
# # Filter all 2-Step estimations to match Iterative, by data series (rowname)
# 
# # 1. Extract just the numbers 'n' from Iterative rownames
# # This removes "Iterative." and leaves "1", "2", etc.
# nums <- sub("Iterative\\.", "", rownames(resultsIter))
# # 2. Reconstruct the corresponding twoStep names
# target_names <- paste0("two_Steps.", nums)
# # 3. Filter TwoStep to only include those rows
# valid_names <- target_names[target_names %in% rownames(results2S)]
# results2S <- results2S[valid_names, ]
# # Conform teh Iterative:
# nums <- sub("two_Steps\\.", "", rownames(results2S))
# target_names <- paste0("Iterative.", nums)
# valid_names <- target_names[target_names %in% rownames(resultsIter)]
# resultsIter <- resultsIter[valid_names, ]
# 
# 
# # Generate simulation Summary Stats (2Step & Iter)
# EstErrorCount <- sum(results2S[,11], na.rm = TRUE)
# FailedConvergeCount <- sum(results2S[,10], na.rm = TRUE)
# HitMaxIterCount <- NROW( results2S[results2S[,1] >= 100,] )
# simSummary2S <- c(NROW(results2S),EstErrorCount,FailedConvergeCount,HitMaxIterCount)
# #
# EstErrorCount <- sum(resultsIter[,11], na.rm = TRUE)
# FailedConvergeCount <- sum(resultsIter[,10], na.rm = TRUE)
# HitMaxIterCount <- NROW( resultsIter[resultsIter[,1] >= 100,] )
# simSummaryIter <- c(NROW(resultsIter),EstErrorCount,FailedConvergeCount,HitMaxIterCount)
# 
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# simSummary2S <- c(simSummary2S, NROW(results2S))
# simSummaryRes(simSummary2S,type="Two-Step ")
# #
# simSummaryIter <- c(simSummaryIter, NROW(resultsIter))
# simSummaryRes(simSummaryIter,type="Iterative")
# 
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# 
# 
# print(fileName)
# 
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# 
# 
# # 2Step:
# stats <- calcStats(results2S,TVpars,GARCHpars)
# 
# #Iterative:
# stats <- calcStats(resultsIter,TVpars,GARCHpars)
# 
# # LogLik se from actual
# avg_Loglik <- mean(results[,9],na.rm = TRUE)
# avgDeviation_Loglik <- sd(results[,9],na.rm = TRUE)
# #
# 
# 
# 
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# 
# # Remove all the failed estimations:
# resultsIter <- resultsIter[resultsIter[,11]==0,]  # Keep rows where Error = 0
# 
# # Remove all that failed to converge
# resultsIter <- resultsIter[resultsIter[,10]==0,]  # Keep rows where NotConverged = 0
# 
# # Remove all that hit MaxIter
# resultsIter <- resultsIter[resultsIter[,1] < 100,]
# 
# # Average Iterations (excl. 2-Step & maxIter):
# avgIterations <- mean(resultsIter[resultsIter[,1] > 2, 1], na.rm = TRUE)
# 
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
# 
# # Filter out all 2-Step estimations that didn't succeed as an Iterative
# results2S <- results2S[rownames(resultsIter), ]
# 
# # Remove all the failed estimations:
# results2S <- results2S[results2S[,11]==1,]  # Keep rows where Error = 0
# 
# # Remove all that failed to converge
# results2S <- results2S[results2S[,10]==0,]  # Keep rows where NotConverged = 0
# 
# 
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #