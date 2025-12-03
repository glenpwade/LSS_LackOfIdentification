#
# Journal of Time Series Analysis:
# Robust Estimation and Inference for Time-Varying Unconditional Volatility
#
# By, Adam Lee | Rickard Sandberg | Genaro Sucarrat
#
#
#  This paper discredits the MTVGARCH model in a simulation comparison to the authors
#  However the comparison is flawed as the MTVGARCH model has not be correctly implemented.
# 
#  The code below will correctly implement the MTVGARCH model in a similar simulation to provide a fair comparison.


# The proposed simulation is as follows:
#
# 1. Create the TVGARCH Model Spec with the provided parameters
# 2. Generate Simulated data sets injected with this process, 1000 Replications, T = 2000, 5000
# 3. Use MTVGARCH modelling processes to identify & estimate the simulated series
# 4. Report the estimation results

library(MTVGARCH)   # Ver. 0.9.5.7

setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

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


# Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above

# Test 1 Series: ####

if(FALSE){
    
    T2000_Data <- readRDS("T2000_Data.RDS")
    e = T2000_Data[,5]
    plot(e,type='l')
    
    Tobs = 2000
    st = (1:Tobs)/Tobs
    shape = tvshape$single
    # Create the tv object with default parameters, speedopt = eta
    TV <- MTVGARCH::tv(st,shape)
    TV$speedopt <- speedopt$gamma
    TV$optimcontrol$parscale <- c(1,2,20,1)
    
    TV <- estimateTV(e,TV)
    summary(TV)
    plot(TV)
    as.numeric(!(TV$Estimated$error))
    
    mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    
    #estCtrl <- list(calcSE = TRUE, verbose = TRUE)
    mod <- estimateTVGARCH(e,mod,estCtrl,autoConverge = TRUE)
    mod@iterations
    
    mod$Estimated$converged
}



# Run Simulation:  ####

# We want to store the parameters and their standard errors for both TV & GARCH
# We'll just use a matrix with TV$pars, GARCH$pars

T2000_Data <- readRDS("T2000_Data.RDS")

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

# Save the results for reporting:  (Row:1 Two-Step, Row:2 Iterative)
pars <- matrix(NA,1,10)

# 2-STEP:
timestamp()
results_2S = foreach(i=1:1100,.combine = rbind,.inorder = TRUE,.packages = "MTVGARCH")%dopar%{

    e = T2000_Data[,i]

    # Do Two-Step estimation:
    ## Bit of a hack, but close enough:
    # Create the TV Specification and set starting params
    TV <- MTVGARCH::tv(st,shape)
    TV$speedopt <- speedopt$gamma
    TV$delta0 = 0.4
    TV$pars["deltaN",1] = 1.0
    TV$pars["speedN",1] = 5
    TV$pars["locN1",1] = 0.4
    TV$optimcontrol$parscale <- c(1,2,20,1)
    
    # Step1: Estimate TV
    TV <- estimateTV(e,TV,estCtrl)
    if(!TV$Estimated$error){
        # Step2: Estimate Garch on filtered data
        GARCH <- garch(garchtype$general)
        GARCH <- estimateGARCH(e,GARCH,estCtrl,TV)    # Will estimate the h(t) after filtering g(t) from e
        pars[1,] <- c(1,2,TV$Estimated$delta0,TV$Estimated$pars[1:3],GARCH$Estimated$pars,as.numeric(!(TV$Estimated$error)) )
    } else {
        pars[1,] <- c(1,2,TV$Estimated$delta0,TV$Estimated$pars[1:3],GARCH$Estimated$pars,as.numeric(!(TV$Estimated$error)) )
    }
    
    # So far both delta0 and omega have been free params
    
    # # Step3: Estimate TV
    # TV@delta0free <- FALSE
    # TV$optimcontrol$parscale <- c(2,20,1)
    # TV$optimcontrol$ndeps <- c(1e-05,1e-05,1e-05)
    # TV <- estimateTV(e,TV,estCtrl,GARCH)          # Finally estimate TV, with only omega free (wrapup for 2 step)

    # Return:
    pars
    
}
timestamp()



# ITERATIVE:

pars <- matrix(NA,1,10)

timestamp()
results_Iter = foreach(i=1:1100,.combine = rbind,.inorder = TRUE,.packages = "MTVGARCH")%dopar%{
    
    e = T2000_Data[,i]
    
    # Create the TV Specification and set starting params
    TV <- MTVGARCH::tv(st,shape)
    TV$speedopt <- speedopt$gamma
    TV$delta0 = 0.4
    TV$pars["deltaN",1] = 1.0
    TV$pars["speedN",1] = 5
    TV$pars["locN1",1] = 0.4
    TV$optimcontrol$parscale <- c(1,2,20,1)
    
    # Do Iterative estimation
    TV <- estimateTV(e,TV,estCtrl)
    TVG <- MTVGARCH::tvgarch(TV,garchType = garchtype$general) 
    TVG <- estimateTVGARCH(e,TVG,estCtrl,autoConverge = TRUE)
    
    tvpars <- TVG$Estimated$tv
    garchpars <- TVG$Estimated$garch
    pars[1,] <- c(2,TVG@iterations,tvpars$delta0,tvpars$pars[1:3],garchpars$pars,as.numeric(TVG$Estimated$converged) )
    
    # Return:
    pars
    
}
timestamp()


# Stop the parallel cluster
stopCluster(cl)

saveRDS(rbind(results_2S,results_Iter),"T2000_Results.RDS")


# Analyse the Results:  ####

calcStats <- function(resultSet){
    
    resultSet <- resultSet[,c(3:9)]
    
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


results <- readRDS("T2000_Results.RDS")

# Any Failures?
NrFails <- NROW(results[results[,10]==0, ])

# Any high iteration/slow converging?
SlowConverges <- results[results[,2] > 10, ]
#
# Conclusion: Slow but Accurate in the end!


results_2S <- results[results[,1]==1, ]  # Col#1 = 1, 2-Step
results_Iter <- results[results[,1]==2, ]  # Col#1 = 2, Iterative

calcStats(results_2S[1:1000,])

calcStats(results_Iter[1:1000,])

summary(results_Iter[1:1000,(3:9)])  # V1-7: d0,d1,spd,loc, omega,alpha,beta

summary(results_Iter[1:1000,2])  # Col2: Iteration Count
hist(results_Iter[1:1000,2],breaks = 20)  # Col2: Iteration Count



