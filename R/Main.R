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
# 2. Generate Simulated data sets injected with this process, 1000 Replications, T = 2000, 4000
# 3. Use MTVGARCH modelling processes to identify & estimate the simulated series
# 4. Report the estimation results

# Model Specification:

library(MTVGARCH)   # Ver. 0.9.5.4

# # TV:
# delta0 <- 10
# delta1 <- 1.5
# speed <- 10        # Note: this is not using our preferred eta-speed scale
# loc1 <- 0.5
# 
# # GARCH:
# omega <- 0.1
# alpha <- 0.1
# beta <- 0.8


# T2000:  ####

Tobs = 2000
st = (1:Tobs)/Tobs
shape = tvshape$single
# Create the tv object with default parameters
TV <- tv(st,shape)


# Now override parameters with desired model specification
TV$speedopt <- speedopt$gamma
TV$delta0 = 10
TV$pars["deltaN",1] = 1.5
TV$pars["speedN",1] = 10
TV$pars["locN1",1] = 0.5

# Create the garch object with default parameters
GARCH = garch(garchtype$general)
# Now override parameters with desired model specification
GARCH$pars["omega",1] = 0.1
GARCH$pars["alpha",1] = 0.1
GARCH$pars["beta",1] = 0.8


# Gen Data ####
# Generate 1000 series of simulated data, injected with this process:
## noiseDist is a named-list describing the error-distribution and parameters
noiseDist <- list()
noiseDist$name = 'Normal'     
noiseDist$mean = 0  
noiseDist$sd = 1
# T2000_Data1 <- generateRefData(nr.series = 1000,nr.obs = 2000,tvObj = TV,garchObj = GARCH)
# saveRDS(T2000_Data,"T2000_Data.RDS")
T2000_Data <- readRDS("T2000_Data.RDS")

# Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above

# Test 1 Series: ####
e = T2000_Data[,1]
TV$optimcontrol$parscale <- c(20,3,20,1)
TV <- estimateTV(e,TV)
summary(TV)
plot(TV)

mtvgarch_mod <- tvgarch(TV,garchType = garchtype$general) 

mtvgarch_mod$tvOptimcontrol$parscale <- c(3,30,1)
mtvgarch_mod$garchOptimcontrol$parscale <- c(1,1,8)

mtvgarch_mod <- estimateTVGARCH(e,mtvgarch_mod)

View(mtvgarch_mod)


# Run Simulation:  ####

# We want to store the parameters and their standard errors for both TV & GARCH
# We'll just use a 1000 x 9 matrix in the same format as the TABLE 1 in LSS paper


# Setup the parallel backend
numCores <- 6
Sys.setenv("MC_CORES" = numCores)
cl <- makeCluster(numCores)
registerDoParallel(cl, cores = numCores)


results = foreach(i=1:1000,.combine = rbind,.inorder = FALSE,.packages = "MTVGARCH")%dopar%{
    
    # Set the estimation controls to suppress console output
    estCtrl <- list(calcSE = TRUE, verbose = FALSE)
    
    e = T2000_Data[,i]
    TV$optimcontrol$parscale <- c(20,3,20,1)
    TV <- estimateTV(e,TV,estCtrl)
    
    mtvgarch_mod <- tvgarch(TV,garchType = garchtype$general) 
    
    mtvgarch_mod$tvOptimcontrol$parscale <- c(3,30,1)
    mtvgarch_mod$garchOptimcontrol$parscale <- c(1,1,8)
    
    # Do first iteration of estimation:
    mtvgarch_mod <- estimateTVGARCH(e,mtvgarch_mod,estCtrl)
    
    # Save these results for reporting:
    # TODO
    
    while(isFALSE(mtvgarch_mod$Estimated$converged)){
        mtvgarch_mod <- estimateTVGARCH(e,mtvgarch_mod,estCtrl)
    }
    # When estimation is complete, write the params and std errors to a 2-row matrix, 1st TV, 2nd GARCH
    
    estTV <- mtvgarch_mod$Estimated$tv
    estGARCH <- mtvgarch_mod$Estimated$garch
    
    pars <- matrix(NA,2,10)
    # TV pars row 1:
    pars[1,] <- c(1,estTV$delta0,estTV$delta0_se,estTV$pars[1,1],estTV$se[1,1],estTV$pars[2,1],estTV$se[2,1],estTV$pars[3,1],estTV$se[3,1],mtvgarch_mod$Estimated$iteration)
    
    
    
}

# Stop the parallel cluster
stopCluster(cl)

