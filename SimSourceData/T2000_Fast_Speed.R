
# The method is as follows:
#
# 1. Create a TVGARCH Model Spec with the parameters used in the paper. Note: The variance level shift is low (delta0 = 0.5, delta1=1.5) and the Persistence is low.
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
# 1. Generate & save simulated data sets, with 5000 Replications and T = 2000
# 2. Each Simulated data set will be injected with the tv & garch process. 
# 3. Use MTVGARCH modelling processes to identify & estimate the simulated series
# 4. Report the estimation results
# 5. Note: 
#           The starting parameters will be set = process parameters.
#           The default optim() control parameters will be used, with the exception of parscale which is set according to the documented recommendation.
#           We will use 'eta' as the speed option.


# Install required Packages - once
if(FALSE){
    install.packages("MTVGARCH")  # Install using R-Studio, using the install from file option - MTVGARCH_0.9.5.7.tar.gz provided
    install.packages("knitr")  
}
library(MTVGARCH)   # Ver. 0.9.5.7
#library(knitr)

# Set a working directory
setwd("C:\\Repos\\LSS_LackOfIdentification")

# Gen Data ####

# Set constants
Reps <- 5000
Tobs <- 2000

if(TRUE){
    
    st = (1:Tobs)/Tobs
    shape = tvshape$single
    # Create the tv object with default parameters
    TV <- MTVGARCH::tv(st,shape)
    
    # Now override parameters with desired model specification
    TV$speedopt <- speedopt$eta
    TV$delta0 = 0.5
    TV$pars["deltaN",1] = 4.0
    TV$pars["speedN",1] = 5.5
    TV$pars["locN1",1] = 0.5
    
    # Create the garch object with default parameters
    GARCH = MTVGARCH::garch(garchtype$general)
    # Now override parameters with desired model specification
    GARCH$pars["omega",1] = 0.05
    GARCH$pars["alpha",1] = 0.05
    GARCH$pars["beta",1] = 0.90
    
    ## noiseDist is a named-list describing the error-distribution and parameters
    noiseDist <- list()
    noiseDist$name = 'Normal'     
    noiseDist$mean = 0  
    noiseDist$sd = 1
    set.seed(1984)
    Ref_Data <- MTVGARCH::generateRefData(nr.series = Reps, nr.obs = Tobs, tvObj = TV, garchObj = GARCH)
    
    fileName = ".\\SimSourceData\\T2000_Fast_Speed.RDS"
    saveRDS(Ref_Data,fileName)
    
}
