#
# Generate a Simulation Data Matrix:  Noise + tvgarch process
#

library(MTVGARCH)   # Ver. 0.9.5.7

setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

set.seed(1984)
Nr.Series <- 12000

## noiseDist is a named-list describing the error-distribution and parameters - common to both
noiseDist <- list()
noiseDist$name = 'Normal'     
noiseDist$mean = 0  
noiseDist$sd = 1

# Gen Data T2000: ####
Tobs = 2000
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

# Generate series of simulated data, injected with this process:
Ref_Data <- MTVGARCH::generateRefData(nr.series = Nr.Series,nr.obs = Tobs,tvObj = TV,garchObj = GARCH)
saveRDS(Ref_Data,"T2000_Data.RDS")


# Gen Data T5000: ####
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

# Generate series of simulated data, injected with this process:
Ref_Data <- MTVGARCH::generateRefData(nr.series = Nr.Series,nr.obs = Tobs,tvObj = TV,garchObj = GARCH)
saveRDS(Ref_Data,"T5000_Data.RDS")


