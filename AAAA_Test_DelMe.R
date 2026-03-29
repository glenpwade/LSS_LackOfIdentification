# START ####
library(MTVGARCH)
library(knitr)
library(foreach)
library(doParallel)

# Set a working directory
#setwd("C:\\Repos\\LSS_LackOfIdentification")
setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")

# Set constants
Reps <- 20
Tobs <- 2000
GARCHparScale <- c(0.05,0.05,0.90)
iterRelTol <- 1e-5  #Convergence tolerance for Iterative estimator


# Check on 1 Series ####
fileName = ".\\SimSourceData\\T2000_Med_VarShift.RDS"
simData <- readRDS(fileName)

# Check a few plots:
plot(simData[,37],type='l')
plot(simData[,77],type='l')
plot(simData[,49],type='l')
plot(simData[,55],type='l')
plot(simData[,61],type='l')

e <- simData[,37]

# Set a seed if you want to control the starting param 'shake':
set.seed(42)

# Create the TV Specification - to be estimated later:

# ---  TV  --- #

st = (1:Tobs)/Tobs
shape = tvshape$single
# Create the TV Specification and set starting params to match the loaded Dataset
TVspec <- tv(st,shape)
TVspec$delta0 = runif(1,0.4,0.6)           # 0.5
TVspec$pars["deltaN",1] = runif(1,3,5)     # 4.0
TVspec$pars["speedN",1] = runif(1,2,3)     # 2.3
TVspec$pars["locN1",1] = runif(1,0.4,0.6)  # 0.5

TVparScale <- c(0.5,4.0,log(10),0.5)
TVspec$optimcontrol$parscale <- TVparScale
TVspec$optimcontrol$ndeps <- c(1e-3,1e-5,1e-5,1e-3)

# 1. Do initial Estimation of g(t) assuming h(t) = 1 AND delta0free = ON (default)
TVspec@delta0free <- TRUE
estCtrl <- list(calcSE = TRUE, verbose = TRUE)
TV <- estimateTV(e,TVspec,estCtrl)
summary(TV)
plot(TV)

# 2. Do initial Estimation of g(t) assuming h(t) = 1 AND delta0free = OFF
TVspec@delta0free <- FALSE
#estCtrl <- list(calcSE = FALSE, verbose = TRUE)
TV <- estimateTV(e,TVspec,estCtrl)
summary(TV)

# ---  GARCH  --- #

GARCHspec <- garch(garchtype$general)
GARCHspec$pars["omega",1] = runif(1,0.04,0.06)           # 0.50
GARCHspec$pars["alpha",1] = runif(1,0.04,0.06)           # 0.50
GARCHspec$pars["beta",1]  = runif(1,0.80,0.99)           # 0.90
GARCHspec$optimcontrol$parscale <- c(0.05,0.05,0.9)
GARCHspec$optimcontrol$ndeps <- c(1e-5,1e-5,1e-5)

# 1. Do initial Estimation of g(t) assuming h(t) = 1 AND omegafree = ON (default)
GARCHspec@omegafree <- TRUE
estCtrl <- list(calcSE = TRUE, verbose = TRUE)
GARCH <- estimateGARCH(e,GARCHspec,estCtrl)
summary(GARCH)

# 1. Do initial Estimation of g(t) assuming h(t) = 1 AND omegafree = OFF
GARCHspec@omegafree <- FALSE
#estCtrl <- list(calcSE = FALSE, verbose = TRUE)
GARCH <- estimateGARCH(e,GARCHspec,estCtrl)
summary(GARCH)

# TVGARCH Estimation:  ####

# 2. Specify a multiplicitive TV GARCH model specification using the TV & GARCH specification above
mod <- MTVGARCH::tvgarch(TV,garchType = garchtype$general)
# 2.1 We need to set the Garch starting pars, before estimating the model:
mod$garchpars["omega",1] = 0.05
mod$garchpars["alpha",1] = 0.05
mod$garchpars["beta",1] = 0.90
mod$garchOptimcontrol$parscale <- c(0.05,0.05,0.90)

# 3. Since we are only doing one - let's see what's going on & calc parameter se's.
estCtrl <- list(calcSE = TRUE, verbose = TRUE)

# 4. Run the 2-Step estimation
mod_2s <- MTVGARCH::estimateTVGARCH_2Step(e,mod,estCtrl)

# 5. Run the iterative estimation
# 5.1 But We don't need to calculate statistics for each iteration
#estCtrl <- list(calcSE = FALSE, verbose = TRUE)
mod_iter <- MTVGARCH::estimateTVGARCH_Iterate(e,mod,estCtrl)



