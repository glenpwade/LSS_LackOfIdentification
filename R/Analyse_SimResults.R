
library(MTVGARCH)   # Ver. 0.9.9.11

# Set a working directory
setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")
#setwd("C:\\Repos\\LSS_LackOfIdentification")

source("R/calcStats.R")


# Analyse the results:  ----

fileName = "T2000_Low_VarShift"

# Modify the start pars:
TVpars <- c(0.5,1.5,log(10),0.5)
GARCHpars <- c(0.1,0.1,0.8)

# 2-Step
resPath = paste0(".\\SimResults\\result_", fileName, "_2Step.RDS")
results2S <- readRDS(resPath)
# Iterative
resPath = paste0(".\\SimResults\\result_", fileName, "_Iter.RDS")
resultsIter <- readRDS(resPath)

# Clean up any bad data

results <- results2S
results <- resultsIter

# Remove all the failed estimations:
results <- results[results[,11]==0,]
# Remove all the Not-Converged estimations:
results <- results[results[,10]==0,]

# Remove all that hit MaxIter
results <- results[results[,1] < 100,]

# Average Iterations (excl. 2-Step & maxIter):
avgIterations <- mean(results[results[,1] > 2, 1])
medIterations <- median(results[results[,1] > 2, 1])

print("T2000_Low_VarShift:")

options("digits" = 5)

# 2Step:
results2S <- results[results[,1]==1,]
stats <- calcStats(results2S,TVpars,GARCHpars)

#Iterative:
resultsIter <- results[results[,1] > 2,]
stats <- calcStats(resultsIter,TVpars,GARCHpars)

# # LogLik se from actual
avg_Loglik <- mean(results[,9])
avgDeviation_Loglik <- sd(results[,9])



