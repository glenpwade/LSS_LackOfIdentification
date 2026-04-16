
library(MTVGARCH)   # Ver. 0.9.9.11

# Set a working directory
setwd("C:\\Source\\Repos\\LSS_LackOfIdentification")
#setwd("C:\\Repos\\LSS_LackOfIdentification")

source("R/calcStats.R")


# Analyse the results:  ----

fileName = "T2000_Slow_Speed"
options("digits" = 4)
#print("T2000_Low_VarShift:")

# Modify the start pars:
TVpars <- c(0.5,4.0,1.5,0.5)
GARCHpars <- c(0.05,0.05,0.9)


{# All Matched ----
# 2-Step
resPath = paste0(".\\SimResults\\result_", fileName, "_2Step.RDS")
results2S <- readRDS(resPath)
# Iterative
resPath = paste0(".\\SimResults\\result_", fileName, "_Iter.RDS")
resultsIter <- readRDS(resPath)

# Clean up any bad data

# results2S

# Remove all the failed estimations:
results2S <- results2S[results2S[,11]==0,]
# Remove all the Not-Converged estimations:
results2S <- results2S[results2S[,10]==0,]
# R
# # LogLik se from actual
avg_Loglik <- mean(results2S[,9])
avgDeviation_Loglik <- sd(results2S[,9])


# resultsIter

# Remove all the failed estimations:
resultsIter <- resultsIter[resultsIter[,11]==0,]
# Remove all the Not-Converged estimations:
resultsIter <- resultsIter[resultsIter[,10]==0,]
# Remove all that hit MaxIter
resultsIter <- resultsIter[resultsIter[,1] < 100,]
# Average Iterations (excl. 2-Step & maxIter):
avgIterations <- mean(resultsIter[resultsIter[,1] > 2, 1])
medIterations <- median(resultsIter[resultsIter[,1] > 2, 1])
# # LogLik se from actual
avg_Loglik <- mean(resultsIter[,9])
avgDeviation_Loglik <- sd(resultsIter[,9])

# Now we need to exclude the 2-Step estimations that failed as Iterative
# Filter all 2-Step estimations to match Iterative, by data series (rowname)

# 1. Extract just the numbers 'n' from Iterative rownames
target_names <- rownames(resultsIter)
# 2. Filter TwoStep to only include those rows
valid_names <- target_names[target_names %in% rownames(results2S)]
results2S <- results2S[valid_names, ]
resultsIter <- resultsIter[valid_names, ]

#
stats <- calcStats(results2S,TVpars,GARCHpars)
#
stats <- calcStats(resultsIter,TVpars,GARCHpars)
}


# Consider the series where 2-Step estimations 'succeeded', but Iterative went on to fail:
# Bad 2Steps ----

# 2-Step
resPath = paste0(".\\SimResults\\result_", fileName, "_2Step.RDS")
results2S <- readRDS(resPath)
# Iterative
resPath = paste0(".\\SimResults\\result_", fileName, "_Iter.RDS")
resultsIter <- readRDS(resPath)

resultsIter <- resultsIter[resultsIter[,1] > 100,]
results2S <- results2S[rownames(resultsIter), ]

# Finally:  Remove any errorred series:  NONE

#
stats <- calcStats(results2S,TVpars,GARCHpars)
#
stats <- calcStats(resultsIter,TVpars,GARCHpars)











