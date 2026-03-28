
# Install required Packages - once
if(FALSE){
    install.packages("knitr") 
    install.packages("MTVGARCH")  # Install using R-Studio, using the install from file option - MTVGARCH_0.9.5.7.tar.gz provided
    install.packages("GPUmatrix") 
    install.packages("foreach") 
    install.packages("doParallel") 
}

# Unload the tvgarch package to avoid masking issues
try(detach("package:tvgarch", unload = TRUE, character.only = TRUE))
# Load required packages
library(knitr)
library(MTVGARCH)   # Ver. 0.9.5.7
library(GPUmatrix)
library(foreach)
library(doParallel)

# Set a working directory
setwd("C:\\Repos\\LSS_LackOfIdentification")

# 1. Setup Parallel Backend (if using %dopar%)
# Note: Ensure your GPU can handle multiple concurrent contexts/processes
cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

# 2. Define the Objective Function
# This function takes a parameter vector and data, converts data to GPU, 
# and performs the calculation.
gpu_objective <- function(params, data_matrix) {
    # Convert to gpu.matrix if not already (uses 'torch' or 'tensorflow' backend)
    gm <- gpu.matrix(data_matrix, type = "torch") 
    
    
    # Example calculation: (X * params)^2 - a simple matrix operation
    # GPUmatrix supports standard operators like %*%, +, -, etc.
    prediction <- gm %*% params
    loss <- sum((prediction - 1)^2) # Example loss
    
    # Return as a standard R numeric for easy collection
    return(as.numeric(loss))
}

# 3. Create Sample Data
X_data <- matrix(rnorm(100000 * 2000), 100000, 2000)
param_list <- list(runif(2000), runif(2000), runif(2000))

# 4. The foreach Loop
# .packages ensures the worker processes load GPUmatrix
results <- foreach(p = param_list, .packages = "GPUmatrix", .combine = 'c') %dopar% {
    gpu_objective(p, X_data)
}

# Clean up
stopCluster(cl)

print(results)


# 
results <- foreach(i = 1:cores, .packages = c("GPUmatrix", "torch")) %dopar% {
    # This will return TRUE if the worker successfully linked to your GPU
    torch::cuda_is_available()
}
print(results) 
# If this returns FALSE, your workers are falling back to CPU mode.




## torch PKG install

# 1. Remove the current torch package
remove.packages("torch")

# 2. Reinstall it
install.packages("torch")

# 3. Force installation of the GPU (CUDA) binaries
# Your RTX 4060 Ti works best with CUDA 11.8 or 12.1
library(torch)
# If you have CUDA 11.8 installed on your system:
torch::install_torch(type = "11.8") 
# OR let R try to auto-detect:
# torch::install_torch()

options(timeout = 600) # Increase to 10 minutes
torch::install_torch(type = "11.8")


torch::cuda_is_available()


library(torch)
# Replace the paths below with where you saved the downloaded files
torch::install_torch(
    TORCH_URL = "C:/Downloads/libtorch-win-shared-with-deps-2.3.1+cu121.zip",
    LANTERN_URL = "C:/Downloads/Windows-cpu.zip"
)


# Replace these with the actual paths to your downloaded GPU zip files
# Note: Use forward slashes (/) or double backslashes (\\)
Sys.setenv(TORCH_URL = "file:///C:/Downloads/libtorch-win-shared-with-deps-2.3.1+cu121.zip")
Sys.setenv(LANTERN_URL = "file:///C:/Downloads/Windows-cpu.zip")
Sys.unsetenv("TORCH_URL")
Sys.unsetenv("LANTERN_URL")
Sys.setenv(LANTERN_BASE_URL = "https://torch-cdn.mlverse.org/binaries/")

torch::get_install_libs_url()


# Force torch to look for a specific CUDA version
Sys.setenv(CUDA = "12.1") 
torch::install_torch(reinstall = TRUE)


# Note: Use forward slashes (/) or double backslashes (\\)
Sys.setenv(TORCH_URL = "C:/Users/glenp/RLibrary/torch/libs/x64")
Sys.setenv(TORCH_INSTALL = 1)

# Now build the Lantern DLL:

# Point to your local LibTorch build directory
Sys.setenv(TORCH_HOME = "C:/Users/glenp/RLibrary/torch")
Sys.setenv(TORCH_PATH = "C:/Users/glenp/RLibrary/torch")

# Force R to build the Lantern wrapper from source
Sys.setenv(BUILD_LANTERN = 1)
Sys.unsetenv("TORCH_INSTALL")

Sys.getenv()





