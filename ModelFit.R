## FPCA function
FPCA = function(sampMat, percent = 0.95){
  # Return
  #       1. eigenvalues
  #       2. eigenfunctions
  #       3. A matrix of projected coefficients, same nrow as sampMat
  ###############################
  
  # TBD
  # list(eigV, eigF, coefMat) # To be returned
}

## load population 
load("data_simu.rda")

##
library(doParallel)

## Simplest Model
Model1 = foreach(B = 1:1000, .combine = 'rbind') %dopar%{
  data_raw = data_simu[sample.int(nrow(data_simu), 200, replace = F), ]  # sample 200 rows from data_simu
  data = cbind(data_raw[,1:5], FPCA(data_raw[,6:ncol(data_raw)])[[3]]) # combine Y,T,W with coefs of X
  
  Ncol = ncol(data)
  glmFit = glm(data[,2] ~ data[,3:Ncol], family = 'binomial')
  p_hat = predict(glmFit, type = "response") # Be careful when some p_hat = 0
  
  tau_HT    = mean(data[,1]*data[,2]/p_hat + data[,1]*(1-data[,2])/(1-p_hat))
  tau_Hajek = sum(data[,1]*data[,2]/p_hat)/sum(data[,2]/p_hat) 
              - sum(data[,1]*(1-data[,2])/(1-p_hat))/sum((1-data[,2])/(1-p_hat))
  c(tau_HT,tau_Hajek)
}

