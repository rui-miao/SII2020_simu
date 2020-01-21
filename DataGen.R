## Simulation Data Generation

propScoreM1 = function(A, W){ # Or a function of Z
  # propensity score model to generate T
  ###############################
  # A --- a vec of length 6
  # W --- a vec of length 3
  ###############################
  
  # Model1 only linear part
  -W[1] + 0.5*W[2] - 0.1*W[3] + 2*A[1] + 0.5*A[2] + 0.5*A[3] + A[4]
}

propScoreM2 = function(A, W){ # Or a function of Z
  # propensity score model to generate T
  ###############################
  # A --- a vec of length 6
  # W --- a vec of length 3
  ###############################
  
  # Model2 only linear part
  # -W[1] + 0.5*W[2] - 0.1*W[3] + ......
}

propScoreM3 = function(Z){ # Or a function of Z
  # propensity score model to generate T
  ###############################
  # Z --- a vec of length 6
  ###############################
  
  # Model3 only linear part
  -Z[1] + 0.5*Z[2] - 0.25*Z[3] - 0.1*Z[4]
}

treatEff = function(T, Z){ # Or a function of Z
  # treatment effect model to generate Y
  ###############################
  # T --- binary variable 0 or 1
  # A --- a vec of length 6
  # W --- a vec of length 3
  ###############################
  
  200 + 10*T + (1.5*T - 0.5) * (27.4*Z[1] + 13.7*Z[2] + 13.7*Z[3] + 13.7*Z[4]) + rnorm(1)
  
}

funcCov = function(A,delta,h){
  # generate a smooth curve on [0,1] by coefficients A and a basis, add gaussian noise with SNR=1/delta
  ###############################
  # A --- a vec of length 6
  # delta --- a small positive real number
  # h --- bandwidth
  ###############################
  t = seq(from = 0, to = 1, by = h)
  X = sqrt(2) * (A[1] * sin(2*pi*t) + A[2] * cos(2*pi*t) 
                 + A[3] * sin(4*pi*t) + A[4] * cos(4*pi*t) 
                 + A[5] * sin(8*pi*t) + A[6] * cos(8*pi*t))
  X = X + rnorm(length(t), mean = 0, sd = delta * sd(X))
}

## Generate a large population, we will sample 1000 times of size n=200 or n=500
data_simu_M1 = t(sapply(1:100000, function(k){
  Z  = rnorm(6)
  A  = 2*Z/(1:6)
  W1 = Z[1] + 2*Z[2]         # Or other function of Z s.t. E W1 = 0 
  W2 = Z[2]^2 -Z[3]^2        # Or other function of Z s.t. E W2 = 0 
  W3 = exp(Z[3]) - exp(0.5)  # Or other function of Z s.t. E W3 = 0 
  W  = c(W1,W2,W3)
  T  = rbinom(n = 1, size = 1, prob = boot::inv.logit(propScoreM1(A,W)))
  Y  = treatEff(T, Z)
  X  = funcCov(A, delta = 0.1, h = 0.02)
  c(Y,T,W,X)
}))

data_simu_M2 = t(sapply(1:100000, function(k){
  Z  = rnorm(6)
  A  = 2*Z/(1:6)
  W1 = Z[1] + 2*Z[2]         # Or other function of Z s.t. E W1 = 0 
  W2 = Z[2]^2 -Z[3]^2        # Or other function of Z s.t. E W2 = 0 
  W3 = exp(Z[3]) - exp(0.5)  # Or other function of Z s.t. E W3 = 0 
  W  = c(W1,W2,W3)
  T  = rbinom(n = 1, size = 1, prob = boot::inv.logit(propScoreM2(A,W)))
  Y  = treatEff(T, Z)
  X  = funcCov(A, delta = 0.1, h = 0.02)
  c(Y,T,W,X)
}))

data_simu_M3 = t(sapply(1:100000, function(k){
  Z  = rnorm(6)
  A  = 2*Z/(1:6)
  W1 = Z[1] + 2*Z[2]         # Or other function of Z s.t. E W1 = 0 
  W2 = Z[2]^2 -Z[3]^2        # Or other function of Z s.t. E W2 = 0 
  W3 = exp(Z[3]) - exp(0.5)  # Or other function of Z s.t. E W3 = 0 
  W  = c(W1,W2,W3)
  T  = rbinom(n = 1, size = 1, prob = boot::inv.logit(propScoreM3(Z)))
  Y  = treatEff(T, Z)
  X  = funcCov(A, delta = 0.1, h = 0.02)
  c(Y,T,W,X)
}))

save(data_simu_M1,data_simu_M2,data_simu_M3, file = "data_simu.rda")
