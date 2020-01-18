## Simulation Data Generation

propScore = function(A, W){ # Or a function of Z
  # propensity score model to generate T
  ###############################
  # A --- a vec of length 6
  # W --- a vec of length 3
  ###############################
  
  # TBD
}

treatEff = function(T, A, W){ # Or a function of Z
  # treatment effect model to generate Y
  ###############################
  # T --- binary variable 0 or 1
  # A --- a vec of length 6
  # W --- a vec of length 3
  ###############################
  
  #TBD
}

funcCov = function(A,delta,h){
  # generate a smooth curve on [0,1] by coefficients A and a basis, add gaussian noise with SNR=1/delta
  ###############################
  # A --- a vec of length 6
  # delta --- a small positive real number
  # h --- bandwidth
  ###############################
  t = seq(from = 0, to = 1, by = h)
  X = A[1] * sin(t) + A[2] * cos(t) + A[3] * sin(2*t) + A[4] * cos(2*t) + A[5] * sin(4*t) + A[6] * cos(4*t)
  X = X + rnorm(length(t), mean = 0, sd = delta * sd(X))
}

## Generate a large population, we will sample 1000 times of size n=200
data_simu = t(sapply(1:100000, function(k){
  Z  = rnorm(6)
  A  = 4*Z/(1:6)
  W1 = Z[1] + 2*Z[2]         # Or other function of Z s.t. E W1 = 0   #TBD
  W2 = Z[3]^2 -Z[4]^2        # Or other function of Z s.t. E W2 = 0   #TBD
  W3 = Z[2] + Z[5]^3         # Or other function of Z s.t. E W3 = 0   #TBD
  W  = c(W1,W2,W3)
  T  = rbinom(n = 1, size = 1, prob = boot::inv.logit(propScore(A,W)))
  Y  = treatEff(T, A, W)
  X  = funcCov(A, delta = 0.1, h = 0.01)
  c(Y,T,W,X)
}))


