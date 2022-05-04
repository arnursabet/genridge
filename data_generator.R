library(mvtnorm)

generateData = function(n=200, p=5, t=15, 
                        exposures_corr=0.5, time_corr=0.95, 
                        beta=NULL, betaInt=NULL, resid=NULL) {
  
  # intList = expand.grid(1:p, 1:p)
  # w = which(intList[,1] < intList[,2])
  # intList = intList[w,]
  
  if (is.null(beta)) 
  {
    beta = list()
    
    for (j in 1:p) 
    {
      beta[[j]] = rep(0, t)
    }
    beta[[1]] = (((1:t)/t)*0.25 - 0.4*((1:t)/t)^2) / 2
    #beta[[1]] = c(.0075, 0.01, 0.015, 0.02, 0.03, 
     #             0.04, 0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0610, 0.0615, 0.0620)
    beta[[5]] = c(rep(0,19), 0.00125, 0.0025, 0.005, .0075, 0.01, 0.015, 0.02, 0.03, 
                  0.04, 0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0610, 0.0615, 0.0620) /2
  }
  
  # if (is.null(betaInt)) {
  #   betaInt = list()
  #   for (j in 1 : choose(p, 2)) {
  #     betaInt[[j]] = matrix(0, t,t)
  #   }
  #   betaInt[[1]] = 0.7*outer(((1:t)/t)*0.6 - 0.8*((1:t)/t)^2,
  #                            ((1:t)/t)*0.7 - 0.7*((1:t)/t)^2) / 7
  # } 
  
  x = list()
  
  for (j in 1:p) 
  {
    x[[j]] = matrix(NA, n, t)
  }
  
  Amat = diag(time_corr, p)
  
  sigma = matrix(NA, p, p)
  
  for (row in 1:p) 
  {
    for (col in 1:p) 
    {
      sigma[row, col] = exposures_corr^(abs(row - col))
    }
  }
  
  exposuresData = array(NA, dim=c(n, p, t))
  exposuresData[,,1] = rmvnorm(n, Sigma=sigma)
  
  for (i in 1:n) 
  {
    for (timePoint in 2:t) 
    {
      exposuresData[i,,timePoint] = Amat %*% exposuresData[i,,timePoint-1] + as.numeric(rmvnorm(1, Sigma=sigma))
    }
  }
  
  for (j in 1:p) 
  {
    for (tt in 1:t) 
    {
      x[[j]][,tt] = (exposuresData[, j, tt] - mean(exposuresData[,j,tt])) / sd(exposuresData[,j,tt])
    }
  }
  
  LinPred = rep(0, n)
  for (j in 1:p)
  {
    LinPred = LinPred + x[[j]] %*% beta[[j]]
  }
  
  # LinPredInt = rep(0,n)
  # for (j in 1 : choose(p,2)) {
  #   j1 = intList[j,1]
  #   j2 = intList[j,2]
  #   for (t1 in 1 : t) {
  #     for (t2 in 1 : t) {
  #       LinPredInt = LinPredInt + x[[j1]][,t1]*x[[j2]][,t2] * betaInt[[j]][t1,t2]
  #     }
  #   }
  # }
  
  if(is.null(resid)) {
    resid = 1
  } 
  
  #y = LinPred + LinPredInt + rnorm(n, sd=sqrt(residError))
  y = LinPred + rnorm(n, sd=sqrt(resid))
  
  
  return(
    list(
      y = y, 
      x=x, 
      beta=beta, 
      betaInt=betaInt
  ))
}
