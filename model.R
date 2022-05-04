# Tasks

# Implement Cross-Validation to select lambda
# Estimate beta hat
# Try different D matrices and estimate beta hat
# Try different data sets and sample sizes
# Visualize MSE by sample size for Linear Regression and Generalized Ridge

# To do 

# change code to make it to estimate beta mse
# remove intercept
# run simulation and figure out which ones work

# restrict to some small sample sizes (over 1000 not needed)
# extend this to more than one exposure

# Plot mean across simulations 
# Calculate MSE of beta hat
# Main interest is beta

# Center outcome and X (subtract mean of column) to make intercept 0

library(mvtnorm)
library(genlasso)

##############################
# Generalized Ridge Regression
##############################

genridge = function(X, Y, lambda, D, int=FALSE) 
{
  Y = Y - mean(Y)

  design = createDesign(X, int)
 
  X = design$design
  D = createD(D, p=design$p)$D
  
  left = solve(crossprod(X, X) + lambda*crossprod(D,D))
  beta.hat = left%*%crossprod(X, Y)
  
  var.beta.hat = sd(Y)^2 * (left %*% crossprod(X, X) %*% left)
  
  return(
    list(
      beta.hat = beta.hat,
      D = D,
      var.beta.hat = var.beta.hat
    )
  )
}

createDesign = function(X, int=FALSE) 
{
  # Dimension
  n = nrow(X[[1]]) # of observations
  t = ncol(X[[1]]) # of time points
  tt = t^2
  p = length(X) # of exposures
  
  # Main
  designMain = matrix(0, nrow=n, ncol=t*p)
  nMainParam = t*p

  # Main
  colCounterLow = 1
  colCounterHigh = t

  for (i in 1:p) 
  {
    designMain[ ,colCounterLow:colCounterHigh] = scale(X[[i]], center = T, scale = F)
    colCounterLow = colCounterHigh + 1
    colCounterHigh = colCounterHigh + t
  }

  if (int) {
    # Interactions
    intList = expand.grid(1:p, 1:p) # grid of exposures
    w = which(intList[,1] < intList[,2])
    intList = intList[w,]
    intComb = choose(p,2) # of exposure interactions
    designInt = array(0, dim=c(n, tt, intComb)) 
    nIntParam = intComb*tt
    nParam = nMainParam + nIntParam
    
    # Interactions
    colCounterLow = t*p + 1
    colCounterHigh = colCounterLow - 1 + tt
    for (i in 1:intComb) 
    {
      j1 = intList[i,1]
      j2 = intList[i,2]
      counter = 1
      for (t1 in 1:t) 
      {
        for (t2 in 1:t) 
        {
          designInt[,counter,i] = designInt[,counter,i] + X[[j1]][,t1]*X[[j2]][,t2]
          designInt[,counter, i] = designInt[,counter, i] - mean(designInt[,counter, i])
          counter = counter + 1
        }
      }
      colCounterLow = colCounterHigh + 1
      colCounterHigh = colCounterHigh + tt
    }
    dim(designInt) = c(n, nIntParam)
    design = cbind(designMain, designInt)
    
    return(
      list(
      design = design,
      designMain = designMain,
      designInt = designInt,
      p = p
    ))
    
  }
  else {
    design = designMain
    return(
      list(
      design = design,
      designMain = designMain,
      p = p
    ))
  }
  
}

createD = function(D, p, int=FALSE) {
  # Dimension
  t = ncol(D) # of time points
  tt = t^2
  
  # D matrix
  nMainParam = t*p
  nD = nrow(D) # of rows in D matrix
  
  
  if(int) {
    # Interactions
    intList = expand.grid(1:p, 1:p) # grid of exposures
    w = which(intList[,1] < intList[,2])
    intList = intList[w,]
    intComb = choose(p,2) # of exposure interactions
    nIntParam = intComb*tt
    nParam = nMainParam + nIntParam
    D2d = getD2d(t, t) # 2d D matrix
    nD2d = nrow(D2d) # of rows in 2d D matrix
    fullD = matrix(0,nrow=(nD*p+nD2d*intComb), ncol=(nParam))
    
    # Interactions
    colCounterLow = t*p + 1
    colCounterHigh = colCounterLow - 1 + tt
    rowCounterLow = nD*p + 1
    rowCounterHigh = rowCounterLow - 1 + nD2d
    for (i in 1:intComb) 
    {
      fullD[rowCounterLow:rowCounterHigh, colCounterLow:colCounterHigh] = D2d
      
      rowCounterLow = rowCounterHigh + 1
      rowCounterHigh = rowCounterHigh + nD2d
      colCounterLow = colCounterHigh + 1
      colCounterHigh = colCounterHigh + tt
    }
    intD = fullD[(nD*p+1):(rowCounterHigh-nD2d), (t*p+1):(colCounterHigh-tt)]
  }
  else {
    nParam = nMainParam
    fullD = matrix(0,nrow=(nD*p), ncol=(nParam))
  }
    
  
  # Main
  colCounterLow = 1
  colCounterHigh = t
  rowCounterLow = colCounterLow
  rowCounterHigh = nD
  for (i in 1:p) 
  {
    fullD[rowCounterLow:rowCounterHigh, colCounterLow:colCounterHigh] = D
    
    colCounterLow = colCounterHigh + 1
    colCounterHigh = colCounterHigh + t
    rowCounterLow = rowCounterHigh + 1
    rowCounterHigh = rowCounterHigh + nD
  }
  
  
  mainD = fullD[1:(nD*p), 1:(t*p)]
  
  
  return(
    list(
      D = fullD,
      mainD = mainD
    )
  )
}



