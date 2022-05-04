############################################################
# Cross Validation for lambda tuning parameter
############################################################

cv.genridge = function(x, 
                       y,
                       lambdas=seq(1000, 100000, 5000), 
                       D, 
                       n_folds=10) 
{
  n_lambdas = length(lambdas)
  n = nrow(x[[1]])
  fold_size = n/n_folds
  t = ncol(D)
  p = length(x)
  nrow_test = n/n_folds
  exp_folds = list()
  rand_ind = sample(1:nrow(y), replace=F)
  y = matrix(y[rand_ind], ncol=1)
  
  for (i in 1:p) 
    exp_folds[i] = list(lapply(split(x[[i]][rand_ind,],rep(c(1:n_folds),each=(fold_size))),matrix,fold_size))

  cv_error = rep(0, n_lambdas)
  for (i in 1:n_lambdas) 
  {
    error = rep(0, n_folds)
    counterLow = 1
    counterUp = fold_size
    for (f in 1:n_folds) 
    {
      x_test = lapply(exp_folds, function(x) { return(x[f]) })
      x_train = lapply(exp_folds, function(x) { return(x[-f]) })
      
      x_test = lapply(x_test, function(x) { return(matrix(unlist(x), ncol=t)) })
      x_train = lapply(x_train, function(x) { return(matrix(unlist(x), ncol=t)) })

      y_test = matrix(y[counterLow:counterUp], ncol=1)
      y_train = matrix(y[-c(counterLow:counterUp),], ncol=1)
      
      counterLow = counterUp + 1
      counterUp = fold_size * (f+1)
      

      beta.hat = genridge(x_train, y_train, lambda = lambdas[i], D, int=FALSE)$beta.hat
      
      beta.hat = matrix(beta.hat, ncol=p)
   
      y_pred = rep(0, nrow_test)
      for (j in 1:length(x_train))
      {
        y_pred = y_pred + x_test[[j]] %*% beta.hat[,j]
      }
      
      error[f] = error(y_test, y_pred)
    }
    cv_error[i] = mean(error)
  }
  
  best_lambda = lambdas[which.min(cv_error)]
  
  return(best_lambda)
}









