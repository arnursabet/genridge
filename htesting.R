
source("data_generator.R")
source("model.R")

chisqFDRandPower = simulation.chisq()



simulation.chisq = function(true_effects, #c(TRUE, FALSE, FALSE, FALSE, TRUE)
                            n=400,
                            p=5,
                            t=37,
                            n_sim=50,
                            k = 0,
                            alpha = 0.05,
                            resid_var = c(0.25, 0.5, 0.75, 1, 4),
                            cv_folds=10,
                            lambdas=NULL,) {
  
  D_k = getDtf(t, k)
  len.resid_var = length(resid_var)
  true_selection = which(isTRUE(true_effects))
  meanFDR = matrix(NA, nrow=len.resid_var, ncol=1)
  meanPower = matrix(NA, nrow=len.resid_var, ncol=sum(truth))
  for(j in 1:len.resid.var) 
  {
    print(j)
    Power = matrix(NA, nrow=n_sim, ncol=2)
    FDR = matrix(NA, nrow=n_sim, ncol=1)
    
    for (i in 1:n_sim) 
    {
      data = generateData(n=n, p=p, t=t, resid=resid_var[j])
      x = data$x
      y = data$y
      
      best.lambda = ifelse(is.null(lambdas), 
                           cv.genridge(x, y, D_k, cv_folds), 
                           cv.genridge(x, y, lambdas, D_k, cv_folds))
      model = genridge(x, y, best.lambda, D_k, int=FALSE)
      
      result = chisqTest(model$var.beta.hat, model$beta.hat, p=p, alpha=alpha)
      reject = as.double(result[1,]) <= alpha
      which.select = which(isTRUE(reject))
      
      fdr.i = ifelse(length(selected) == 0, 
                     0, 
                     1 - sum(which.select %in% true_selection) / max(1, length(which.select)))
      Power[i,] = calculatePower(reject, true_effects)
      FDR[i,] = fdr.i
    }
    meanPower[j,] = apply(Power, 2, mean)
    meanFDR[j,] = apply(FDR, 2, mean)
  }
  
  return(
    list(
      power = meanPower,
      FDR = meanFDR
    ))
}

calculatePower = function(reject, truth) {
  power = c()
  power[1] = reject[1] / truth[1]
  power[2] = reject[5] / truth[5]
  
  return(power)
}


extract.var = function(mat, p=2) {
  t = ncol(mat)/p
  out.array = array(NA, dim=c(t, t, p))
  counterLow = 1
  counterUp = t
  
  for (i in 1:p) {
    out.array[,,i] = mat[counterLow:counterUp, counterLow:counterUp]
    counterLow = counterUp + 1
    counterUp = counterUp + t
  }
  return(out.array)
}

extract.beta = function(vec, p=2) {
  t = length(vec)/p
  out.array = array(NA, dim=c(t, 1, p))
  
  counterLow = 1
  counterUp = t
  
  for(i in 1:p) {
    out.array[,,i] = vec[counterLow:counterUp]
    counterLow = counterUp + 1
    counterUp = counterUp + t
  }
  
  return(out.array)
}


chisqTest = function(var.beta.hat, beta.hat, p=5, alpha=0.05) {
  var.p = extract.var(var.beta.hat, p=p)
  beta.p = extract.beta(beta.hat, p=p)
  
  p.values = numeric(p)
  reject = numeric(p)
  for (j in 1:p) {
    beta.j = matrix(beta.p[,,j], ncol=1)
    var.j = var.p[,,j]
    T.j = t(beta.j) %*% solve(var.j) %*% beta.j
    
    p.val = pchisq(T.j, df=length(beta.j), lower.tail = F)
    
    
    reject[j] = ifelse(p.val <= alpha, "reject", "not reject")
    p.values[j] = signif(p.val, 5)
  }
  names(p.values) = as.character(c(1:p))
  
  return(rbind(p.values, reject))
}





