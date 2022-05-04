source("model.R")
source("data_generator.R")
library(knockoff)
library(doParallel)

powerKnock = FDRandPowerKnock$power
fdrKnock = FDRandPowerKnock$FDR
powerChisq = chisqFDRandPower$power
fdrChisq = chisqFDRandPower$FDR
dev.off()
par(mfrow=c(2,1))
bplot = barplot(powerKnock, main="Power detecting effects via knockoff filter", beside = T, 
        ylim=c(0, 1), names.arg=c("Exposure 1", "Exposure 5"),
        col=c("darkblue","orange", "red"), ylab="Power", legend=c(0.5, 1, 4), 
        args.legend=list(title="Residual variance", cex=1.2), cex.names = 1.5, 
        cex.axis = 1.5, cex.main=1.5, cex.lab=1.5)
text(x=bplot, y = powerKnock, label = powerKnock, pos = 1, cex = 1, col = "white")

bplot2 = barplot(powerChisq, main="Power detecting effects via Chi-squared", beside = T, 
                ylim=c(0, 1), names.arg=c("Exposure 1", "Exposure 5"),
                col=c("darkblue","orange", "red"), ylab="Power", legend=c(0.5, 1, 4), 
                args.legend=list(title="Residual variance", cex=1.2), cex.names = 1.5, 
                cex.axis = 1.5, cex.main=1.5, cex.lab=1.5)
text(x=bplot2, y = powerChisq, label = powerChisq, pos = 1, cex = 1, col = "white")

bplot3 = barplot(fdrKnock, main="Power detecting effects via knockoff filter", beside = T, 
                ylim=c(0, 1),
                col=c("darkblue","orange", "red"), ylab="Power", legend=c(0.5, 1, 4), 
                args.legend=list(title="Residual variance", cex=1.2), cex.names = 1.5, 
                cex.axis = 1.5, cex.main=1.5, cex.lab=1.5)
text(x=bplot3, y = fdrKnock, label = round(fdrKnock,2), pos = 1, cex = 1, col = "white")

bplot4 = barplot(fdrChisq, main="Power detecting effects via knockoff filter", beside = T, 
                 ylim=c(0, 1),
                 col=c("darkblue","orange", "red"), ylab="Power", legend=c(0.5, 1, 4), 
                 args.legend=list(title="Residual variance", cex=1.2), cex.names = 1.5, 
                 cex.axis = 1.5, cex.main=1.5, cex.lab=1.5)
text(x=bplot4, y = fdrChisq, label = fdrChisq, pos = 1, cex = 1, col = "white")

FDRandPowerKnock = simulationKnockOff()

simulationKnockOff = function() {
  
  n_sim=50
  #lambda = 100000
  #d = getDtf(37, 0)
  truth = c(1, 5)
  resid.error = c(0.5, 1, 4)
  meanFDR = matrix(NA, nrow=3, ncol=1)
  meanPower = matrix(NA, nrow=3, ncol=2)
  for(j in 1:length(resid.error)) {
    print(j)
    Power = matrix(NA, nrow=n_sim, ncol=2)
    FDR = matrix(NA, nrow=n_sim, ncol=1)
    for (i in 1:n_sim) {
      data = generateData(n=400, p=5, t=37, resid=resid.error[j])
      x = data$x
      y = data$y
      
      X = createDesign(x, int=F)$design
      #X = x[[2]]
      Y = y
      
      #result = knockoff::create.fixed(X, method="sdp")
      rho = 0.25
      mu = rep(0,ncol(X)); Sigma = toeplitz(rho^(0:(ncol(X)-1)))
      
      X_k = knockoff::create.gaussian(X, mu, Sigma)
      #X_k = knockoff::create.second_order(X, method = "asdp")
      
      #X = result$X
      #X_k = result$Xk
      
      selected = knockoff.procedure(X, X_k, Y, 
                                           statistic = my_knockoff_stat, 
                                           fdr=0.2)
      fdr.i = ifelse(length(selected) == 0, 0, 1 - sum(selected %in% truth) / max(1, length(selected)))
      
      FDR[i,] = fdr.i
      #selected = knockoff_result$selected
      if (truth[1] %in% selected) {
        Power[i, 1] = 1.0
      } else {
        Power[i, 1] = 0
      }
      
      if (truth[2] %in% selected) {
        Power[i, 2] = 1.0
      } else {
        Power[i, 2] = 0
      }
      
    }
    meanPower[j,] = apply(Power, 2, mean)
    meanFDR[j,] = apply(FDR, 2, mean)
  }
  
  
  return(list(
    power = meanPower,
    FDR = meanFDR))
  
}


my_knockoff_stat = function(X, X_k, Y) {
  X = cbind(X, X_k)
  X = scale(X)
  Y = scale(Y)
  #print(dim(X))
  lambda=100000
  D = createD(getD1d(37),p=10, int=F)$mainD
  #my.image.plot(D)
  left = solve(crossprod(X, X) + lambda*crossprod(D,D))
  beta.hat = matrix(left%*%crossprod(X, Y), ncol=2)
  #print(dim(beta.hat))
  beta.hat.p = extract.beta(beta.hat[,1], p=5)
  beta.hat.tilde = extract.beta(beta.hat[,2], p=5)
  #print(beta.hat.p[,,1])
  #print(beta.hat.p[,,5])
  test.stat = c()
  for (i in 1:5) {
    test.stat[i] = sum(beta.hat.p[,,i]^2)^2 - sum(beta.hat.tilde[,,i]^2)^2
  }
  return(test.stat)
}

knockoff.procedure <- function(X, Xk, y, statistic=my_knockoff_stat, fdr=0.2) {
  n = nrow(X); p = ncol(X)
  
  # Compute statistics
  W = statistic(X, Xk, y)
  t = knockoff.threshold(W, fdr=fdr, offset=0)

  selected = sort(which(W >= t))
  
  # Package up the results.
  # structure(list(call = match.call(),
  #                X = X,
  #                Xk = Xk,
  #                y = y,
  #                statistic = W,
  #                threshold = t,
  #                selected = selected),
  #                class = 'knockoff.result')
  return(selected)
}

knockoffThreshold <- function(W, fdr=0.20) {
  
  ts = sort(c(0,abs(W)))
  #print(ts)
  ratio = sapply(ts, function(t) (sum(W <= -t)) / max(1, sum(W >= t)))
  #print(ratio)
  
  ok = which(ratio <= fdr)
  #print(ok)
  
  ifelse(length(ok) > 0, ts[ok[1]], Inf)

}




