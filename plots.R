par(mar=c(5, 5, 4, 4) + 0.1)
dev.off()

plot(x=sample_sizes, y=results[[1]][[1]][,2], 
     type="l", ylim=c(0.01, 0.12),
     main="MAE by sample size", xlab="sample size", ylab="MAE", col=2, lty=1)

for (i in 2:ncol(D0.errors))
{
  lines(x=sample_sizes, y=results[[i]][[1]][,2], col=i+1, lty=1)
}

legend("topright", legend = c(paste("D, ord", c(0:5))),
       lwd = c(rep(1, 6)), col = c(c(2:7)), ncol=3, lty=c(rep(1, 6)))

par(new=TRUE)



############# Linear Regression

plotEffectsLR = function(simulation.lm, p, t, ss) 
{
  exp_ind = (p*t+1):((p+1)*t)
  beta_hats = simulation.lm[[2]]
  beta_true = simulation.lm[[3]][[1]][exp_ind, drop=F]
  
  plot(x=1:t, y=rep(0, t), ylim=c(-0.06, 0.06), type="l", lty=2, main="Simulation with OLS | p=5", xlab="Time t", ylab="")
  
  mean_effect = matrix(NA, t, n_sim)
  for(sim in 1:n_sim) 
  {
    b.hat = beta_hats[[sim]][[ss]][exp_ind, drop=F]
    lines(x=1:t, y=b.hat, col=alpha("darkgrey", 0.4))
    mean_effect[, sim] = b.hat
  }
  mean_effect = apply(mean_effect, 1, mean)
  lines(x=1:t, y=mean_effect, col="red", lwd=2)
  lines(x=1:t, y=beta_true, col="black", lwd=2)
  legend("bottom", legend=c("Truth", "Average Estimate", "Estimate"), col=c("black", "red", "darkgrey"), lwd=c(2,2,2),lty = c(1,1,1))
}


# Generalized Ridge Regression

plotEffectsGR = function(results, p, t, D.order, ss) 
{
  exp_ind = (p*t+1):((p+1)*t)
  beta_hats = results[[D.order+1]][[2]]
  beta_true = results[[D.order+1]][[3]][[1]][exp_ind, drop=F]
  plot(x=1:t, y=rep(0, t), type="l", lty=2, ylim=c(-0.1, 0.1), 
       main=paste("Main Effects p = ", p+1," | D order ", D.order), 
       ylab="", xlab="Time t", cex.main=1.6, cex.axis=1.5, cex.lab=1.5)
  
  effects.mat = matrix(NA, t, 50)

  for (sim in 1:length(beta_hats))
  {
    b.hat = beta_hats[[sim]][[ss]][exp_ind, drop=F]
    lines(x=1:t, y=b.hat, col=alpha("darkgrey", 0.4), lty=1)
    effects.mat[, sim] = b.hat
  }
  mean.effect = apply(effects.mat, 1, mean)
  lines(x=1:t, y=mean.effect, lty=1, type="l", col="red", lwd=3)
  lines(x=1:t, y=beta_true, col="black", lwd=3)
  legend("top", legend=c("Truth", "Average Estimate", "Estimate"), col=c("black", "red", "darkgrey"), lwd=c(2,2,2),lty = c(1,1,1), cex=1.2)
}

p = 0
t = 37
D.order = 0
ss = 1

dev.off()
plotEffectsLR(simulation.lm, p, t, ss)

dev.off()
plotEffectsGR(simulation.genridge, p, t, D.order, ss)

par(mfrow=c(1,2))
plotEffectsLR(simulation.ridge, p, t, ss)
plotEffectsLR(simulation.lm, p, t, ss)

dev.off()
par(mfrow=c(2,3))
for (i in 0:5) 
{
  plotEffectsGR(simulation.genridge, p, t, i, ss)
}
plotEffectsLR(simulation.lm, p, t, ss)



dev.off()
graph = getGraph(getD2dSparse(3,15))
plot(graph)



plot(x=1:t, y=simulation.genridge[[1]][[3]][[1]][1:37,], 
     type="l", main="True Effects", xlab="Time t", ylab="", ylim=c(-0.08, 0.03), col="red")
lines(x=1:t, y=simulation.genridge[[1]][[3]][[1]][149:185,], col="blue")
lines(x=1:t, y=rep(0,37), type="l", lty=2, col="darkgrey")
legend("bottomleft", legend=c(paste("Exposure", c(1,5))), col=c("red", "blue"), lty=c(1,1),cex=1.53)
