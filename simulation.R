source("model.R")
source("model_selection.R")
source("data_generator.R")

experiment = function(n, t, p, lambdas, D) 
{
  data = generateData(n=n, p=p, t=t,
                      exposures_corr=0.5,
                      time_corr=0.9)
  x = data$x
  y = data$y
  beta.list = data$beta
  beta = matrix(unlist(beta.list), ncol=1)
  
  if (length(lambdas) == 1) 
  {
    best.lambda = lambdas[1]
  } 
  else
  {
    best.lambda = cv.genridge(x, y, lambdas, D, n_folds=5)
  }
  
  beta.hat = genridge(x, y, best.lambda, D, int=F)$beta.hat
  mae = error(beta, beta.hat)
  
  return(
    list(
      "error" = mae,
      "beta.hat" = beta.hat,
      "best.lambda" = best.lambda,
      "beta" = beta
    )
  )
}


simulation1 = function(n_sim, n_vector, t, p, lambdas, D) 
{
  n_points = length(n_vector)
  errors.matrix = matrix(NA, n_points, n_sim)
  betahats = list()
  beta = list()
  for (i in 1:n_sim) 
  {
    if(i%%10 == 0) {print(i)}
    errors.vector = rep(0, n_points)
    lambdas.vector = rep(0, n_points)
    betahats_sim = list()
    for(n in 1:n_points) 
    {
      obj = experiment(n_vector[n], t, p=p, lambdas, D)
      errors.vector[n] = obj$error
      #print(errors.vector)
      lambdas.vector[n] = obj$best.lambda
      betahats_sim[[n]] = obj$beta.hat
      if (i == 1) {beta[[1]] = obj$beta}
    }
    betahats[[i]] = betahats_sim
    errors.matrix[, i] = errors.vector
  }
  results = list()
  results[[1]] = cbind(n_vector, apply(errors.matrix, 1, mean), lambdas.vector)
  results[[2]] = betahats
  results[[3]] = beta
  results[[4]] = errors.matrix
  
  return(results)
}


n_sim = 50
sample_sizes = c(400)
lambdas = c(10, 100, 500, 1000, 2500, 5000, 10000, 100000, 200000)
time.points = 37
exposures = 5
ord.poly = c(0:5)
I = diag(time.points)

# Linear Regression
simulation.lm = simulation1(n_sim = n_sim, 
                            n_vector = sample_sizes,
                            t = time.points, 
                            p = exposures,
                            lambdas = c(0), 
                            D = I) # can be any D matrix since lambda is 0

# Ridge Regression 
simulation.ridge = simulation1(n_sim = n_sim, 
                               n_vector = sample_sizes, 
                               t = time.points, 
                               p = exposures,
                               lambdas = lambdas,
                               D = I)

# Generalized Ridge Regression with D matrix
simulation.genridge = list()

for (i in ord.poly) 
{
  print(i)
  D = getDtf(time.points, ord=i) # Create D matrix
  simulation.obj = simulation1(n_sim = n_sim, 
                               n_vector = sample_sizes, 
                               t = time.points, 
                               p = exposures,
                               lambdas=lambdas, 
                               D = D)
  simulation.genridge[[i+1]] = simulation.obj
}






