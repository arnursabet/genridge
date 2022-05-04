error = function(true, predicted) { 
  return(mean(abs(true-predicted))) 
}


my.image.plot = function(m) {
  image.plot(t(apply(m, 2, rev)))
}
