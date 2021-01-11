emd<-function(inv_data,wasserstein_tyep='2',method='networkflow'){
  #method = c("networkflow", "shortsimplex", "revsimplex","primaldual")
  inv_data = as.matrix(inv_data)
  N=nrow(inv_data)
  pp=ncol(inv_data)
  design<-sobol(N, dim = pp, init = TRUE, scrambling = 1)

  mass=matrix(1/N,nrow=1,ncol=N)
  mass_a=as.vector(mass)
  mass_b=as.vector(mass)
  Rcpp::sourceCpp('R/costfunction.cpp')
  cost_ab=costfunction(inv_data,design,as.integer(wasserstein_tyep))
  ot_plan=transport(mass_a,mass_b,cost_ab,method)
  result=design[ot_plan[ ,2],]

  return (result)
}
