#'@import T4transport

sample_sinkhorn<-function(inv_data,p=2,lambda=0.02,maxiter=100000000,abstol=1e-40){
  #method = c("networkflow", "shortsimplex", "revsimplex","primaldual")
  inv_data = as.matrix(inv_data)
  N=nrow(inv_data)
  pp=ncol(inv_data)

  design<-sobol(N, dim = pp, init = TRUE, scrambling = 1)


  one=matrix(1,nrow=N,ncol=1)

  skh=sinkhorn(inv_data ,design,p=p ,lambda=lambda,maxiter=maxiter,abstol=abstol )
  ot_plan=skh$plan

  mapping=diag(as.vector(ot_plan%*%one))
  mapping=solve(mapping)%*%ot_plan%*%design
  #mapping=N*ot_plan%*%design


  return (mapping)

}
