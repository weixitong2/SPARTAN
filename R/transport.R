#'@import transport
emd<-function(inv_data,wasserstein_tyep='2',method='networkflow', fullreturn=FALSE, control = list(), threads=1){
  #method = c("networkflow", "shortsimplex", "revsimplex","primaldual")
  inv_data = as.matrix(inv_data)
  N=nrow(inv_data)
  pp=ncol(inv_data)
  design<-sobol(N, dim = pp, init = TRUE, scrambling = 1)

  mass=matrix(1/N,nrow=1,ncol=N)
  mass_a=as.vector(mass)
  mass_b=as.vector(mass)

  cost_ab=costfunction_r(inv_data,design,wasserstein_tyep)
  ot_plan=transport(mass_a,mass_b,cost_ab,method,fullreturn=FALSE, control = list(), threads=1)
  result=design[ot_plan[ ,2]]

  return (result)
}


costfunction_r<-function(A,B,K='2'){
  n=nrow(A)
  costm=matrix(0,nrow=n,ncol=n)
  if (K !='1'){
    for (i in c(1:n)){
      for (j in c(1:n)){
        c=A[i,]-B[j,]
        costm[i,j]=norm(c,K)^as.integer(K)
      }
    }
  }


  if (K=='1') {
    for (i in c(1:n)){
      for (j in c(1:n)){
        c=A[i,]-B[j,]
        costm[i,j]=sum(abs(c))
      }

    }
  }

  return (costm)
}
