#'@import T4transport
#'#'@import randtoolbox
#'@import nloptr
#'@import Rfast
#'@import MaxPro
#'@import transport
#######################################################################################################
spartan_sinkhorn <- function(data, subsize, wasserstein_tyep='2', method='networkflow', randommethod = 'maxprol',itr=30){
  library(rngWELL)
  library(nloptr)
  library(MaxPro)
  library(randtoolbox)

  pp = ncol(data)
  N = nrow(data)

  #ot map
  Inv=function(x){
    l<-length(x)
    (1:l/l-0.5/l)[rank(x)]
  }

  idt=inv_data=apply(data,2,Inv)

  design<-sobol(N, dim = pp, init = TRUE, scrambling = 1)
  itr_data<-inv_data

  ###sinkhorn
  normal_data<-sample_sinkhorn(data,p=2,lambda=0.02,maxiter=100000000,abstol=1e-40)

  #抽样方法
  ###maxpro
  if(randommethod == 'maxprol'){
    design<-MaxProLHD(subsize, pp)$Design}

  ###sobol
  if(randommethod == 'sobol'){
    design<-sobol(subsize,pp, init = TRUE, scrambling = 1)}
  lhd.ind<-nabor::knn(normal_data, design, k=1)$nn.idx[,1]
  #返回值

  new_lhd = unique(lhd.ind)
  size = subsize - length(new_lhd)
  to_sample = setdiff(seq(N), new_lhd)
  if(subsize - length(new_lhd) != 0){
    sample <- sample(to_sample,size,replace=F)
    new_lhd = c(new_lhd,sample)
  }
  lhd.ind <- new_lhd

  return(list(sample_index=lhd.ind, sample_point=data[lhd.ind,]))
}
