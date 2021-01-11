#'@import rngWELL
#'@import randtoolbox
#'@import nloptr
#'@import Rfast
#'@import MaxPro
#'@import transport
#'@param data is raw data
#'@param subsize is the size which you want to sample
#'@param method = c('randomprojection', 'PPMM')
#'@param randommethod = c('maxprol','sobol' )
#'@param itr=30, the number is used to solve optimal transport which means the number of iterations
#'@export
#'@examples
#'x = matrix(seq(100),50,2)
#'spartan(x,10)
#'$sample_index
#'[1]  2 13 12 30  9 41 11 46 34 39
#'
#'$sample_point
#'[,1] [,2]
#'[1,]    2   52
#'[2,]   13   63
#'[3,]   12   62
#'[4,]   30   80
#'[5,]    9   59
#'[6,]   41   91
#'[7,]   11   61
#'[8,]   46   96
#'[9,]   34   84
#'[10,]   39   89
#'
spartan <- function(data, subsize, method = 'randomprojection', randommethod = 'maxprol',itr=30){
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

  design<-randtoolbox::sobol(N, dim = pp, init = TRUE, scrambling = 1)
  itr_data<-inv_data

###rp
if(method == 'randomprojection'){
  for(i in 1:itr){

    a<-matrix(rnorm(pp^2),ncol=pp)
    diag(a)<-abs(diag(a))
    b<-qr.Q(qr(a))

    itr_data<-itr_data%*%t(b)


    des<-design%*%t(b)
    for(ii in 1:pp){
      itr_data[,ii]<-des[,ii][order(des[,ii])[rank(itr_data[,ii])]]
    }

    itr_data<-apply(itr_data%*%b,2,Inv)
  }
  new_data<-itr_data
}

  ###PPMM
  if(method == 'PPMM'){
    new_data<-itrPDF_EDR(data,itr)
  }


  #抽样方法
  ###maxpro
  if(randommethod == 'maxprol'){
  design<-MaxPro::MaxProLHD(subsize, pp)$Design}

  ###sobol
  if(randommethod == 'sobol'){
  design<-randtoolbox::sobol(subsize,pp, init = TRUE, scrambling = 1)}

  lhd.ind<-nabor::knn(new_data, design, k=1)$nn.idx[,1]
  new_lhd = unique(lhd.ind)
  size = subsize - length(new_lhd)
  to_sample = setdiff(seq(N), new_lhd)
  if(subsize - length(new_lhd) != 0){
    sample <- sample(to_sample,size,replace=F)
    new_lhd = c(new_lhd,sample)
  }
  lhd.ind <- new_lhd
  #返回值
  return(list(sample_index=lhd.ind, sample_point=data[lhd.ind,]))
}
