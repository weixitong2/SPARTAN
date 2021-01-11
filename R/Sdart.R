#'@import data.table
#'@import randtoolbox

set.seed(100)
Inv <- function(x){
  l<-length(x)
  (1:l/l-0.5/l)[frank(x)]
}

spin <- function(x){
  eee<-eigen(cov(x))
  vec<-eee$vectors
  val<-eee$values
  t(vec%*%diag(1/sqrt(val))%*%t(vec)%*%t(x))
}

itrPDF <- function(inv_data, itr=50){
  N=nrow(inv_data)
  pp=ncol(inv_data)
  design<-randtoolbox::sobol(N, dim = pp, init = TRUE, scrambling = 1)
  itr_data<-inv_data

  for(i in 1:itr){

    a<-matrix(rnorm(pp^2),ncol=pp)
    diag(a)<-abs(diag(a))
    b<-qr.Q(qr(a))


    itr_data<-itr_data%*%t(b)

    des<-design%*%t(b)
    for(ii in 1:pp){
      mm<-min(des[,ii])
      itr_data[,ii]<-fsort(des[,ii]-mm)[frank(itr_data[,ii])]+mm
    }

    itr_data<-apply(itr_data%*%b,2,Inv)
  }
  itr_data

}
