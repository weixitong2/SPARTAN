#Please install the Microsoft R open for fast calculation, https://mran.microsoft.com/download
#'@import randtoolbox
set.seed(100)
Inv <- function(x){
  l<-length(x)
  (1:l/l-0.5/l)[frank(x)]
}

##################################################################
matpower = function(a,alpha){
  a = round((a + t(a))/2,7); tmp = eigen(a)
  return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))}

##################################################################
hMat = function(data_mat, cov_mat){
  signrt = matpower(cov_mat,-1/2)
  H_temp = Hpi.diag(data_mat%*%signrt)
  H_spin = H_temp%*%matpower(cov_mat,1/2)
  return((H_spin+t(H_spin))/2)
}

##################################################################
# SAVE_dir
SAVE = function(x,y){
  n1=nrow(x); n2=nrow(y); pp=ncol(x)
  data_bind<-rbind(x,y)
  signrt=matpower(cova(data_bind),-1/2)

  cm<-colMeans(data_bind)
  v1<-cova((x- rep(1,n1) %*% t(cm))%*%signrt)
  v2<-cova((y- rep(1,n2) %*% t(cm))%*%signrt)

  savemat=((v1-diag(pp))%*%(v1-diag(pp))+(v2-diag(pp))%*%(v2-diag(pp)))/2
  dir_temp<-as.vector(signrt%*%eigen(savemat)$vectors[,1])
  return(dir_temp/sqrt(crossprod(dir_temp)[1]))
}

#################################################################

itrPDF_EDR <- function(inv_data, itr=50){
  inv_data = as.matrix(inv_data)
  N=nrow(inv_data)
  pp=ncol(inv_data)
  design<-randtoolbox::sobol(N, dim = pp, init = TRUE, scrambling = 1)
  itr_data<-inv_data

  for(i in 1:itr){
    b<- SAVE(itr_data, design)
    ori_proj<-itr_data%*%b
    des_proj<-design%*%b

    mm<-min(des_proj)
    ori_proj_new<-fsort(des_proj-mm)[frank(ori_proj)]+mm
    delta<-ori_proj_new-ori_proj

    itr_data<-itr_data+delta%*%t(b)
  }
  itr_data<-apply(itr_data,2,Inv)
  itr_data
}
#################################################################


