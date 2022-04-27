#  beta_F:Benchmark;
#  Dim(X)=1, dim(U)=3;
#  Return :
#   . Beta: a MC*5 matrix storing the estimates of beta
#   . bias: the averaged bias of MC simulations
#   . epsd: empirical SD of MC simulations, based on Beta in (i)
#   . asysd: sqrt of diagonal elements of the averaged Cov matrix over MC simulations
#   . cover: 95% coverages of beta

expit<-function(x){             
  y=exp(x)/(1+exp(x))
  return(y)
}

full_MAR <- function(xd,inno,alpha,n,seed){
  Beta=matrix(0,MC,5)
  Cov_beta=matrix(0,5,5)
  coverage<-rep(0,5)
  upper_list = vector()
  lower_list = vector()
  for (mc in 1:MC){
    set.seed(seed+mc)
    ww=matrix(0,nrow=n,ncol=5)
    xx=NULL
    wwt=ww
    zz1=NULL
    zz2=NULL
    zz3=NULL
    rr=NULL
    yy=NULL
    P=NULL
    S1=0
    S2=0
    for (i in 1:n){
      X=(xd=="n")*rnorm(1)+(xd=="t")*rt(1,df=5)/sqrt(5/3)+(xd=="g")*(rgamma(1,shape=5,scale=1)-5)/sqrt(5)
      Z1=rnorm(1)
      Z2=rnorm(1)
      Z3=rnorm(1)
      W=c(1,X,Z1,Z2,Z3)
      wwt[i,]=W
      innov=(inno=="n")*rnorm(1)+(inno=="t")*rt(1,df=5)/sqrt(5/3)+(inno=="g")*(rgamma(1,shape=5,scale=1)-5)/sqrt(5)
      y=t(W)%*%beta+innov
      p=expit(t(alpha)%*%c(1,y,Z1,Z2,Z3))
      R=rbinom(1,1,p)
      X=X*R-999*(1-R)
      W=c(1,X,Z1,Z2,Z3)
      ww[i,]=W
      yy=c(yy,y)
      xx=c(xx,X)
      zz1=c(zz1,Z1)
      zz2=c(zz2,Z2)
      zz3=c(zz3,Z3)
      rr=c(rr,R)
    }
    
    for (i in 1:n){
      W=wwt[i,]
      y=yy[i]
      S1=S1+W%*%t(W)
      S2=S2+W*y
    }
    beta_hat <- solve(S1)%*%S2
    Beta[mc,]=beta_hat
    
    Gn=0
    Vn=0
    for (i in 1:n){
      W=wwt[i,]
      y=yy[i]
      U <- W*(as.numeric(y-t(W)%*%beta_hat))
      Gn=Gn+(1/n)*W%*%t(W)
      Vn=Vn+(1/n)*U%*%t(U)
    }
    asycov <- solve(Gn)%*%Vn%*%solve(Gn)
    asyse <- sqrt(diag(asycov)/n)
    lower <- beta_hat-1.96*asyse
    upper <- beta_hat+1.96*asyse
    lower_list <- append(lower_list,lower)
    upper_list <- append(upper_list,upper)
    coverage <- coverage + (beta>=lower)*(beta<=upper)
    Cov_beta <- Cov_beta + asycov
  }
  asysd=sqrt(diag(Cov_beta/MC))
  bias=apply(Beta,2,mean)-beta
  epsd=apply(Beta,2,sd)/sqrt(MC)
  cover=coverage/MC
  return(list(Beta=Beta,bias=bias,epsd=epsd,asysd=asysd,cover=cover, lower_list=lower_list, upper_list=upper_list))
}
