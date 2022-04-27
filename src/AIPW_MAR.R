# beta_AIPW-SIM;
# Inverse prob. and estimated by a logistic model, 
# Argu. using SIM under MAR;
# dim(X)=1, dim(U)=3;
# Return:
#   . Beta: a MC*5 matrix storing the estimates of beta
#   . bias: the averaged bias of MC simulations
#   . epsd: empirical SD of MC simulations, based on Beta in (i)
#   . asysd: sqrt of diagonal elements of the averaged Cov matrix over MC simulations
#   . cover: 95% coverages of beta
library(np)              

expit<-function(x){             
  y=exp(x)/(1+exp(x))
  return(y)
}

ip_aug_MAR <- function(xd,inno,alpha,n,seed,trim){
  Beta=matrix(0,MC,5)
  Cov_beta=matrix(0,5,5)
  coverage<-rep(0,5)
  upper_list = vector()
  lower_list = vector()
  for (mc in 1:MC){
    set.seed(seed+mc)
    ww=matrix(0,nrow=n,ncol=5)
    xx=NULL
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
    alpha.fit<-glm(rr~yy+zz1+zz2+zz3,family=binomial())
    alpha_hat<-alpha.fit$coeff
    P_hat<-predict(alpha.fit,type="response")
    P_hat[P_hat <= trim] <- trim
    
    
    S1=0
    S2=0
    Phi=matrix(0,nrow=n,ncol=2)
    theta=c(1,-beta[3],-beta[4],-beta[5])
    Q=cbind(yy,zz1,zz2,zz3)
    sindex=Q%*%theta
    h=1*sd(sindex)*n^(-1/3)
    G1 <- npksum(txdat=sindex,tydat=rr,bws=h)$ksum
    G2 <- npksum(txdat=sindex,tydat=rr,bws=0.4*h)$ksum
    G1[G1==0] <- trim
    G2[G2==0] <- trim
    Phi[,1] <- npksum(txdat=sindex,tydat=rr*xx,bws=h)$ksum/G1
    Phi[,2] <- npksum(txdat=sindex,tydat=rr*(xx^2),bws=0.4*h)$ksum/G2
    for (i in 1:n){
      W=ww[i,]
      y=yy[i]
      R=rr[i]
      Z1=zz1[i]
      Z2=zz2[i]
      Z3=zz3[i]
      phi1=c(1,Phi[i,1],Z1,Z2,Z3)*y
      phi2=c(1,Phi[i,1],Z1,Z2,Z3)%*%t(c(1,Phi[i,1],Z1,Z2,Z3))
      phi2[2,2]=Phi[i,2]
      S1=S1+(R/P_hat[i])*W%*%t(W)+(1-R/P_hat[i])*phi2
      S2=S2+(R/P_hat[i])*W*y+(1-R/P_hat[i])*phi1
    }
    beta_hat=solve(S1)%*%S2
    Beta[mc,]=beta_hat
    
    #print(mc)
  }
  
  bias0=apply(Beta,2,mean)-beta
  epsd0=apply(Beta,2,sd)/sqrt(MC)
  Beta_sort <- apply(Beta,2,sort)
  Beta_sort <- Beta_sort[(0.01*MC+1):(0.99*MC),]
  bias=apply(Beta_sort,2,mean)-beta
  epsd=apply(Beta_sort,2,sd)/sqrt(0.98*MC)
  for (mc in 1:MC){
    beta_hat <- Beta[mc,]
    lower <- beta_hat-1.96*epsd*sqrt(0.98*MC)
    upper <- beta_hat+1.96*epsd*sqrt(0.98*MC)
    lower_list <- append(lower_list,lower)
    upper_list <- append(upper_list,upper)
    coverage <- coverage + (beta>=lower)*(beta<=upper)
  }
  cover=coverage/MC
  return(list(Beta=Beta,bias0=bias0,bias=bias,epsd0=epsd0,epsd=epsd,cover=cover, lower_list=lower_list, upper_list=upper_list))
}
