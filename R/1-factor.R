
# TP: a matrix with columns the TPs for various tests
# TN: a matrix with columns the TNs for various tests
# FN: a matrix with columns the FNs for various tests
# FP: a matrix with columns the FPs for various tests
wbinomprod=function(x,gl,TP,FN,FP,TN)
{ n1=TP+FN
  n2=TN+FP
  K=length(TP)
  nq=nrow(x)
  f1=f2=matrix(NA,nq,K)
  for(j in 1:K)
  { f1[,j]=gl$w %*% dbinom(TP[j],size=n1[j],prob=x[,,j]) 
    f2[,j]=gl$w %*% dbinom(TN[j],size=n2[j],prob=x[,,(K+j)])
  }
  apply(cbind(f1,f2),1,prod)
}



###################################################

factloglik.norm<-function(param,TP,FN,FP,TN,
                              gl,mgrid,qcond1,tau2par1,
                          qcond2,tau2par2)
{ K=ncol(TP)
  K2=2*K
  p=param[1:K2]
  si=param[(K2+1):(2*K2)]
  tau=param[(2*K2+1):(3*K2)]
  if(any(p<=0) | any(p>=1)) return(1.e10)
  if(any(si<=0)) return(1.e10)
  if(any(tau< -0.95) | any(tau>0.95)) return(1.e10)
  
  mu=log(p/(1-p))
  th=c(sapply(tau[1:K],tau2par1),sapply(tau[(K+1):K2],tau2par2))
  
  nq=length(gl$nodes)
  uu=array(NA,c(nq,nq,K2))
  for(j in 1:K)
  { u=qcond1(mgrid$y,mgrid$x,th[j]) 
    u[u>1]=0.9999
    x=qnorm(u,mu[j],si[j])
    t=exp(x) 
    uu[,,j]=t/(1+t) 
  }
  for(j in (K+1):K2)
  { u=qcond2(mgrid$y,mgrid$x,th[j]) 
    u[u>1]=0.9999
    x=qnorm(u,mu[j],si[j])
    t=exp(x) 
    uu[,,j]=t/(1+t) 
  }
  
  N=nrow(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=wbinomprod(uu,gl,TP[i,],FN[i,],FP[i,],TN[i,])
    prob[i]= temp%*%as.matrix(gl$w)
  }
  -sum(log(prob))
}

factloglik.beta<-function(param,TP,FN,FP,TN,
                          gl,mgrid,qcond1,tau2par1,qcond2,tau2par2)
{ K=ncol(TP)
  K2=2*K
  p=param[1:K2]
  g=param[(K2+1):(2*K2)]
  tau=param[(2*K2+1):(3*K2)]
  if(any(p<=0) | any(p>=1)) return(1.e10)
  if(any(g<=0)| any(g>=1)) return(1.e10)
  if(any(tau< -0.95) | any(tau>0.95)) return(1.e10)


  th=c(sapply(tau[1:K],tau2par1),sapply(tau[(K+1):K2],tau2par2))
  
  a=p/g-p
  b=(1-p)*(1-g)/g

  nq=length(gl$nodes)
  x=array(NA,c(nq,nq,K2))
  for(j in 1:K)
  { u=qcond1(mgrid$y,mgrid$x,th[j]) 
    u[u>1]=0.9999
    x[,,j]=qbeta(u,a[j],b[j])
  }
  for(j in (K+1):K2)
  { u=qcond2(mgrid$y,mgrid$x,th[j])
    u[u>1]=0.9999
    x[,,j]=qbeta(u,a[j],b[j])
  }
  N=nrow(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=wbinomprod(x,gl,TP[i,],FN[i,],FP[i,],TN[i,])
    prob[i]= temp%*%as.matrix(gl$w)
  }
  -sum(log(prob))
}

FactorCopulaREMADA.beta=function(TP,FN,FP,TN,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2)
{ SE=TP/(TP+FN)
  SP=TN/(TN+FP)
  z=cbind(SE,SP)
  p<-apply(z,2,mean)
  si2<-apply(z,2,var)
  g=si2/p/(1-p)
  mvn.taus=bvn.cpar2tau(factanal(covmat = cov(z),factors = 1)$l[,1])
  mvn.taus[mvn.taus> 0.8]=0.8
  mvn.taus[mvn.taus< -0.8]=-0.8
  inipar=c(p,g,mvn.taus)
  est=nlm(factloglik.beta,inipar,TP,FN,FP,TN,gl,mgrid,
          qcond1,tau2par1,qcond2,tau2par2,hessian=TRUE,iterlim=1000)
  est
}



FactorCopulaREMADA.norm=function(TP,FN,FP,TN,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2)
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  SE=rTP/(rTP+rFN)
  SP=rTN/(rTN+rFP)
  z=cbind(SE,SP)
  logitz=log(z/(1-z))
  p=apply(z,2,mean)
  si<-sqrt(apply(logitz,2,var))
  temp=factanal(covmat = cov(logitz),factors = 1)
  cpar=temp$l[,1]
  mvn.taus=bvn.cpar2tau(cpar)
  mvn.taus[mvn.taus> 0.8]=0.8
  mvn.taus[mvn.taus< -0.8]=-0.8
  inipar=c(p,si,mvn.taus)
  est=nlm(factloglik.norm,inipar,TP,FN,FP,TN,gl,mgrid,
          qcond1,tau2par1,qcond2,tau2par2,hessian=TRUE,iterlim=1000)
  est
}


bvn.cpar2tau=function(rho) 
{ (2/pi) * asin(rho)
}
