rgammaShifted=function (n,shape,scale,thres) 
{ rgamma(n, shape, 1/scale) + thres }


sim1fact=function (n, parobj, qcond1, qcond2) 
{ vv = runif(n)
  d = length(parobj)
  dhalf=d/2
  yy = matrix(0, n, d)
  for (j in 1:dhalf) 
  { qq = runif(n)
    th = parobj[j]
    yy[, j] = qcond1(qq, vv, th)
  }
  for (j in (dhalf+1):d) 
  { qq = runif(n)
    th = parobj[j]
    yy[, j] = qcond2(qq, vv, th)
  }
  yy
}

rFactorCopulaREMADA.norm=function(N,p,si,taus,qcond1,tau2par1,
                                  qcond2,tau2par2)
{ K2=length(taus)
  K=K2/2
  theta=c(sapply(taus[1:K],tau2par1),sapply(taus[(K+1):K2],tau2par2))
  u=sim1fact(N,theta,qcond1,qcond2)
  n = round(rgammaShifted(N, shape = 1.2, scale = 100, thres = 30))
  n1=rbinom(N,size=n,prob=0.4)
  n2=n-n1
  mu=log(p/(1-p))
  uu=u
  for(j in 1:K2)
  {  x=qnorm(u[,j],mu[j],si[j])
     t=exp(x) 
     uu[,j]=t/(1+t)
  }

  TP=TN=matrix(NA,N,K)
  for(j in 1:K)
  { TP[,j]=rbinom(N,size=n1,prob=uu[,j])
    TN[,j]=rbinom(N,size=n2,prob=uu[,(K+j)])
  }
  FN=n1-TP
  FP=n2-TN
  list(TP=TP,TN=TN,FN=FN,FP=FP)
}

  
rFactorCopulaREMADA.beta=function(N,p,g,taus,qcond1,tau2par1,
                                  qcond2,tau2par2)
{ K2=length(taus)
  K=K2/2
  theta=c(sapply(taus[1:K],tau2par1),sapply(taus[(K+1):K2],tau2par2))
  u=sim1fact(N,theta,qcond1,qcond2)
  n = round(rgammaShifted(N, shape = 1.2, scale = 100, thres = 30))
  n1=rbinom(N,size=n,prob=0.4)
  n2=n-n1
  
  a=p/g-p
  b=(1-p)*(1-g)/g
  
  
  uu=u
  for(j in 1:K2)
  {  uu[,j]=qbeta(u[,j],a[j],b[j]) }
  
  TP=TN=matrix(NA,N,K)
  for(j in 1:K)
  { TP[,j]=rbinom(N,size=n1,prob=uu[,j])
  TN[,j]=rbinom(N,size=n2,prob=uu[,(K+j)])
  }
  FN=n1-TP
  FP=n2-TN
  list(TP=TP,TN=TN,FN=FN,FP=FP)
}
  
  
  


