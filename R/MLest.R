tribinomprod=function(x1,x2,x3,TP,FN,FP,TN,NEP,NEN)
{ n1=TP+FN
  n2=TN+FP
  n3=n1+n2+NEP+NEN
  # 12, 13, 23|1
  f1=dbinom(TP,size=n1,prob=x1)
  f2=dbinom(TN,size=n2,prob=x2)
  f3=dbinom(TP+FN+NEP,size=n3,prob=x3)
  f1*f2*f3    
}


###################################################
vineloglik.norm<-function(param,TP,FN,FP,TN,NEP,NEN,gl,mgrid,
                          qcondcop12,qcondcop13,qcondcop23,
                          tau2par12,tau2par13,tau2par23)
{ p=param[1:3]
  si=param[4:6]
  tau12=param[7]
  tau13=param[8]
  tau23=param[9]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(si[1]<=0 | si[2]<=0 | si[3]<=0) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau13< -0.95 | tau13>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  th23=tau2par23(tau23)
  w1=mgrid$x
  w2=mgrid$y
  w3=mgrid$z
  u1=w1
  u2=qcondcop12(mgrid$y,mgrid$x,th12)
  t=qcondcop23(w3,w2,th23)
  u3=qcondcop13(t,w1,th13)
  mu=log(p/(1-p))
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  x3=qnorm(u3,mu[3],si[3])
  t1=exp(x1)
  t2=exp(x2)
  t3=exp(x3)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  x3=t3/(1+t3)
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tribinomprod(x1,x2,x3,TP[i],FN[i],FP[i],TN[i],NEP[i],NEN[i])
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  -sum(log(prob))
}




VineCopulaREMADA.norm=function(TP,FN,FP,TN,gl,mgrid,
                               qcondcop12,qcondcop13,qcondcop23,
                               tau2par12,tau2par13,tau2par23,
                               NEP=rep(0,length(TP)),NEN=rep(0,length(TP)))
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  rNEP=NEP + 0.5*(NEP==0)
  rNEN=NEN + 0.5*(NEN==0)
  SE=rTP/(rTP+rFN)
  SP=rTN/(rTN+rFP)
  PR=(rTP+rFN+rNEP)/(rTP+rFN+rTN+rFP+rNEP+rNEN)
  z=cbind(SE,SP,PR)
  logitz=log(z/(1-z))
  p=apply(z,2,mean)
  si<-sqrt(apply(logitz,2,var))
  stau=cor(logitz,method="kendall")
  inipar=c(p,si,stau[1,2],stau[1,3],-0.1)
  est=nlm(vineloglik.norm,inipar,TP,FN,FP,TN,NEP,NEN,gl,mgrid,
          qcondcop12,qcondcop23,qcondcop13,
          tau2par12,tau2par23,tau2par13,hessian=T)
  est
}

###################################################
vineloglik.beta<-function(param,TP,FN,FP,TN,NEP,NEN,gl,mgrid,
                          qcondcop12,qcondcop13,qcondcop23,
                          tau2par12,tau2par13,tau2par23)
{ p=param[1:3]
  g=param[4:6]
  tau12=param[7]
  tau13=param[8]
  tau23=param[9]
  if(p[1]<= 0 | p[1]>=1) return(1.e10)
  if(g[1]<= 0 | g[1]>=1) return(1.e10)
  if(p[2]<= 0 | p[2]>=1) return(1.e10)
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(g[3]<= 0 | g[3]>=1) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau13< -0.95 | tau13>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  th23=tau2par23(tau23)
  w1=mgrid$x
  w2=mgrid$y
  w3=mgrid$z
  u1=w1
  u2=qcondcop12(mgrid$y,mgrid$x,th12)
  t=qcondcop23(w3,w2,th23)
  u3=qcondcop13(t,w1,th13)
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tribinomprod(x1,x2,x3,TP[i],FN[i],FP[i],TN[i],NEP[i],NEN[i])
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  -sum(log(prob))
}




VineCopulaREMADA.beta=function(TP,FN,FP,TN,gl,mgrid,
                               qcondcop12,qcondcop13,qcondcop23,
                               tau2par12,tau2par13,tau2par23,
                               NEP=rep(0,length(TP)),NEN=rep(0,length(TP)))
{ temp1=TP+FN
  temp2=TN+FP
  z=cbind(TP/temp1,TN/temp2,(temp1+NEP)/(temp1+temp2+NEP+NEN))
  p<-apply(z,2,mean)
  si2<-apply(z,2,var)
  g=si2/p/(1-p)
  stau=cor(z,method="kendall")
  inipar=c(p,g,stau[1,2],stau[1,3],-0.1)
  est=nlm(vineloglik.beta,inipar,TP,FN,FP,TN,NEP,NEN,gl,mgrid,
          qcondcop12,qcondcop23,qcondcop13,
          tau2par12,tau2par23,tau2par13,hessian=T)
  est
}



