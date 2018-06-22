tetrabinomprod=function(x1,x2,x3,x4,TP1,FN1,FP1,TN1,
                        TP2,FN2,FP2,TN2)
{ n11=TP1+FN1
  n21=TN1+FP1
  n12=TP2+FN2
  n22=TN2+FP2
  f1=dbinom(TP1,size=n11,prob=x1)
  f2=dbinom(TN1,size=n21,prob=x2)
  f3=dbinom(TP2,size=n12,prob=x3)
  f4=dbinom(TN2,size=n22,prob=x4)
  f1*f2*f3*f4    
}



###################################################


quadvineloglik.norm<-function(param,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                              gl,mgrid,qcond1,pcond1,tau2par1,
                            qcond2,pcond2,tau2par2)
{ p=param[1:4]
  si=param[5:8]
  tau12=param[9]
  tau23=param[10]
  tau34=param[11]
  tau13.2=param[12]
  tau24.3=param[13]
  tau14.23=param[14]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(p[4]<= 0 | p[4]>=1) return(1.e10)

  if(si[1]<=0 | si[2]<=0 | si[3]<=0 | si[4]<=0) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau34< -0.95 | tau34>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  if(tau13.2< -0.95 | tau13.2>=0.95) return(1.e10)
  if(tau24.3< -0.95 | tau24.3>=0.95) return(1.e10)
  if(tau14.23< -0.95 | tau14.23>=0.95) return(1.e10)
  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par.bvn(tau23)
  th[1,4]=tau2par2(tau34)
  th[2,3]=tau2par.bvn(tau13.2)
  th[2,4]=tau2par.bvn(tau24.3)
  th[3,4]=tau2par.bvn(tau14.23)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  p4=mgrid$w

  u1=p1
  q11=p1
  q22=p2
  u2=qcond1(p2,p1,th[1,2])
  q12=u2
  v12=pcond1(u1,u2,th[1,2])
  
  q33=p3
  q23=qcondbvn(q33,v12,th[2,3]) 
  q13=qcondbvn(q23,u2,th[1,3])
  u3=q13
  v13=pcondbvn(u2,u3,th[1,3])

  v23=pcondbvn(v12,q23,th[2,3]) 
  q44=p4
  q34=qcondbvn(q44,v23,th[3,4])
  
  q24=qcondbvn(q34,v13,th[2,4])
  q14=qcond2(q24,u3,th[1,4])
  u4=q14

  mu=log(p/(1-p))
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  x3=qnorm(u3,mu[3],si[3])
  x4=qnorm(u4,mu[4],si[4])

  t1=exp(x1)
  t2=exp(x2)
  t3=exp(x3)
  t4=exp(x4)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  x3=t3/(1+t3)
  x4=t4/(1+t4)
  N=length(TP1)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tetrabinomprod(x1,x2,x3,x4,TP1[i],FN1[i],FP1[i],TN1[i],
                        TP2[i],FN2[i],FP2[i],TN2[i])
    prob[i]= tensor(tensor(tensor(temp,gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  -sum(log(prob))
}

   


quadVineCopulaREMADA.norm=function(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                               gl,mgrid,qcond1,pcond1,tau2par1,
                               qcond2,pcond2,tau2par2)
{ rTP1=TP1 + 0.5*(TP1==0)
  rFN1=FN1 + 0.5*(FN1==0)
  rFP1=FP1 + 0.5*(FP1==0)
  rTN1=TN1 + 0.5*(TN1==0)
  SE1=rTP1/(rTP1+rFN1)
  SP1=rTN1/(rTN1+rFP1)
  
  rTP2=TP2 + 0.5*(TP2==0)
  rFN2=FN2 + 0.5*(FN2==0)
  rFP2=FP2 + 0.5*(FP2==0)
  rTN2=TN2 + 0.5*(TN2==0)
  SE2=rTP2/(rTP2+rFN2)
  SP2=rTN2/(rTN2+rFP2)
  
  
  
  z=cbind(SE1,SP1,SE2,SP2)
  
  logitz=log(z/(1-z))
  p=apply(z,2,mean)
  si<-sqrt(apply(logitz,2,var))
  stau=cor(logitz,method="kendall")
  inipar=c(p,si,stau[1,2],stau[2,3],stau[3,4],rep(0,3))
  est=nlm(quadvineloglik.norm,inipar,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
          gl,mgrid,qcond1,pcond1,tau2par1,
          qcond2,pcond2,tau2par2,hessian=T,iterlim=1000,print.level=2)
  est
}




###################################################
quadvineloglik.beta<-function(param,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                            gl,mgrid,qcond1,pcond1,tau2par1,
                            qcond2,pcond2,tau2par2)
{ p=param[1:4]
  g=param[5:8]
  tau12=param[9]
  tau23=param[10]
  tau34=param[11]
  tau13.2=param[12]
  tau24.3=param[13]
  tau14.23=param[14]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(p[4]<= 0 | p[4]>=1) return(1.e10)
  if(g[1]<= 0 | g[1]>=1) return(1.e10)
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  if(g[3]<= 0 | g[3]>=1) return(1.e10)
  if(g[4]<= 0 | g[4]>=1) return(1.e10)
  
  
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau34< -0.95 | tau34>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  if(tau13.2< -0.95 | tau13.2>=0.95) return(1.e10)
  if(tau24.3< -0.95 | tau24.3>=0.95) return(1.e10)
  if(tau14.23< -0.95 | tau14.23>=0.95) return(1.e10)
  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par.bvn(tau23)
  th[1,4]=tau2par2(tau34)
  th[2,3]=tau2par.bvn(tau13.2)
  th[2,4]=tau2par.bvn(tau24.3)
  th[3,4]=tau2par.bvn(tau14.23)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  p4=mgrid$w
  u1=p1
  q11=p1
  q22=p2
  u2=qcond1(p2,p1,th[1,2])
  q12=u2
  v12=pcond1(u1,u2,th[1,2])
  
  q33=p3
  q23=qcondbvn(q33,v12,th[2,3]) 
  q13=qcondbvn(q23,u2,th[1,3])
  u3=q13
  v13=pcondbvn(u2,u3,th[1,3])
  
  v23=pcondbvn(v12,q23,th[2,3]) 
  q44=p4
  q34=qcondbvn(q44,v23,th[3,4])
  
  q24=qcondbvn(q34,v13,th[2,4])
  q14=qcond2(q24,u3,th[1,4])
  u4=q14
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3])
  x4=qbeta(u4,a[4],b[4])
  
  N=length(TP1)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tetrabinomprod(x1,x2,x3,x4,TP1[i],FN1[i],FP1[i],TN1[i],
                        TP2[i],FN2[i],FP2[i],TN2[i])
    prob[i]= tensor(tensor(tensor(temp,gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  -sum(log(prob))
}




quadVineCopulaREMADA.beta=function(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                                 gl,mgrid,qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
{ rTP1=TP1 + 0.5*(TP1==0)
  rFN1=FN1 + 0.5*(FN1==0)
  rFP1=FP1 + 0.5*(FP1==0)
  rTN1=TN1 + 0.5*(TN1==0)
  SE1=rTP1/(rTP1+rFN1)
  SP1=rTN1/(rTN1+rFP1)

  rTP2=TP2 + 0.5*(TP2==0)
  rFN2=FN2 + 0.5*(FN2==0)
  rFP2=FP2 + 0.5*(FP2==0)
  rTN2=TN2 + 0.5*(TN2==0)
  SE2=rTP2/(rTP2+rFN2)
  SP2=rTN2/(rTN2+rFP2)
  
  z=cbind(SE1,SP1,SE2,SP2)
  p<-apply(z,2,mean)
  si2<-apply(z,2,var)
  g=si2/p/(1-p)
  stau=cor(z,method="kendall")
  inipar=c(p,g,stau[1,2],stau[2,3],stau[3,4],rep(0,3))
  est=nlm(quadvineloglik.beta,inipar,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
          gl,mgrid,qcond1,pcond1,tau2par1,
          qcond2,pcond2,tau2par2,hessian=T,iterlim=1000,print.level=2)
  est
}

quadvineloglik.norm.beta<-function(param,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                                 gl,mgrid,qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
{ p=param[1:4]
  si=param[5:6]
  g=param[7:8]
  tau12=param[9]
  tau23=param[10]
  tau34=param[11]
  tau13.2=param[12]
  tau24.3=param[13]
  tau14.23=param[14]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(p[4]<= 0 | p[4]>=1) return(1.e10)
  if(g[1]<= 0 | g[1]>=1) return(1.e10)
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  if(si[1]<=0 | si[2]<=0) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau34< -0.95 | tau34>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  if(tau13.2< -0.95 | tau13.2>=0.95) return(1.e10)
  if(tau24.3< -0.95 | tau24.3>=0.95) return(1.e10)
  if(tau14.23< -0.95 | tau14.23>=0.95) return(1.e10)
  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par.bvn(tau23)
  th[1,4]=tau2par2(tau34)
  th[2,3]=tau2par.bvn(tau13.2)
  th[2,4]=tau2par.bvn(tau24.3)
  th[3,4]=tau2par.bvn(tau14.23)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  p4=mgrid$w

  u1=p1
  q11=p1
  q22=p2
  u2=qcond1(p2,p1,th[1,2])
  q12=u2
  v12=pcond1(u1,u2,th[1,2])
  
  q33=p3
  q23=qcondbvn(q33,v12,th[2,3]) 
  q13=qcondbvn(q23,u2,th[1,3])
  u3=q13
  v13=pcondbvn(u2,u3,th[1,3])
  
  v23=pcondbvn(v12,q23,th[2,3]) 
  q44=p4
  q34=qcondbvn(q44,v23,th[3,4])
  
  q24=qcondbvn(q34,v13,th[2,4])
  q14=qcond2(q24,u3,th[1,4])
  u4=q14

  p1=p[1:2]
  mu=log(p1/(1-p1))
  p2=p[3:4]
  a=p2/g-p2
  b=(1-p2)*(1-g)/g

  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  t1=exp(x1)
  t2=exp(x2)
  x1=t1/(1+t1)
  x2=t2/(1+t2)

  x3=qbeta(u3,a[1],b[1])
  x4=qbeta(u4,a[2],b[2])

  N=length(TP1)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tetrabinomprod(x1,x2,x3,x4,TP1[i],FN1[i],FP1[i],TN1[i],
                        TP2[i],FN2[i],FP2[i],TN2[i])
    prob[i]= tensor(tensor(tensor(temp,gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  -sum(log(prob))
}

quadVineCopulaREMADA.norm.beta=function(TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
                                      gl,mgrid,qcond1,pcond1,tau2par1,
                                      qcond2,pcond2,tau2par2)
{ rTP1=TP1 + 0.5*(TP1==0)
  rFN1=FN1 + 0.5*(FN1==0)
  rFP1=FP1 + 0.5*(FP1==0)
  rTN1=TN1 + 0.5*(TN1==0)
  SE1=rTP1/(rTP1+rFN1)
  SP1=rTN1/(rTN1+rFP1)

  rTP2=TP2 + 0.5*(TP2==0)
  rFN2=FN2 + 0.5*(FN2==0)
  rFP2=FP2 + 0.5*(FP2==0)
  rTN2=TN2 + 0.5*(TN2==0)
  SE2=rTP2/(rTP2+rFN2)
  SP2=rTN2/(rTN2+rFP2)
  
  z=cbind(SE1,SP1,SE2,SP2)
  p<-apply(z,2,mean)

  si2<-apply(z,2,var)
  g=si2/p/(1-p)

  logitz=log(z/(1-z))
  si<-sqrt(apply(logitz,2,var))

  stau=cor(z,method="kendall")
  inipar=c(p,c(si[1:2],g[3:4]),stau[1,2],stau[2,3],stau[3,4],rep(0,3))
  est=nlm(quadvineloglik.norm.beta,inipar,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,
        gl,mgrid,qcond1,pcond1,tau2par1,
        qcond2,pcond2,tau2par2,hessian=T,print.level = 2,iterlim=1000)
  est
}










