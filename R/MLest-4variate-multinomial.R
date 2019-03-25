
multinomvineloglik.norm<-function(param,TP,FN,FP,TN,NEP,NEN,
gl,mgrid,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
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
  if(1-p[1]-p[3]<= 0 | 1-p[1]-p[3]>=1) return(1.e10)
  if(1-p[2]-p[4]<= 0 | 1-p[2]-p[4]>=1) return(1.e10)

  if(si[1]<=0 | si[2]<=0 | si[3]<=0 | si[4]<=0) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau34< -0.95 | tau34>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  if(tau13.2< -0.95 | tau13.2>=0.95) return(1.e10)
  if(tau24.3< -0.95 | tau24.3>=0.95) return(1.e10)
  if(tau14.23< -0.95 | tau14.23>=0.95) return(1.e10)
  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par2(tau23)
  th[1,4]=tau2par1(tau34)
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
  #j=3
  q33=p3
  #ell=2
  q23=qcondbvn(q33,v12,th[2,3]) 
  q13=qcond2(q23,u2,th[1,3])
  u3=q13
  v13=pcond2(u2,u3,th[1,3])
  #ell=2
  v23=pcondbvn(v12,q23,th[2,3]) 
  #j=4
  q44=p4
  #ell=3
  q34=qcondbvn(q44,v23,th[3,4])
  #ell=2
  q24=qcondbvn(q34,v13,th[2,4])
  q14=qcond1(q24,u3,th[1,4])
  u4=q14

  # connect TP with TN
  mu1=log(p[1]/(1-p[1]-p[3]))
  mu2=log(p[2]/(1-p[2]-p[4]))
    
    
  # connect NEP with NEN
  mu3=log(p[3]/(1-p[1]-p[3]))
  mu4=log(p[4]/(1-p[2]-p[4]))
  

  x1=qnorm(u1,mu1,si[1])
  x2=qnorm(u2,mu2,si[2])
  x3=qnorm(u3,mu3,si[3])
  x4=qnorm(u4,mu4,si[4])

  t1=1+exp(x1)+exp(x3)
  t2=1+exp(x2)+exp(x4)
  
  p1=exp(x1)/t1
  p2=exp(x2)/t2
  
  p3=exp(x3)/t1
  p4=exp(x4)/t2
  
  N=length(TP)
  s=0
  for(i in 1:N)
  { temp=multinomprod(p1,p2,p3,p4,TP[i],FN[i],FP[i],TN[i],
                      NEP[i],NEN[i])
  prob=tensor(tensor(tensor(temp,gl$w,4,1),gl$w,
                     3,1),gl$w,2,1)%*%gl$w
  s=s+log(prob)
  }
  -s
}




multinomVineCopulaREMADA.norm=function(TP,FN,FP,TN,NEP,NEN,
                               gl,mgrid,qcond1,pcond1,tau2par1,
                               qcond2,pcond2,tau2par2)
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  rNEP=NEP + 0.5*(NEP==0)
  rNEN=NEN + 0.5*(NEN==0)
  
  p1=rTP/(rTP+rFN+rNEP)
  p2=rTN/(rTN+rFP+rNEN)
  p3=rNEP/(rTP+rFN+rNEP)
  p4=rNEN/(rTN+rFP+rNEN)
  
  
  
 
  p=apply(cbind(p1,p2,p3,p4),2,mean)
  
  z1=log(p1/(1-p1-p3))
  z2=log(p2/(1-p2-p4))
  z3=log(p3/(1-p1-p3))
  z4=log(p4/(1-p2-p4))
  z=cbind(z1,z2,z3,z4)
  si<-sqrt(apply(z,2,var))
  stau=cor(z,method="kendall")
  inipar=c(p,si,stau[1,2],stau[2,3],stau[3,4],rep(0,3))
  est=nlm(multinomvineloglik.norm,inipar,TP,FN,FP,TN,NEP,NEN,
          gl,mgrid,qcond1,pcond1,tau2par1,
          qcond2,pcond2,tau2par2,hessian=T,print.level = 1,iterlim=1000)
  est
}


multinomvineloglik.beta<-function(param,TP,FN,FP,TN,NEP,NEN,
                            gl,mgrid,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
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
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  if(g[3]<= 0 | g[3]>=1) return(1.e10)
  if(g[4]<= 0 | g[4]>=1) return(1.e10)
  if(1-p[1]-p[3]<= 0 | 1-p[1]-p[3]>=1) return(1.e10)
  if(1-p[2]-p[4]<= 0 | 1-p[2]-p[4]>=1) return(1.e10)


  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau34< -0.95 | tau34>=0.95) return(1.e10)
  if(tau23< -0.95 | tau23>=0.95) return(1.e10)
  if(tau13.2< -0.95 | tau13.2>=0.95) return(1.e10)
  if(tau24.3< -0.95 | tau24.3>=0.95) return(1.e10)
  if(tau14.23< -0.95 | tau14.23>=0.95) return(1.e10)



  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par2(tau23)
  th[1,4]=tau2par1(tau34)
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
  #j=3
  q33=p3
  #ell=2
  q23=qcondbvn(q33,v12,th[2,3]) 
  q13=qcond2(q23,u2,th[1,3])
  u3=q13
  v13=pcond2(u2,u3,th[1,3])
  #ell=2
  v23=pcondbvn(v12,q23,th[2,3]) 
  #j=4
  q44=p4
  #ell=3
  q34=qcondbvn(q44,v23,th[3,4])
  #ell=2
  q24=qcondbvn(q34,v13,th[2,4])
  q14=qcond1(q24,u3,th[1,4])
  u4=q14

  p1=p[1]
  p3=p[3]/(1-p[1])

  p2=p[2]
  p4=p[4]/(1-p[2])

  a1=p1/g[1]-p1
  b1=(1-p1)*(1-g[1])/g[1]

  a2=p2/g[2]-p2
  b2=(1-p2)*(1-g[2])/g[2]

  a3=p3/g[3]-p3
  b3=(1-p3)*(1-g[3])/g[3]

  a4=p4/g[4]-p4
  b4=(1-p4)*(1-g[4])/g[4]


  x1=qbeta(u1,a1,b1)
  x2=qbeta(u2,a2,b2) 
  x3=qbeta(u3,a3,b3)
  x4=qbeta(u4,a4,b4)

  p1=x1
  p3=x3*(1-x1)

  p2=x2
  p4=x4*(1-x2)
  
  N=length(TP)
  s=0
  for(i in 1:N)
  { temp=multinomprod(p1,p2,p3,p4,TP[i],FN[i],FP[i],TN[i],
                      NEP[i],NEN[i])
  prob=tensor(tensor(tensor(temp,gl$w,4,1),gl$w,
                     3,1),gl$w,2,1)%*%gl$w
  s=s+log(prob)
  }
  -s
}



multinomVineCopulaREMADA.beta=function(TP,FN,FP,TN,NEP,NEN,
                                 gl,mgrid,qcond1,pcond1,tau2par1,
                                 qcond2,pcond2,tau2par2)
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  rNEP=NEP + 0.5*(NEP==0)
  rNEN=NEN + 0.5*(NEN==0)

  p1=rTP/(rTP+rFN+rNEP)
  p2=rTN/(rTN+rFP+rNEN)
  p3=rNEP/(rTP+rFN+rNEP)
  p4=rNEN/(rTN+rFP+rNEN)
  p=apply(cbind(p1,p2,p3,p4),2,mean)

  x1=p1
  x3=p3/(1-p1)
  x2=p2
  x4=p4/(1-p2)
  x=cbind(x1,x2,x3,x4)
  meanx <- apply(x, 2, mean)
  varx <- apply(x, 2, var)
  g = varx/meanx/(1 - meanx)

  stau = cor(x, method = "kendall")
  
  inipar = c(p, g, stau[1, 2], stau[2, 3], stau[3, 4], rep(0,3))
  est=nlm(multinomvineloglik.beta,inipar,TP,FN,FP,TN,NEP,NEN,
  gl,mgrid,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2,hessian=T,
  print.level = 1,iterlim=1000)
  est
}
