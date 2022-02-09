doublelik.norm=function(param,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,gl,mgrid1,mgrid2,
                        qcondcop12,qcondcop13,
                        tau2par12,tau2par13,
                        qcond,tau2par)
{ p=param[1:3]
  si=param[4:6]
  tau12=param[7]
  tau13=param[8]
  tau=param[9]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(si[1]<=0 | si[2]<=0 | si[3]<=0) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau13< -0.95 | tau13>=0.95) return(1.e10)
  if(tau< -0.95 | tau>=0.95) return(1.e10)
  N=length(TP1)
  lik.cohort=tvineloglik.norm(param[-9],TP1,FN1,FP1,TN1,
  rep(0,N),rep(0,N),gl,mgrid1,
  qcondcop12,qcondcop13,tau2par12,tau2par13)
  lik.caseControl=loglik.norm(c(p[1],p[2],si[1],si[2],tau),
                              TP2,FN2,FP2,TN2,gl,mgrid2,qcond,tau2par)
  out=lik.caseControl + lik.cohort
  out
}


doublelik.beta=function(param,TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,gl,mgrid1,mgrid2,
                        qcondcop12,qcondcop13,
                        tau2par12,tau2par13,qcond,tau2par)
{ p=param[1:3] 
  g=param[4:6]
  tau12=param[7]
  tau13=param[8]
  tau=param[9]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(g[1]<=0 | g[1]>=1) return(1.e10)
  if(g[2]<=0 | g[2]>=1) return(1.e10)
  if(g[3]<= 0 | g[3]>=1) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau13< -0.95 | tau13>=0.95) return(1.e10)
  if(tau< -0.95 | tau>=0.95) return(1.e10)
  N=length(TP1)
  lik.cohort=tvineloglik.beta(param[-9],TP1,FN1,FP1,TN1,
  rep(0,N),rep(0,N),gl,mgrid1,qcondcop12,qcondcop13,
  tau2par12,tau2par13)
  lik.caseControl=loglik.beta(c(p[1],p[2],g[1],g[2],tau),
  TP2,FN2,FP2,TN2,gl,mgrid2,qcond,tau2par)
  out= lik.caseControl + lik.cohort
  out
}

hybridCopulaREMADA.norm=function(TP,FN,FP,TN,type,gl,mgrid1,mgrid2,
               qcondcop12,qcondcop13,
               tau2par12,tau2par13,qcond,tau2par)
{ TP1=TP[type==1]
  TN1=TN[type==1]
  FP1=FP[type==1]
  FN1=FN[type==1]
  TP2=TP[type==2]
  TN2=TN[type==2]
  FP2=FP[type==2]
  FN2=FN[type==2]
  # get initials 
  rTP1=TP1 + 0.5*(TP1==0)
  rFN1=FN1 + 0.5*(FN1==0)
  rFP1=FP1 + 0.5*(FP1==0)
  rTN1=TN1 + 0.5*(TN1==0)
  SE1=rTP1/(rTP1+rFN1)
  SP1=rTN1/(rTN1+rFP1)
  PR1=(rTP1+rFN1)/(rTP1+rFN1+rTN1+rFP1)
  z1=cbind(SE1,SP1,PR1)
  logitz1=log(z1/(1-z1))
  p1=apply(z1,2,mean)
  si1<-sqrt(apply(logitz1,2,var))
  stau1=cor(logitz1,method="kendall")
  inipar1=c(p1,si1,stau1[1,2],stau1[1,3])
  rTP2=TP2 + 0.5*(TP2==0)
  rFN2=FN2 + 0.5*(FN2==0)
  rFP2=FP2 + 0.5*(FP2==0)
  rTN2=TN2 + 0.5*(TN2==0)
  SE2=rTP2/(rTP2+rFN2)
  SP2=rTN2/(rTN2+rFP2)
  z2=cbind(SE2,SP2)
  logitz2=log(z2/(1-z2))
  p2=apply(z2,2,mean)
  si2<-sqrt(apply(logitz2,2,var))
  stau2=cor(logitz2,method="kendall")[1,2]
  inipar2=c(p2,si2,stau2)
  temp=(inipar1[c(1,2,4,5)]+inipar2[1:4])/2
  inipar=c(temp[1:2],inipar1[3],temp[3:4],inipar1[6],stau1[1,2],stau1[1,3],stau2)
  est=nlm(doublelik.norm,inipar,
          TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,gl,mgrid1,mgrid2,
          qcondcop12,qcondcop13,tau2par12,tau2par13,
          qcond,tau2par,print.level=1,hessian=T) 
  est
}



hybridCopulaREMADA.beta=function(TP,FN,FP,TN,type,gl,mgrid1,mgrid2,
                    qcondcop12,qcondcop13,
                    tau2par12,tau2par13,qcond,tau2par)
{ TP1=TP[type==1]
  TN1=TN[type==1]
  FP1=FP[type==1]
  FN1=FN[type==1]
  TP2=TP[type==2]
  TN2=TN[type==2]
  FP2=FP[type==2]
  FN2=FN[type==2]
  # get initials
  rTP1=TP1 + 0.5*(TP1==0)
  rFN1=FN1 + 0.5*(FN1==0)
  rFP1=FP1 + 0.5*(FP1==0)
  rTN1=TN1 + 0.5*(TN1==0)
  SE1=rTP1/(rTP1+rFN1)
  SP1=rTN1/(rTN1+rFP1)
  PR1=(rTP1+rFN1)/(rTP1+rFN1+rTN1+rFP1)
  z1=cbind(SE1,SP1,PR1)
  p1<-apply(z1,2,mean)
  si2<-apply(z1,2,var)
  g1=si2/p1/(1-p1)
  stau1=cor(z1,method="kendall")
  inipar1=c(p1,g1,stau1[1,2],stau1[1,3])
  rTP2=TP2 + 0.5*(TP2==0)
  rFN2=FN2 + 0.5*(FN2==0)
  rFP2=FP2 + 0.5*(FP2==0)
  rTN2=TN2 + 0.5*(TN2==0)
  SE2=rTP2/(rTP2+rFN2)
  SP2=rTN2/(rTN2+rFP2)
  z2=cbind(SE2,SP2)
  p2<-apply(z2,2,mean)
  si2<-apply(z2,2,var)
  g2=si2/p2/(1-p2)
  stau2=cor(z2,method="kendall")[1,2]
  inipar2=c(p2,g2,stau2)
  temp=(inipar1[c(1,2,4,5)]+inipar2[1:4])/2
  inipar=c(temp[1:2],inipar1[3],temp[3:4],inipar1[6],stau1[1,2],stau1[1,3],stau2)
  est=nlm(doublelik.beta,inipar,
        TP1,FN1,FP1,TN1,TP2,FN2,FP2,TN2,gl,mgrid1,mgrid2,
        qcondcop12,qcondcop13,tau2par12,tau2par13,
        qcond,tau2par,print.level=2,hessian=T) 
  est
}









