dvinesim=function(nsim,param,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
{ tau12=param[1]
  tau23=param[2]
  tau34=param[3]
  tau13.2=param[4]
  tau24.3=param[5]
  tau14.23=param[6]
  
  p = matrix(runif(nsim * 4), nsim, 4)
  
  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par.bvn(tau23)
  th[1,4]=tau2par2(tau34)
  th[2,3]=tau2par.bvn(tau13.2)
  th[2,4]=tau2par.bvn(tau24.3)
  th[3,4]=tau2par.bvn(tau14.23)
  
  u1=p[,1]
  q11=p[,1]
  q22=p[,2]
  u2=qcond1(p[,2],p[,1],th[1,2])
  q12=u2
  v12=pcond1(u1,u2,th[1,2])
  
  q33=p[,3]
  
  q23=qcondbvn(q33,v12,th[2,3]) 
  q13=qcondbvn(q23,u2,th[1,3])
  u3=q13
  v13=pcondbvn(u2,u3,th[1,3])
  
  v23=pcondbvn(v12,q23,th[2,3]) 
  
  q44=p[,4]
  
  q34=qcondbvn(q44,v23,th[3,4])
  
  q24=qcondbvn(q34,v13,th[2,4])
  q14=qcond2(q24,u3,th[1,4])
  u4=q14
  cbind(u1,u2,u3,u4)
}


