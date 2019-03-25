# D-vine simulation as in Joe (2011, Chapter 7, Vine Copula Handbook),

dtrvinesim=function(nsim,trparam,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
{ tau12=trparam[1]
  tau23=trparam[2]
  tau34=trparam[3]
  #tau13.2=param[4]
  #tau24.3=param[5]
  #tau14.23=param[6]

  p = matrix(runif(nsim * 4), nsim, 4)

  th=matrix(0,4,4)
  th[1,2]=tau2par1(tau12)
  th[1,3]=tau2par2(tau23)
  th[1,4]=tau2par1(tau34)
  #th[2,3]=tau2par.bvn(tau13.2)
  #th[2,4]=tau2par.bvn(tau24.3)
  #th[3,4]=tau2par.bvn(tau14.23)

  u1=p[,1]
  q11=p[,1]
  q22=p[,2]
  u2=qcond1(p[,2],p[,1],th[1,2])
  q12=u2
  v12=pcond1(u1,u2,th[1,2])

  q33=p[,3]

  q23=q33 #qcond(q33,v12,th[2,3]) because qcond(u,v,indep)=u
  q13=qcond2(q23,u2,th[1,3])
  u3=q13
  v13=pcond2(u2,u3,th[1,3])

  v23=v12 #pcond(v12,q23,th[2,3]) because pcond(u,v,indep)=u
  
  q44=p[,4]

  q34=q44 #qcond(q44,v23,th[3,4]) because qcond(u,v,indep)=u
 
  q24=q34 #qcond(q34,v13,th[2,4]) because qcond(u,v,indep)=u
  q14=qcond1(q24,u3,th[1,4])
  u4=q14
  cbind(u1,u2,u3,u4)
}







rmultinomVineCopulaREMADA.beta=function(N,p,g,taus,qcond1,
pcond1,tau2par1,qcond2,
pcond2,tau2par2)
{ dat=dtrvinesim(N,taus,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]
  u4=dat[,4]
  
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
  
  
  n1=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n2=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  
  
  
  prob1=cbind(p1,p3,1-p1-p3)
  prob2=cbind(p2,p4,1-p2-p4)
  out1=rmultinomial(N, size = n1, prob = prob1)
  out2=rmultinomial(N, size = n2, prob = prob2)
  cbind(out1,out2)
}


rmultinomVineCopulaREMADA.norm=function(N,p,si,taus,qcond1,
                                    pcond1,tau2par1,qcond2,
                                    pcond2,tau2par2)
{ dat=dtrvinesim(N,taus,qcond1,pcond1,tau2par1,qcond2,pcond2,tau2par2)
u1=dat[,1]
u2=dat[,2]
u3=dat[,3]
u4=dat[,4]

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


n1=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
n2=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))



prob1=cbind(p1,p3,1-p1-p3)
prob2=cbind(p2,p4,1-p2-p4)
out1=rmultinomial(N, size = n1, prob = prob1)
out2=rmultinomial(N, size = n2, prob = prob2)
cbind(out1,out2)
}

