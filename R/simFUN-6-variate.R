# D-vine simulation as in Joe (2011, Chapter 7, Vine Copula Handbook),
dvine6dsim=function(nsim,tau,qcond,tau2par)
{ p = matrix(runif(nsim * 6), nsim, 6)
  th=sapply(tau,tau2par)
  u1=p[,1]
  u2=qcond(p[,2],u1,th[1])
  u3=qcond(p[,3],u2,th[2])
  u4=qcond(p[,4],u3,th[3])
  u5=qcond(p[,5],u4,th[4])
  u6=qcond(p[,6],u5,th[5])
  cbind(u1,u2,u3,u4,u5,u6)
}








rmultinom6dVineCopulaREMADA.beta=function(N,p,g,taus,qcond,tau2par)
{ dat=dvine6dsim(N,taus,qcond,tau2par)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]
  u4=dat[,4]
  u5=dat[,5]
  u6=dat[,6]
  
  

  p1=p[1]
  p2=p[2]/(1-p1)
  p3=p[3]/(1-p2)/(1-p1)
  p4=p[4]
  p5=p[5]/(1-p4)
  p6=p[6]/(1-p5)/(1-p4)
  
  
  
  a1=p1/g[1]-p1
  b1=(1-p1)*(1-g[1])/g[1]
  
  a2=p2/g[2]-p2
  b2=(1-p2)*(1-g[2])/g[2]
  
  a3=p3/g[3]-p3
  b3=(1-p3)*(1-g[3])/g[3]
  
  a4=p4/g[4]-p4
  b4=(1-p4)*(1-g[4])/g[4]
  
  a5=p5/g[5]-p5
  b5=(1-p5)*(1-g[5])/g[5]
  
  a6=p6/g[6]-p6
  b6=(1-p6)*(1-g[6])/g[6]
  
  x1=qbeta(u1,a1,b1)
  x2=qbeta(u2,a2,b2) 
  x3=qbeta(u3,a3,b3)
  x4=qbeta(u4,a4,b4)
  x5=qbeta(u5,a5,b5)
  x6=qbeta(u6,a6,b6)
  
  
  p1=x1
  p2=x2*(1-x1)
  p3=x3*(1-x2)*(1-x1)
  
  p4=x4
  p5=x5*(1-x4)
  p6=x6*(1-x5)*(1-x4)
  
  n1=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n2=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  
  prob1=cbind(p1,p2,p3,1-p1-p2-p3)
  prob2=cbind(p4,p5,p6,1-p4-p5-p6)
  out1=rmultinomial(N, size = n1, prob = prob1)
  out2=rmultinomial(N, size = n2, prob = prob2)
  cbind(out1,out2)
}


rmultinom6dVineCopulaREMADA.norm=function(N,p,si,taus,qcond,tau2par)
{ dat=dvine6dsim(N,taus,qcond,tau2par)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]
  u4=dat[,4]
  u5=dat[,5]
  u6=dat[,6]



  mu1=log(p[1]/(1-sum(p[1:3])))
  mu2=log(p[2]/(1-sum(p[1:3])))
  mu3=log(p[3]/(1-sum(p[1:3])))
  mu4=log(p[4]/(1-sum(p[4:6])))
  mu5=log(p[5]/(1-sum(p[4:6])))
  mu6=log(p[6]/(1-sum(p[4:6])))

  x1=qnorm(u1,mu1,si[1])
  x2=qnorm(u2,mu2,si[2])
  x3=qnorm(u3,mu3,si[3])
  x4=qnorm(u4,mu4,si[4])
  x5=qnorm(u5,mu5,si[5])
  x6=qnorm(u6,mu6,si[6])

  t1=1+exp(x1)+exp(x2)+exp(x3)
  t2=1+exp(x4)+exp(x5)+exp(x6)

  p1=exp(x1)/t1
  p2=exp(x2)/t1
  p3=exp(x3)/t1
  p4=exp(x4)/t2
  p5=exp(x5)/t2
  p6=exp(x6)/t2


  n1=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n2=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))

  prob1=cbind(p1,p2,p3,1-p1-p2-p3)
  prob2=cbind(p4,p5,p6,1-p4-p5-p6)
  out1=rmultinomial(N, size = n1, prob = prob1)
  out2=rmultinomial(N, size = n2, prob = prob2)
  cbind(out1,out2)
}


