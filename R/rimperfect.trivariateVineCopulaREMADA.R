# D-vine simulation as in Joe (2011, Chapter 7, Vine Copula Handbook),
dvine3dsim=function(nsim,tau,qcond1,tau2par1,qcond2,tau2par2)
{ p = matrix(runif(nsim* 3), nsim, 3)
  th=c(tau2par1(tau[1]),tau2par2(tau[2]))
  u1=p[,1]
  u2=qcond1(p[,2],u1,th[1])
  u3=qcond2(p[,3],u2,th[2])
  cbind(u1,u2,u3)
}






rimperfect.trivariateVineCopulaREMADA.norm=function(N,p,si,taus,select.random,qcond1,tau2par1,
                                                    qcond2,tau2par2)
{ dat=dvine3dsim(N,taus,qcond1,tau2par1,qcond2,tau2par2)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]
  
  mu=log(p[select.random]/(1-p[select.random])) 
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  x3=qnorm(u3,mu[3],si[3])
  t1=exp(x1)
  t2=exp(x2)
  t3=exp(x3)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  x3=t3/(1+t3)
  
  x = matrix(NA,3,N)
  x[1,]=x1
  x[2,]=x2
  x[3,]=x3
  
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  
  p11=lst$x1*lst$x2*lst$x3+(1-lst$x1)*(1-lst$x4)*(1-lst$x5)
  p10=lst$x1*lst$x2*(1-lst$x3)+(1-lst$x1)*(1-lst$x4)*lst$x5
  p01=lst$x1*(1-lst$x2)*lst$x3+(1-lst$x1)*lst$x4*(1-lst$x5)
  p00=1-p11-p10-p01
  n=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  prob=cbind(p11,p10,p01,p00)
  out=rmultinomial(N, size = n, prob = prob)
  list(y11=out[,1],y10=out[,2],y01=out[,3],y00=out[,4])
  
}

rimperfect.trivariateVineCopulaREMADA.beta=function(N,p,g,taus,select.random,qcond1,
tau2par1,qcond2,tau2par2)
{ dat=dvine3dsim(N,taus,qcond1,tau2par1,qcond2,tau2par2)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]



  a=p[select.random]/g-p[select.random]
  b=(1-p[select.random])*(1-g)/g

  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  x = matrix(NA,3,N)
  x[1,]=x1
  x[2,]=x2
  x[3,]=x3
  
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  p11=lst$x1*lst$x2*lst$x3+(1-lst$x1)*(1-lst$x4)*(1-lst$x5)
  p10=lst$x1*lst$x2*(1-lst$x3)+(1-lst$x1)*(1-lst$x4)*lst$x5
  p01=lst$x1*(1-lst$x2)*lst$x3+(1-lst$x1)*lst$x4*(1-lst$x5)
  p00=1-p11-p10-p01
  n=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  prob=cbind(p11,p10,p01,p00)
  out=rmultinomial(N, size = n, prob = prob)
  list(y11=out[,1],y10=out[,2],y01=out[,3],y00=out[,4])
}

