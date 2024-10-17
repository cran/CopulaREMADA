imperfect.trivariateVineCopulaREMADA.loglik.norm.comprehensive<-function(param,
y11,y10,y01,y00,gl,mgrid,qcond,tau2par,select.random)
{ p=param[1:5]
  si=param[6:8]
  tau=param[9:10]
  if(any(tau>0.95) | any(tau< -0.95)) return(1.e10)
  p=expit(p)
  if(p[3]+p[5]-1<0) return(1e10)
  si=exp(si)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  
  th=sapply(tau,tau2par)
  u1=p1
  u2=qcond(p2,u1,th[1])
  u3=qcond(p3,u2,th[2])

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
  x = array(NA, dim = c(3,dim(x1)))
  x[1,,,]=x1
  x[2,,,]=x2
  x[3,,,]=x3
  
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
  }
  #out
  -sum(log(prob[prob!=0]))
}

imperfect.trivariateVineCopulaREMADA.norm.comprehensive=function(y11,y10,y01,y00,
                               gl,mgrid,qcond,tau2par,select.random,start)
{ start=c(logit(start[1:5]),log(start[6:8]),start[9:10])
  est=nlm(imperfect.trivariateVineCopulaREMADA.loglik.norm.comprehensive,start,y11,y10,y01,y00,
          gl,mgrid,qcond,tau2par,select.random,
          hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),exp(est$e[6:8]),est$e[9:10]),
     SE=c(se.p(est$e[1:5],se[1:5]),se.si(est$e[6:8],se[6:8]),se[9:10]))
}


imperfect.trivariateVineCopulaREMADA.loglik.beta.comprehensive<-function(param,
                                      y11,y10,y01,y00,
                                     gl,mgrid,qcond,tau2par,select.random)
{ p=param[1:5]
  g=param[6:8]
  tau=param[9:10]
  if(any(tau>0.95) | any(tau< -0.95)) return(1.e10)
  
  p=expit(p)
  g=expit(g)
  
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  
  th=sapply(tau,tau2par)
  u1=p1
  u2=qcond(p2,u1,th[1])
  u3=qcond(p3,u2,th[2])
  
  a=p[select.random]/g-p[select.random]
  b=(1-p[select.random])*(1-g)/g
  
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  x = array(NA, dim = c(3,dim(x1)))
  x[1,,,]=x1
  x[2,,,]=x2
  x[3,,,]=x3
  
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
  }
  #out
  -sum(log(prob[prob!=0]))
}



imperfect.trivariateVineCopulaREMADA.beta.comprehensive=function(y11,y10,y01,y00,
                                          gl,mgrid,qcond,tau2par,select.random,start)
{ start=c(logit(start[1:8]),start[9:10])
  est=nlm(imperfect.trivariateVineCopulaREMADA.loglik.beta.comprehensive,start,
          y11,y10,y01,y00,gl,mgrid,qcond,tau2par,select.random,
          hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),expit(est$e[6:8]),est$e[9:10]),
     SE=c(se.p(est$e[1:5],se[1:5]),se.p(est$e[6:8],se[6:8]),se[9:10]))
}


imperfect.trivariateVineCopulaREMADA.loglik.norm<-function(param,
  y11,y10,y01,y00,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random)
{ p=param[1:5]
  si=param[6:8]
  tau=param[9:10]
  if(tau[1]>0.95 | tau[1]< 0) return(1.e10)
  if(tau[2]>0 | tau[2]< -0.95) return(1.e10)
  p=expit(p)
  if(p[3]+p[5]-1<0) return(1e10)
  si=exp(si)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z


  th=c(tau2par1(tau[1]),tau2par2(tau[2]))
  u1=p1
  u2=qcond1(p2,u1,th[1])
  u3=qcond2(p3,u2,th[2])

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
  x = array(NA, dim = c(3,dim(x1)))
  x[1,,,]=x1
  x[2,,,]=x2
  x[3,,,]=x3
  
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
  }
  #out
  -sum(log(prob[prob!=0]))
}




imperfect.trivariateVineCopulaREMADA.norm=function(y11,y10,y01,y00,
                gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random,start)
{ start=c(logit(start[1:5]),log(start[6:8]),start[9:10])
  est=nlm(imperfect.trivariateVineCopulaREMADA.loglik.norm,start,y11,y10,y01,y00,
          gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,
          select.random,hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),exp(est$e[6:8]),est$e[9:10]),
       SE=c(se.p(est$e[1:5],se[1:5]),se.si(est$e[6:8],se[6:8]),se[9:10]))
}


imperfect.trivariateVineCopulaREMADA.loglik.beta<-function(param,
y11,y10,y01,y00,gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random)
{ p=param[1:5]
  g=param[6:8]
  tau=param[9:10]
  if(tau[1]>0.95 | tau[1]< 0) return(1.e10)
  if(tau[2]>0 | tau[2]< -0.95) return(1.e10)

  p=expit(p)
  g=expit(g)

  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z

  th=c(tau2par1(tau[1]),tau2par2(tau[2]))
  u1=p1
  u2=qcond1(p2,u1,th[1])
  u3=qcond2(p3,u2,th[2])

  a=p[select.random]/g-p[select.random]
  b=(1-p[select.random])*(1-g)/g

  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  
  x = array(NA, dim = c(3,dim(x1)))
  x[1,,,]=x1
  x[2,,,]=x2
  x[3,,,]=x3
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
    #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
  }
  #out
  -sum(log(prob[prob!=0]))
}



imperfect.trivariateVineCopulaREMADA.beta=function(y11,y10,y01,y00,
                                                   gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,select.random,start)
{ start=c(logit(start[1:8]),start[9:10])
  est=nlm(imperfect.trivariateVineCopulaREMADA.loglik.beta,start,y11,y10,y01,y00,
          gl,mgrid,qcond1,tau2par1,qcond2,tau2par2,
          select.random,hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),expit(est$e[6:8]),est$e[9:10]),
     SE=c(se.p(est$e[1:5],se[1:5]),se.p(est$e[6:8],se[6:8]),se[9:10]))
}

