imperfect.quadrivariateVineCopulaREMADA.loglik.norm.comprehensive<-function(param,
                                                  y11,y10,y01,y00,
                                                  gl,mgrid,qcond,tau2par,select.random)
{ p=param[1:5]
  si=param[6:9]
  tau=param[10:12]
  if(any(tau< -0.95 | tau>=0.95)) return(1.e10)
  p=expit(p)
  si=exp(si)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  p4=mgrid$w
  th=sapply(tau,tau2par)
  u1=p1
  u2=qcond(p2,u1,th[1])
  u3=qcond(p3,u2,th[2])
  u4=qcond(p4,u3,th[3])
  mu=log(p[select.random]/(1-p[select.random]))
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
  x = array(NA, dim = c(4,dim(x1)))
  x[1,,,,]=x1
  x[2,,,,]=x2
  x[3,,,,]=x3
  x[4,,,,]=x4
  
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]=tensor(tensor(tensor(temp,gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
   #if(prob!=0) -log(prob) else 0
  }
  -sum(log(prob[prob!=0]))
  #out
}



imperfect.quadrivariateVineCopulaREMADA.norm.comprehensive=function(y11,y10,y01,y00,
                               gl,mgrid,qcond,tau2par,select.random,start)
{ start=c(logit(start[1:5]),log(start[6:9]),start[10:12])
  est=nlm(imperfect.quadrivariateVineCopulaREMADA.loglik.norm.comprehensive,start,y11,y10,y01,y00,
          gl,mgrid,qcond,tau2par,select.random,
          hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),exp(est$e[6:9]),est$e[10:12]),
     SE=c(se.p(est$e[1:5],se[1:5]),se.si(est$e[6:9],se[6:9]),se[10:12]))
}




imperfect.quadrivariateVineCopulaREMADA.loglik.beta.comprehensive<-function(param,
                                                                            y11,y10,y01,y00,
                                                                            gl,mgrid,qcond,tau2par,select.random)
{ p=param[1:5]
  g=param[6:9]
  tau=param[10:12]
  if(any(tau< -0.95 | tau>=0.95)) return(1.e10)
  p=expit(p)
  g=expit(g)
  p1=mgrid$x
  p2=mgrid$y
  p3=mgrid$z
  p4=mgrid$w
  th=sapply(tau,tau2par)
  u1=p1
  u2=qcond(p2,u1,th[1])
  u3=qcond(p3,u2,th[2])
  u4=qcond(p4,u3,th[3])
 
  a=p[select.random]/g-p[select.random]
  b=(1-p[select.random])*(1-g)/g
  
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3])
  x4=qbeta(u4,a[4],b[4])
  x = array(NA, dim = c(4,dim(x1)))
  x[1,,,,]=x1
  x[2,,,,]=x2
  x[3,,,,]=x3
  x[4,,,,]=x4

  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]

  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]=tensor(tensor(tensor(temp,gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
  }
  -sum(log(prob[prob!=0]))
  #out
}



imperfect.quadrivariateVineCopulaREMADA.beta.comprehensive=function(y11,y10,y01,y00,
                                                                    gl,mgrid,qcond,tau2par,select.random,start)
{ start=c(logit(start[1:9]),start[10:12])
  est=nlm(imperfect.quadrivariateVineCopulaREMADA.loglik.beta.comprehensive,start,y11,y10,y01,y00,
          gl,mgrid,qcond,tau2par,select.random,
          hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:9]),est$e[10:12]),
       SE=c(se.p(est$e[1:9],se[1:9]),se[10:12]))
}



