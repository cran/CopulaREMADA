imperfect.fivevariateVineCopulaREMADA.loglik.norm.comprehensive<-function(param,
y11,y10,y01,y00,gl,mgrid,qcond,tau2par)
{ p=param[1:5]
  si=param[6:10]
  tau=param[11:14]
  if(any(tau< -0.95 | tau>=0.95)) return(1.e10)
  p=expit(p)
  si=exp(si)
  
  p1=mgrid$x1
  p2=mgrid$x2
  p3=mgrid$x3
  p4=mgrid$x4
  p5=mgrid$x5
 
  th=tau2par(tau)
  u1=p1
  u2=qcond(p2,u1,th[1])
  u3=qcond(p3,u2,th[2])
  u4=qcond(p4,u3,th[3])
  u5=qcond(p5,u4,th[4])
  
  
  mu=log(p/(1-p))
  
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  x3=qnorm(u3,mu[3],si[3])
  x4=qnorm(u4,mu[4],si[4])
  x5=qnorm(u5,mu[5],si[5])
  
  t1=exp(x1)
  t2=exp(x2)
  t3=exp(x3)
  t4=exp(x4)
  t5=exp(x5)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  x3=t3/(1+t3)
  x4=t4/(1+t4)
  x5=t5/(1+t5)
  lst=list(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5)
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]=tensor(tensor(tensor(tensor(temp,gl$w,
    5,1),gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
  }
  -sum(log(prob[prob!=0]))
  #out
}




imperfect.fivevariateVineCopulaREMADA.norm.comprehensive=function(y11,y10,y01,y00,
                               gl,mgrid,qcond,tau2par,start)
{ start=c(logit(start[1:5]),log(start[6:10]),start[11:14])
  est=nlm(imperfect.fivevariateVineCopulaREMADA.loglik.norm.comprehensive,start,y11,y10,y01,y00,
          gl,mgrid,qcond,tau2par,hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),exp(est$e[6:10]),est$e[11:14]),
       SE=c(se.p(est$e[1:5],se[1:5]),se.si(est$e[6:10],se[6:10]),se[11:14]))
}


imperfect.fivevariateVineCopulaREMADA.loglik.beta.comprehensive<-function(param,
                                      y11,y10,y01,y00,
                                     gl,mgrid,qcond,tau2par)
{ p=param[1:5]
  g=param[6:10]
  tau=param[11:14]
  if(any(tau< -0.95 | tau>=0.95)) return(1.e10)
  p=expit(p)
  g=expit(g)
  
  p1=mgrid$x1
  p2=mgrid$x2
  p3=mgrid$x3
  p4=mgrid$x4
  p5=mgrid$x5
  
  th=sapply(tau,tau2par)
  u1=p1
  u2=qcond(p2,u1,th[1])
  u3=qcond(p3,u2,th[2])
  u4=qcond(p4,u3,th[3])
  u5=qcond(p5,u4,th[4])
  
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  x4=qbeta(u4,a[4],b[4]) 
  x5=qbeta(u5,a[5],b[5]) 
  lst=list(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5)
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]=tensor(tensor(tensor(tensor(temp,gl$w,
                                     5,1),gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
    #if(prob!=0) -log(prob) else 0
    }
  -sum(log(prob[prob!=0]))
  #out
}



imperfect.fivevariateVineCopulaREMADA.beta.comprehensive=function(y11,y10,y01,y00,
                                          gl,mgrid,qcond,tau2par,start)
{ start=c(logit(start[1:10]),start[11:14])
  est=nlm(imperfect.fivevariateVineCopulaREMADA.loglik.beta.comprehensive,
          start,y11,y10,y01,y00,gl,mgrid,qcond,tau2par,
          hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:10]),est$e[11:14]),
       SE=c(se.p(est$e[1:10],se[1:10]),se[11:14]))
}
