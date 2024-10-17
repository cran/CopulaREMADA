imperfect.CopulaREMADA.loglik.norm<-function(param,y11,y10,y01,y00,gl,mgrid,qcond,
                                        tau2par,select.random)
{ p=param[1:5]
  si=param[6:7]
  tau=param[8]
  if(tau< -0.95 | tau>0.95) return(1.e10)
  p=expit(p)
  if(p[3]+p[5]-1<0) return(1e10)
  si=exp(si)
  
  u1=mgrid$x
  th=tau2par(tau)
  u2=qcond(mgrid$y,mgrid$x,th)
  mu=log(p[select.random]/(1-p[select.random]))
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  t1=exp(x1)
  t2=exp(x2)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  x = array(NA, dim = c(2,dim(x1)))
  x[1,,]=x1
  x[2,,]=x2
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= prob[i]=gl$w %*% temp %*% as.matrix(gl$w)
    #if(prob!=0) -log(prob) else 0
  }
  -sum(log(prob[prob!=0]))
}




imperfectCopulaREMADA.norm=function(y11,y10,y01,y00,
                               gl,mgrid,qcond,tau2par,select.random,start)
{ start=c(logit(start[1:5]),log(start[6:7]),start[8])
  est=nlm(imperfect.CopulaREMADA.loglik.norm,start,y11,y10,y01,y00,
          gl,mgrid,qcond,tau2par,select.random,
          hessian=T,print.level =2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),exp(est$e[6:7]),est$e[8]),
     SE=c(se.p(est$e[1:5],se[1:5]),se.si(est$e[6:7],se[6:7]),se[8]))
}


imperfect.CopulaREMADA.loglik.beta<-function(param,y11,y10,y01,y00,gl,mgrid,qcond,
                                             tau2par,select.random)
{ p=param[1:5]
  g=param[6:7]
  tau=param[8]
  if(tau< -0.95 | tau>0) return(1.e10)
  p=expit(p)
  g=expit(g)
 
  u1=mgrid$x
  th=tau2par(tau)
  u2=qcond(mgrid$y,mgrid$x,th)
  
  
  
  a=p[select.random]/g-p[select.random]
  b=(1-p[select.random])*(1-g)/g
  
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  
  x = array(NA, dim = c(2,dim(x1)))
  x[1,,]=x1
  x[2,,]=x2
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  for(jj in 1:length(select.random)) lst[[select.random[jj]]]=x[jj,,]
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= prob[i]=gl$w %*% temp %*% as.matrix(gl$w)
    #if(prob!=0) -log(prob) else 0
  }
  -sum(log(prob[prob!=0]))
  #out
}



imperfectCopulaREMADA.beta=function(y11,y10,y01,y00,gl,mgrid,qcond,
                                            tau2par,select.random,start)
{ start=c(logit(start[1:5]),logit(start[6:7]),start[8])
  est=nlm(imperfect.CopulaREMADA.loglik.beta,start,y11,y10,y01,y00,
            gl,mgrid,qcond,tau2par,select.random,
            hessian=T,print.level =2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),expit(est$e[6:7]),est$e[8]),
     SE=c(se.p(est$e[1:5],se[1:5]),se.p(est$e[6:7],se[6:7]),se[8]))
}
