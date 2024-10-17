expit=function(x){1/(1+exp(-x))}
logit=function(p){log(p/(1-p))}

se.p=function(x,se) {se*der.expit(x)}

der.expit=function(x)
{ temp=expit(x)
temp*(1-temp)
}
se.si=function(x,se) {se*exp(x)}


imperfect.multinomprod=function (lst,y11,y10,y01,y00) 
{ p11=lst$x1*lst$x2*lst$x3+(1-lst$x1)*(1-lst$x4)*(1-lst$x5)
  p10=lst$x1*lst$x2*(1-lst$x3)+(1-lst$x1)*(1-lst$x4)*lst$x5
  p01=lst$x1*(1-lst$x2)*lst$x3+(1-lst$x1)*lst$x4*(1-lst$x5)
  p00=1-p11-p10-p01
  f=y11*log(p11)+y10*log(p10)+y01*log(p01)+y00*log(p00)
  exp(f)
}


imperfect.REMADA.loglik.norm<-function(param,y11,y10,y01,y00,gl,select.random)
{ p=param[1:5]
  si=param[6]
  p=expit(p)
  if(p[3]+p[5]-1<0) return(1e10)
  si=exp(si)
  u=gl$n
  mu=log(p[select.random]/(1-p[select.random])) 
  x=qnorm(u,mu,si)
  t=exp(x)
  x=t/(1+t)
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  lst[[select.random]]=x
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N)   
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= temp%*%gl$w
  }
  -sum(log(prob[prob!=0]))
}






imperfectREMADA.norm=function(y11,y10,y01,y00,
                               gl,select.random,start)
{ start=c(logit(start[1:5]),log(start[6]))
  est=nlm(imperfect.REMADA.loglik.norm,start,y11,y10,y01,y00,gl,select.random,
          hessian=T,print.level=2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),exp(est$e[6])),
     SE=c(se.p(est$e[1:5],se[1:5]),se.si(est$e[6],se[6])))
}


imperfect.REMADA.loglik.beta<-function(param,y11,y10,y01,y00,gl,select.random)
{ p=param[1:5]
  g=param[6]
  p=expit(p)
  g=expit(g)
  u=gl$n
  a=p[select.random]/g-p[select.random]
  b=(1-p[select.random])*(1-g)/g
  x=qbeta(u,a,b)
  lst=list(x1=1,x2=1,x3=1,x4=1,x5=1)
  lst[[select.random]]=x
  for(ii in (1:5)[-select.random]) lst[[ii]]=p[ii]
  N=length(y11)
  prob=rep(NA,N)
  for(i in 1:N) 
  { temp=imperfect.multinomprod(lst,y11[i],y10[i],y01[i],y00[i])
    prob[i]= temp%*%gl$w
  }
  -sum(log(prob[prob!=0]))
}



imperfectREMADA.beta=function(y11,y10,y01,y00,
                                          gl,select.random,start)
{ start=logit(start)
  est=nlm(imperfect.REMADA.loglik.beta,start,y11,y10,y01,y00,
        gl,select.random,hessian=T,print.level= 2,iterlim=1000)
  se=sqrt(diag(solve(est$hessian)))
  list(LogLikelihood=-est$m,Estimates=c(expit(est$e[1:5]),expit(est$e[6])),
     SE=c(se.p(est$e[1:5],se[1:5]),se.p(est$e[6],se[6])))
}
