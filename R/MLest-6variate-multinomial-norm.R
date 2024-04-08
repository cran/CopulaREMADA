sixmultinomprod=function (x1,x2,x3,x4,x5,x6, 
                       y001,y011,y101,y111,
                       y000,y010,y100,y110) 
{ f1=y101*log(x1)+y011*log(x2)+y111*log(x3)+y001*log(1-x1-x2-x3)
  f2=y100*log(x4)+y010*log(x5)+y110*log(x6)+y000*log(1-x4-x5-x6)
  exp(f1+f2)
}




sixmultinomvineloglik.norm<-function(param,
y001,y011,y101,y111,
y000,y010,y100,y110,gl,mgrid,
qcond1,qcond2,qcond3,qcond4,qcond5,
tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,sel1,sel2,sel3)
{ p=param[1:6]
  si=param[7:12]
  tau=param[13:17]
  if(any(p<=0 | p>=1)) return(1.e10)
  if(1-sum(p[1:3])<= 0 | 1-sum(p[1:3])>=1) return(1.e10)
  if(1-sum(p[4:6])<= 0 | 1-sum(p[4:6])>=1) return(1.e10)
  if(any(si<=0)) return(1.e10)
  
  if(any(tau[sel1]>=0.95 | tau[sel1]<=0)) return(1.e10)
  if(any(tau[sel2]<= -0.95 | tau[sel2]>=0)) return(1.e10)
  if(any(tau[sel3]<= -0.95 | tau[sel3]>=0.95)) return(1.e10)
  
  th=c(tau2par1(tau[1]),tau2par2(tau[2]),tau2par3(tau[3]),
       tau2par4(tau[4]),tau2par5(tau[5]))
  
  
  p1=mgrid$x1
  p2=mgrid$x2
  p3=mgrid$x3
  p4=mgrid$x4
  p5=mgrid$x5
  p6=mgrid$x6
  
  u1=p1
  u2=qcond1(p2,u1,th[1])
  u3=qcond2(p3,u2,th[2])
  u4=qcond3(p4,u3,th[3])
  u5=qcond4(p5,u4,th[4])
  u6=qcond5(p6,u5,th[5])
  
  mu1=log(p[1]/(1-sum(p[1:3])))
  mu2=log(p[2]/(1-sum(p[1:3])))
  mu3=log(p[3]/(1-sum(p[1:3])))
  mu4=log(p[4]/(1-sum(p[4:6])))
  mu5=log(p[5]/(1-sum(p[4:6])))
  mu6=log(p[6]/(1-sum(p[4:6])))
 
  
  x1=exp(qnorm(u1,mu1,si[1]))
  x2=exp(qnorm(u2,mu2,si[2]))
  x3=exp(qnorm(u3,mu3,si[3]))
  x4=exp(qnorm(u4,mu4,si[4]))
  x5=exp(qnorm(u5,mu5,si[5]))
  x6=exp(qnorm(u6,mu6,si[6]))
  
  t1=1+x1+x2+x3
  t2=1+x4+x5+x6
  
  p1=x1/t1
  p2=x2/t1
  p3=x3/t1
  p4=x4/t2
  p5=x5/t2
  p6=x6/t2
  
  N=length(y001)
  prob=rep(NA,N)
  for(i in 1:N)
  #out=foreach(i=1:N,.combine='+') %dopar%
  { temp=sixmultinomprod(p1,p2,p3,p4,p5,p6, 
                         y001[i],y011[i],y101[i],y111[i],
                         y000[i],y010[i],y100[i],y110[i])
    
    prob[i]=tensor(tensor(tensor(tensor(tensor(temp,gl$w,
    6,1),gl$w,5,1),gl$w,4,1),gl$w,3,1),gl$w,2,1)%*%gl$w
  
    #if(prob!=0) -log(prob) else 0
  }
  #out
  -sum(log(prob[prob!=0]))
}




multinom6dVineCopulaREMADA.norm=function(y001,y011,y101,y111,
                                y000,y010,y100,y110,gl,mgrid,
                                qcond1,qcond2,qcond3,qcond4,qcond5,
                    tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,sel1,sel2,sel3)
{ y001.c=y001 + 0.5*(y001==0)
  y011.c=y011 + 0.5*(y011==0)
  y101.c=y101 + 0.5*(y101==0)
  y111.c=y111 + 0.5*(y111==0)
  y000.c=y000 + 0.5*(y000==0)
  y010.c=y010 + 0.5*(y010==0)
  y100.c=y100 + 0.5*(y100==0)
  y110.c=y110 + 0.5*(y110==0)
  
  n1=y001.c+y011.c+y101.c+y111.c
  n2=y000.c+y010.c+y100.c+y110.c
  
  p011=y011.c/n1
  p101=y101.c/n1
  p111=y111.c/n1
  p010=y010.c/n2
  p100=y100.c/n2
  p110=y110.c/n2
  
  p=apply(cbind(p101,p011,p111,p100,p010,p110),2,mean)
  
  z1=log(p101/(1-p011-p101-p111))
  z2=log(p011/(1-p011-p101-p111))
  z3=log(p111/(1-p011-p101-p111))
  z4=log(p100/(1-p010-p100-p110))
  z5=log(p010/(1-p010-p100-p110))
  z6=log(p110/(1-p010-p100-p110))
  z=cbind(z1,z2,z3,z4,z5,z6)
  si<-sqrt(apply(z,2,var))
  stau=cor(z,method="kendall")
  stau=c(stau[1,2],stau[2,3],stau[3,4],0.1,stau[5,6])
  
  inipar=c(p,si,stau)
  est=nlm(sixmultinomvineloglik.norm,inipar,y001,y011,y101,y111,
          y000,y010,y100,y110,gl,mgrid,
          qcond1,qcond2,qcond3,qcond4,qcond5,
          tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,sel1,sel2,sel3,
          hessian=T,print.level =2,iterlim=1000)
  est
}


