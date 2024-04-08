







sixmultinomvineloglik.beta<-function(param,
                                     y001,y011,y101,y111,
                                     y000,y010,y100,y110,
              gl,mgrid,qcond1,qcond2,qcond3,qcond4,qcond5,
              tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,sel1,sel2,sel3)
{ p=param[1:6]
  g=param[7:12]
  tau=param[13:17]
  if(any(p<=0 | p>=1)) return(1.e10)
  if(1-sum(p[1:3])<= 0 | 1-sum(p[1:3])>=1) return(1.e10)
  if(1-sum(p[4:6])<= 0 | 1-sum(p[4:6])>=1) return(1.e10)
  if(any(g<=0 | g>=1)) return(1.e10)
  
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




multinom6dVineCopulaREMADA.beta=function(y001,y011,y101,y111,
                                          y000,y010,y100,y110,
                                          gl,mgrid,
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

  x1=p101
  x2=p011/(1-x1)
  x3=p111/(1-x2)/(1-x1)
  x4=p100
  x5=p010/(1-x4)
  x6=p110/(1-x5)/(1-x4)

  x=cbind(x1,x2,x3,x4,x5,x6)
  meanx <- apply(x, 2, mean)
  varx <- apply(x, 2, var)
  g = varx/meanx/(1 - meanx)
  stau = cor(x, method = "kendall")

  stau=c(stau[1,2],stau[2,3],stau[3,4],0.1,stau[5,6])
  inipar=c(p,g,stau)
  est=nlm(sixmultinomvineloglik.beta,inipar,y001,y011,y101,y111,
        y000,y010,y100,y110,
        gl,mgrid,qcond1,qcond2,qcond3,qcond4,qcond5,
        tau2par1,tau2par2,tau2par3,tau2par4,tau2par5,sel1,sel2,sel3,
        hessian=T,print.level = 2,iterlim=1000)
  est
}
