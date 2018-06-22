# Functions for conditional cdf for 1-parameter bivariate copula families and t


# conditional cdfs of bivariate copula C_{2|1}(v,u,cpar) 
# inverse of conditional cdfs C_{2|1}^{-1}(p,u,cpar)
# 0<v<1, 0<u<1, 0<p<1, cpar=copula parameter
# Most functions here should work if u,v,p,cpar are vectors of the same length,
#  or if only one of the three is a vector and the other two are scalars.
# The boundary constraints on the functions are not checked on.

# C_{2|1}(v|u} for bivariate  normal copula
# cpar = copula parameter with -1<cpar<1
pcondbvn=function(v,u,cpar)
{ val=pnorm((qnorm(v)-cpar*qnorm(u))/sqrt(1-cpar^2))
  val[v <= 0 | u <= 0 | u >= 1]=0
  val[v == 1]=1
  val
}







# Frank
# cpar = copula parameter: cpar>0 or cpar<0 (latter for negative dependence)
pcondfrk=function(v,u,cpar)
{ #if(cpar==0.) return(v)
  cpar[cpar==0]=1.e-10
  cpar1=1.-exp(-cpar);
  tem=1.-exp(-cpar*u);
  ccdf=(1.-tem)/(cpar1/(1.-exp(-cpar*v))-tem);
  ccdf
}

# MTCJ
# cpar = copula parameter >0
pcondcln=function(v,u,cpar)
{ tem=v^(-cpar)-1
  tem=tem*(u^cpar)+1
  ccdf=tem^(-1-1/cpar)
  ccdf
}

pcondcln90=function(v,u,cpar)
{ cpar=-cpar
  u=1-u
  pcondcln(v,u,cpar)
}


pcondcln180=function(v,u,cpar)
{ u=1-u
  v=1-v
  1-pcondcln(v,u,cpar)
}


pcondcln270=function(v,u,cpar)
{ cpar=-cpar
  v=1-v
  1-pcondcln(v,u,cpar)
}

