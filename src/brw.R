brw<-function(n, n0=7,sig=0.02){
alpha1=1
alpha2=1
k=-6
tau=n0

x=rep(0,n)
x[1]=n0
for(i in 2:n) { 
    x[i]=x[i-1]+exp(k)*(exp(-alpha1*(x[i-1]-tau))-exp(alpha2*(x[i-1]-tau)))+rnorm(1,0,sig)
}
return(x)
}

bottom<-function(n=200, n0=c(7,2.6), pr=1) {
  s0=0.9
  alpha=6
  beta=7
  a=3.5
  b=0.25
  c=0.6
  p=rep(0,n)
  v=brw(n, n0[1])
  p[1]=n0[2]
  
  for (i in 2:n){
    p[i]=c*alpha*v[i-1]/(beta*p[i-1]+v[i-1])*p[i-1]+pr*rnorm(1,sd=0.02)
  }
  den=matrix(c(v,p),n,2)
  return(den)
}


top<-function(n=200, n0=c(8,0.4)){
  s0=0.9
  alpha=6
  beta=7
  a=3.5
  b=0.25
  c=0.6
  v=rep(0,n)
  p=brw(n, n0[2])
  v[1]=n0[1]
  for (i in 2:n){
    v1=s0*(1-alpha*p[i-1]/(beta*p[i-1]+v[i-1]))*v[i-1]
    v[i]=a/(1+b*v1)*v1+rnorm(1,sd=0.1)
  }
  den=matrix(c(v,p),n,2)
  return(den)
}


iso<-function(range=c(7,1)){
  s0=0.9
  alpha=6
  beta=7
  a=3.5
  b=0.25
  c=0.6
  
  min=range[1]
  span=range[2]
  x=min+span*(0:100)/100
  y=x*(-b*s0*x+a*s0-1)/(b*s0*(beta-alpha)*x+beta-a*s0*(beta-alpha))
  return(matrix(c(x,y),ncol=2))
}

pp<-function(para=c(0.9,6,7,3.5,.25,.3),n=100, n0=c(7,0.5), pry=F, prd=F, sed=35, v=1) {
  if (length(para)!=6) {return(NULL)}
  s0=para[1]
  alpha=para[2]
  beta=para[3]
  a=para[4]
  b=para[5]
  c=para[6]
  
  if (pry==T) {brw(n, n0=s0, sig=0.01*v)->s.a} else {s.a=rep(s0,n)}
  if (prd==T) {brw(n, c, sig=v*0.001)->s.b} else {s.b=rep(c,n)}
  
  p=rep(0,n)
  v=rep(0,n)
  p[1]=n0[2]
  v[1]=n0[1]
  for (i in 2:n){
    v1=s.a[i]*(1-alpha*p[i-1]/(beta*p[i-1]+v[i-1]))*v[i-1]
    v[i]=a/(1+b*v1)*v1
    p[i]=s.b[i]*alpha*v[i-1]/(beta*p[i-1]+v[i-1])*p[i-1]
  }
  den=matrix(c(v,p),n,2)
  rzlt=list(data=den,para=para, s.a=s.a, s.b=s.b)
  return(rzlt)
}

#iso for pp
iso<-function(para=c(0.9,6,7,3.5,.25,.6), range=c(7,1)){
  s0=para[1]
  alpha=para[2]
  beta=para[3]
  a=para[4]
  b=para[5]
  c=para[6]
  
  min=range[1]
  span=range[2]
  x=min+span*(0:100)/100
  y=x*(-b*s0*x+a*s0-1)/(b*s0*(beta-alpha)*x+beta-a*s0*(beta-alpha))
  ypd=(c*alpha-1)/beta*x
  
  dat=matrix(c(x,y,ypd),ncol=3)
  rzlt=list(data=dat,para=para)
  return(rzlt)
}
