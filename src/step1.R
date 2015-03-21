#####################################################
#######################Step 1########################
#####################################################

sspp<-function(n=100, ch.lp='s1', ch.mc=c('brw','l','cps','ocps'), rec=c('int','bot','top')){
	para=c(0.8,6,7,3.5,.25,.3)
	#Prey specif. param.
	s0=para[1]
	alpha=para[2]
	beta=para[3]
	a=para[4]
	b=para[5]
	c=para[6]
  
	# if (ch.lp!='a' & ch.lp!='k' & ch.lp!=alpha & ch.lp!=beta){ 
	# if(ch.mc=='brw'){
	# str=paste(ch.lp,'=pnorm(brw(n, n0=',ch.lp,', sig=0.02, seed=',seed[1],'),', ch.lp,', 10',sep='')
	# } else if(ch.mc='cps') {
	# str=paste(ch.lp,'=',sep='')
	# }
	# } else {
	# if(ch.mc=='brw'){
	# str=paste(ch.lp,'=brw(n, n0=',ch.lp,', sig=0.02, seed=',seed[1],')',sep='')
	# }
	# }
	
	# eval(parse(text=str))
	
	# A stationary process in the [-1,1] interval
	ex=pnorm(brw(n,9,sig=0.3),9,5)
	ex.s=(ex-mean(ex))/max(abs(ex-mean(ex)))
	
	N=rep(0,n)
	P=N
	N[1]=7
	P[1]=0.5
	switch(rec[1],
	int={
	s0=s0+0.03*ex.s
	lf.m=s0
	for (i in 2:n){
	v1=s0[i]*(1-alpha*P[i-1]/(beta*P[i-1]+N[i-1]))*N[i-1]
	N[i]=a/(1+b*v1)*v1
    P[i]=c*alpha*N[i-1]/(beta*P[i-1]+N[i-1])*P[i-1]
	}},
	bot={
	s0=s0+0.1*ex.s
	lf.m=s0
	for (i in 2:n){
	v1=s0[i]*N[i-1]
    N[i]=a/(1+b*v1)*v1
    P[i]=c*alpha*N[i-1]/(beta*P[i-1]+N[i-1])*P[i-1]
	}},
	top={
	c=c+0.1*ex.s
	lf.m=c
	for (i in 2:n){
	v1=s0*(1-alpha*P[i-1]/(beta*P[i-1]+N[i-1]))*N[i-1]
    N[i]=a/(1+b*v1)*v1
    P[i]=c[i]*alpha*2/(beta*P[i-1]+2)*P[i-1]
	}}
	)
	# Here, we set the type of interactions!!
	# bot=1
	# if(rec[1]=='bot'){
		# bot=0
		# s1=s1+0.2*ex.s
		# lf.m=s1
		# for (i in 2:n) {
			# M1[]=c(s1[i]*(1-m1[i]), s1[i]*m1[i], a[i]/(1+(a[i]-1)/k[i]*N[2,i-1]), p1[i]*(1-bot*alpha[i]*P[2,i-1]/(beta[i]*P[2,i-1]+N[2,i-1])))
			# M2[]=c(s2[i]*(1-m2[i]), s2[i]*m2[i], c[i]*alpha[i]*N[2,i-1]/(beta[i]*P[2,i-1]+N[2,i-1]), p2[i])
			# W=matrix(c(rnorm(1,0,.01),rnorm(1,0,.01)),2,1)
			# V=matrix(c(rnorm(1,0,0.001),rnorm(1,0,0.001)),2,1)
			# N[,i]=M1%*%N[,i-1]+W
			# P[,i]=M2%*%P[,i-1]+V
	# }
	# } else if (rec[1]=='top'){
		# s2=s2+0.1*ex.s
		# lf.m=s2
		# for (i in 2:n) {
			# M1[]=c(s1[i]*(1-m1[i]), s1[i]*m1[i], a[i]/(1+(a[i]-1)/k[i]*N[2,i-1]), p1[i]*(1-bot*alpha[i]*P[2,i-1]/(beta[i]*P[2,i-1]+N[2,i-1])))
			# M2[]=c(s2[i]*(1-m2[i]), s2[i]*m2[i], c[i]*alpha[i]*0.5/(beta[i]*P[2,i-1]+0.5), p2[i])
			# W=matrix(c(rnorm(1,0,0.01),rnorm(1,0,0.01)),2,1)
			# V=matrix(c(rnorm(1,0,0.001),rnorm(1,0,0.001)),2,1)
			# N[,i]=M1%*%N[,i-1]+W
			# P[,i]=M2%*%P[,i-1]+V
		# }
	# } else {
		# s1=s1+0.03*ex.s
		# lf.m=s1
		# for (i in 2:n) {
			# M1[]=c(s1[i]*(1-m1[i]), s1[i]*m1[i], a[i]/(1+(a[i]-1)/k[i]*N[2,i-1]), p1[i]*(1-bot*alpha[i]*P[2,i-1]/(beta[i]*P[2,i-1]+N[2,i-1])))
			# M2[]=c(s2[i]*(1-m2[i]), s2[i]*m2[i], c[i]*alpha[i]*N[2,i-1]/(beta[i]*P[2,i-1]+N[2,i-1]), p2[i])
			# W=matrix(c(rnorm(1,0,0.01),rnorm(1,0,0.01)),2,1)
			# V=matrix(c(rnorm(1,0,0.001),rnorm(1,0,0.001)),2,1)
			# N[,i]=M1%*%N[,i-1]+W
			# P[,i]=M2%*%P[,i-1]+V
		# }
	# }
	dent=list(N=N,P=P,lf=lf.m)
	return(dent)
	
}


simu<-function(i=1, rec=c('int','top','bot','non','non.t','non.b'), range=30, drax=T, size=c(1,1)){
switch(rec[1],
	int={t=sspp(n=210+range,rec='int');
	dat=data.frame(N=(t$N),P=(t$P))},
	top={t=sspp(n=210+range,rec='top')
		dat=data.frame(N=(t$N),P=(t$P))},
	bot={t=sspp(n=210+range,rec='bot')
		dat=data.frame(N=(t$N),P=(t$P))},
	non={t=sspp(n=210+range,rec='int')
		 s=sspp(n=210+range,rec='int')
		 dat=data.frame(N=t$N,P=s$P)},
	non.t={t=sspp(n=210+range,rec='top')
		 s=sspp(n=210+range,rec='top')
		 dat=data.frame(N=t$N,P=s$P)},
	non.b={t=sspp(n=210+range,rec='bot')
		 s=sspp(n=210+range,rec='bot')
		 dat=data.frame(N=t$N,P=s$P)}
)
dat.t=dat[1:range+200,]  # training data
dat.tp=dat[200+range+1:10,] # prediction data

if(drax==T){
	obs.error=data.frame(rnorm(range+10, 0, size[1]), rnorm(range+10, 0, size[2]))
	dat.t=dat.t+obs.error[1:range,]
	dat.tp=dat.tp+obs.error[range+1:10,]
}

fi=VARselect(dat.t, lag.max=4, type='const')
k=fi$selection[3]
if(k==4 | k==1){k=2}
f=ca.jo(dat.t, type='trace', ecdet='const',K=k, spec='transitory')
# f.fit=vec2var(f,1)

# f.fitt=fitted(f.fit)
# f.er=resid(f.fit)
# wx=boot(f.er, jo.boot, R=100, fitt=f.fitt, k=k, init=dat.t[1:k,])
# teststat=rep(0,2)
# for(i in 1:2) {
#   boot.ci(cajo.boot, index=i, type='perc')->xv
#   teststat[i]=xv[[4]][4]
# }
# teststat=apply(wx$t, 2, mean)

if(f@teststat[2]>f@cval[2,2]) {z2=1} else {z2=0}
if(f@teststat[1]<f@cval[1,2]) {z1=1} else {z1=0}

# if(z2==0){
# r=0
# dt.dat.t=diff(ts(dat.t))
# fit=VAR(dt.dat.t,p=k-1, type='const')
# fit.v=VAR(dat.t,p=k, type='const')
# classi(fit,k=k-1,y=dt.dat.t)->stater
# classi(fit.v,k=k, y=dat.t)->statev
# signx=0
# } else 

if (z1==1 & z2==0){
signx=1
stater=classv(f,k=k, y=dat.t)
fit.v=VAR(dat.t,p=k, type='const')
classi(fit.v,k=k, y=dat.t)->statev
} else {
  signx=2
  fit=VAR(dat.t,p=k, type='const')
  classi(fit,k, y=dat.t)->statev
  stater=statev
}
# predict(fit,n.ahead=pd)->wrec
# predict(fit.v, n.ahead=pd)->wvar
# 
# if (signx==0){
# if(pd==1){
# level=c(dat.t[range,1],dat.t[range,2])
# } else {
# level=c(dat.t[range,1], dat.t[range,1]+wrec$fcst$N[1:(pd-1),1],
# 		dat.t[range,2], dat.t[range,2]+wrec$fcst$P[1:(pd-1),1])
# }} else {level=0}
# 
# goal=(c(wrec$fcst$N[1:pd,1], wrec$fcst$P[1:pd,1])+level-
# 	  c(dat.tp[1:pd,1], dat.tp[1:pd,2]))^2-
# 	  (c(wvar$fcst$N[1:pd,1], wvar$fcst$P[1:pd,1])-
# 	  c(dat.tp[1:pd,1], dat.tp[1:pd,2]))^2
# goal=sum(c(goal[1:pd]/mean(dat.t[,1]),goal[pd+1:pd]/mean(dat.t[,2])))
# if(goal<0){perf=1} else {perf=0}

## Correlation coefficient M
 statec=classc(dat.t[,1], dat.t[,2])

## Regression method
statel=classl(dat.t[,1], dat.t[,2])

rzl=cbind(i, z1,z2,stater, statev, 0 , signx, statec, statel)
return(rzl)
}




















