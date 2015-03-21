classi.t<-function(w,k=2,p=0.05){
if(!(class(w)=="varest")){
stop("\nPlease provide object of class 'varest' as 'z'\n")
}
summary(w)->w
flag=c(0,0)
col.name=names(w$varresult)
# top down evidence
coeff=w$varresult[[col.name[1]]]$coefficients[,4]
sequ=c((1:k)*2)
for(i in sequ){
	if(coeff[i]<p){flag[1]=1}
}
# bottom up evidence
coeff=w$varresult[[col.name[2]]]$coefficients[,4]
for(i in sequ-1){
	if(coeff[i]<p){flag[2]=1}
}

# Classification
clas=NULL
if (flag[1]==1){
	if(flag[2]==1){
	clas='int'
	} else {clas='top'}
} else if (flag[2]==1){
	clas='bot'
} else {clas='non'}
return(clas)
}

classv<-function(w, k=2, p=0.05, y=dat.t){
if(!(class(w)=="ca.jo")){
stop("\nPlease provide object of class 'ca.jo' as 'z'\n")
}
flag=c(0,0)
w.fit=vec2var(w,1)

w.fitt=fitted(w.fit)
w.er=resid(w.fit)
seq.1=(1:k)*2
seq.2=(1:k)*2-1
init=y[1:k,]
wx=boot(w.er, vec.boot, R=100, fitt=w.fitt, k=k, seq.1=seq.1, seq.2=seq.2, init=init)

for(i in 1:k) {
boot.ci(wx, index=i, type='perc')->xv
if (xv[[4]][4]*xv[[4]][5]>0){flag[1]=1;break}
}
for(i in (k+1):(2*k)){
boot.ci(wx, index=i, type='perc')->xv
if (xv[[4]][4]*xv[[4]][5]>0){flag[2]=1;break}
}
# Class
clas=NULL
if (flag[1]==1){
	if(flag[2]==1){
	clas='int'
	} else {clas='top'}
} else if (flag[2]==1){
	clas='bot'
} else {clas='non'}
return(clas)
}

classc<-function(x, y, p=0.05){
fit=cor.test(x, y)
p.v=fit$p.value
estm=fit$estimate

# Class
clas=NULL
if (p.v<p){
	if (estm>0) {
	clas='bot'
	} else { clas='top'}
} else {clas='ind'}
return(clas)
}

classl<-function(x, y, p=0.05){
l=length(x)
fit.xy=summary(lm(x[2:l]~y[1:(l-1)]))
fit.yx=summary(lm(y[2:l]~x[1:(l-1)]))
flag=c(0,0)
if (fit.xy$coefficients[2,4]<p ){flag[1]=1}
if (fit.yx$coefficients[2,4]<p ){flag[2]=1}
# Class
clas=NULL
if (flag[1]==1){
	if(flag[2]==1){
	clas='int'
	} else {clas='top'}
} else if (flag[2]==1){
	clas='bot'
} else {clas='non'}
return(clas)
}

classi.GC<-function(w,k=2,p=0.05){
if(!(class(w)=="varest")){
stop("\nPlease provide object of class 'varest' as 'z'\n")
}
flag=c(0,0)
col.names=names(w$varresult)
causality(w, cause=col.names[1], boot=T, boot.runs=500)->wx
if(wx$Granger$p.value<p){flag[1]=1}
causality(w, cause=col.names[2], boot=T, boot.runs=500)->xw
if(xw$Granger$p.value<p){flag[2]=1}
# Classification
clas=NULL
if (flag[1]==1){
	if(flag[2]==1){
	clas='int'
	} else {clas='top'}
} else if (flag[2]==1){
	clas='bot'
} else {clas='non'}
return(clas)
}

classi<-function(w,k=2,p=0.05,y){
if(!(class(w)=="varest")){
stop("\nPlease provide object of class 'varest' as 'z'\n")
}
flag=c(0,0)
w.fitt=fitted(w)
w.er=residuals(w)
seq.up=2*(1:k)
seq.dn=2*(0:(k-1))+1
kora.boot=boot(w.er, var.boot, R=100, k=k, n1=seq.up, n2=seq.dn, fitt=w.fitt, init=y[1:k,])
for(i in 1:k) {
boot.ci(kora.boot, index=i, type='perc')->xv
if (xv[[4]][4]*xv[[4]][5]>0){flag[1]=1;break}
}
for(i in (k+1):(2*k)){
boot.ci(kora.boot, index=i, type='perc')->xv
if (xv[[4]][4]*xv[[4]][5]>0){flag[2]=1;break}
}
# Classification
clas=NULL
if (flag[1]==1){
	if(flag[2]==1){
	clas='int'
	} else {clas='top'}
} else if (flag[2]==1){
	clas='bot'
} else {clas='non'}
return(clas)
}

var.boot<-function(data, indices, k, n1, n2, fitt, init){
	data.b=ts(rbind(init, fitt+data[indices,]))
	mod=VAR(data.b, p=k, type='const')
	coef(mod)->wx
	rzl=c(wx[[1]][n1,1],wx[[2]][n2,1])
	return(rzl)
}


vec.boot<-function(data, indices, k, fitt, seq.1, seq.2, init){
  data.b=rbind(ts(init), fitt+data[indices,])
  mod=ca.jo(data.b, type='trace', ecdet='const',K=k, spec='transitory')
  fit=vec2var(mod,1)
  rzl=fit[[2]][[1]]
  for(i in 1:k){
  rzl=cbind(rzl, fit[[2]][[i]])
  }
  return(c(rzl[1,seq.1], rzl[2, seq.2]))
}

jo.boot<-function(data, indices, k, fitt, init){
  data.b=rbind(ts(init), fitt+data[indices,])
  mod=ca.jo(data.b, type='trace', ecdet='const',K=k, spec='transitory')
  return(mod@teststat)
}





