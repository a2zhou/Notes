proc.time()[3]->t
library(vars)
library(urca)
library(snow)
library(snowfall)
library(rlecuyer)
library(boot)
source('step1.R')
source('brw.R')
source('step2.class.R')


n=1000
pd=2
x=30
knife=T
obs.e.r=10

sfInit(parallel=T, cpus=3, type='SOCK')
sfLibrary(vars)
sfLibrary(urca)
sfLibrary(boot)
sfExportAll()
sfClusterSetupRNG()

# obs.error
t.int=sspp(n=2000, rec='int')
var.P2=var(t.int$P[200:2000])
var.N2=var(t.int$N[200:2000])
# 5% obs.error rate
# 10% obs.error rate
# 20% obs.error rate
size.int=sqrt(obs.e.r/100*c(var.N2, var.P2))

rzl<-sfLapply(1:n, simu, rec='int', range=x, drax=knife, size=size.int)
z1.int=sapply(rzl, function(x){x[1,2]})
z2.int=sapply(rzl, function(x){x[1,3]})
stater.int=sapply(rzl, function(x){x[1,4]})
statev.int=sapply(rzl, function(x){x[1,5]})
perf.int=sapply(rzl, function(x){x[1,6]})
signx.int=sapply(rzl, function(x){x[1,7]})
statec.int=sapply(rzl, function(x){x[1,8]})
statel.int=sapply(rzl, function(x){x[1,9]})
rm(list='rzl')


# obs.error
t.bot=sspp(n=2000, rec='bot')
var.P2=var(t.bot$P[200:2000])
var.N2=var(t.bot$N[200:2000])
# 5% obs.error rate
# 10% obs.error rate
# 20% obs.error rate
size.bot=sqrt(obs.e.r/100*c(var.N2, var.P2))

rzl<-sfLapply(1:n, simu, rec='bot', range=x, drax=knife, size=size.bot)
z1.bot=sapply(rzl, function(x){x[1,2]})
z2.bot=sapply(rzl, function(x){x[1,3]})
stater.bot=sapply(rzl, function(x){x[1,4]})
statev.bot=sapply(rzl, function(x){x[1,5]})
perf.bot=sapply(rzl, function(x){x[1,6]})
signx.bot=sapply(rzl, function(x){x[1,7]})
statec.bot=sapply(rzl, function(x){x[1,8]})
statel.bot=sapply(rzl, function(x){x[1,9]})
rm(list='rzl')


# obs.error
t.top=sspp(n=2000, rec='top')
var.P2=var(t.top$P[200:2000])
var.N2=var(t.top$N[200:2000])
# 5% obs.error rate
# 10% obs.error rate
# 20% obs.error rate
size.top=sqrt(obs.e.r/100*c(var.N2, var.P2))

rzl<-sfLapply(1:n, simu, rec='top', range=x, drax=knife, size=size.top)
z1.top=sapply(rzl, function(x){x[1,2]})
z2.top=sapply(rzl, function(x){x[1,3]})
stater.top=sapply(rzl, function(x){x[1,4]})
statev.top=sapply(rzl, function(x){x[1,5]})
perf.top=sapply(rzl, function(x){x[1,6]})
signx.top=sapply(rzl, function(x){x[1,7]})
statec.top=sapply(rzl, function(x){x[1,8]})
statel.top=sapply(rzl, function(x){x[1,9]})
rm(list='rzl')

# obs.error

rzl<-sfLapply(1:n, simu, rec='non', range=x, drax=knife, size=size.int)
z1.non=sapply(rzl, function(x){x[1,2]})
z2.non=sapply(rzl, function(x){x[1,3]})
stater.non=sapply(rzl, function(x){x[1,4]})
statev.non=sapply(rzl, function(x){x[1,5]})
perf.non=sapply(rzl, function(x){x[1,6]})
signx.non=sapply(rzl, function(x){x[1,7]})
statec.non=sapply(rzl, function(x){x[1,8]})
statel.non=sapply(rzl, function(x){x[1,9]})
rm(list='rzl')


# obs.error

rzl<-sfLapply(1:n, simu, rec='non.t', range=x, drax=knife, size=size.top)
z1.non.t=sapply(rzl, function(x){x[1,2]})
z2.non.t=sapply(rzl, function(x){x[1,3]})
stater.non.t=sapply(rzl, function(x){x[1,4]})
statev.non.t=sapply(rzl, function(x){x[1,5]})
perf.non.t=sapply(rzl, function(x){x[1,6]})
signx.non.t=sapply(rzl, function(x){x[1,7]})
statec.non.t=sapply(rzl, function(x){x[1,8]})
statel.non.t=sapply(rzl, function(x){x[1,9]})
rm(list='rzl')


# obs.error

rzl<-sfLapply(1:n, simu, rec='non.b', range=x, drax=knife, size=size.bot)
z1.non.b=sapply(rzl, function(x){x[1,2]})
z2.non.b=sapply(rzl, function(x){x[1,3]})
stater.non.b=sapply(rzl, function(x){x[1,4]})
statev.non.b=sapply(rzl, function(x){x[1,5]})
perf.non.b=sapply(rzl, function(x){x[1,6]})
signx.non.b=sapply(rzl, function(x){x[1,7]})
statec.non.b=sapply(rzl, function(x){x[1,8]})
statel.non.b=sapply(rzl, function(x){x[1,9]})
rm(list='rzl')

sfStop()
te=proc.time()[3]-t
print(te)

length(stater.int[which(stater.int=='int')])/n
length(statev.int[which(statev.int=='int')])/n

length(stater.non[which(stater.non=='non')])/n
length(statev.non[which(statev.non=='non')])/n

length(stater.bot[which(stater.bot=='bot')])/n
length(statev.bot[which(statev.bot=='bot')])/n

length(stater.non.b[which(stater.non.b=='non')])/n
length(statev.non.b[which(statev.non.b=='non')])/n

length(stater.top[which(stater.top=='top')])/n
length(statev.top[which(statev.top=='top')])/n

length(stater.non.t[which(stater.non.t=='non')])/n
length(statev.non.t[which(statev.non.t=='non')])/n


