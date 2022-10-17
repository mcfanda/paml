library(lmerTest)
library(simr)
NC<-100
NU<-50

fixed<-c(.5,.1)
varcov<-c(3,14)
int<-rnorm(NC,fixed[1],sqrt(varcov[1]))
b<-rnorm(NC,fixed[2],sqrt(varcov[2]))
res<-list()
for (i in seq_len(NC)) {
  x<-rnorm(NU)
  y<-int[i]+b[i]*x+rnorm(NU)
  res[[length(res)+1]]<-cbind(cluster=i,x=x,y=y)
}
data<-as.data.frame(do.call(rbind,res))
data$cluster<-factor(data$cluster)
model<-lme4::lmer(y~x+(1+x|cluster),data=data)
summary(model)
summary(as_lmerModLmerTest(model))
paml::modelpower(model,method="test",df="sat")

mmodel<-model
vv<-VarCorr(model)
vv$cluster
VarCorr(mmodel)<-matrix(c(3,0,0,4),ncol=2)

simr::powerSim(mmodel,test=simr::fixed("x",method="sa"),nsim=500)

paml::modelpower(mmodel,method="sim",df="sat")
paml::modelpower(model,method="sim",df="sat")

paml::modelpower(model,method="test",df="sat")
paml::modelpower(mmodel,method="test",df="sat")
summary(mmodel)


a<-as_lmerModLmerTest(model)
summary(a)
b<-as_lmerModLmerTest(mmodel)
summary(b)
summary(as_lmerModLmerTest(model))
summary(as_lmerModLmerTest(mmodel))
#pbkrtest::get_Lb_ddf(model,c(0,1))
#pbkrtest::get_Lb_ddf(mmodel,c(0,1))

paml::modelpower(model,method="test",df="sat")


