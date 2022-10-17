library(gamlj)
library(paml)
library(simr)
library(pbkrtest)
data("beers_bars")
data<-beers_bars
data$bar<-factor(data$bar)
zdata<-.levelstandard(data,"beer","beer","bar","case")

model1<-lmer(smile~zw_beer+(1+zw_beer|bar),data=zdata)

summary(model1)
vv<-VarCorr(model1)
aa<-vv$bar
aa<-matrix(c(3,-0.1,-0.1,8.113),ncol=2)
model2<-simr::makeLmer(smile~zw_beer+(1+zw_beer|bar), fixef = c(5,1), VarCorr=aa, sigma=sigma(model1),data=zdata)
summary(as_lmerModLmerTest(model2))
get_Lb_ddf(model2, c(0,1))
pdata<-zdata
pdata$smile<-predict(model2)
p<-ggplot(pdata,aes(y=smile,x=zw_beer,color=bar))
p<-p+geom_point()
p<-p+geom_smooth(method = "lm")
p



op<-simr::powerSim(model2,test=simr::fixed("zw_beer",method="sa"),nsim=500)
op
sp<-paml::modelpower(model2,method="sim",df="kr",nsim = 1000)
sp
paml::modelpower(model2,method="test",df="kr")
paml::modelpower(model2,method="test",df="sat")
nv<-2
lmat<-diag(1,nv,nv)
res<-as.data.frame(apply(lmat,1, function (l) {
  unlist(lmerTest::calcSatterth(model2,l))
}))
res
res$t<-sqrt(res$Fstat)

summary(as_lmerModLmerTest(model2))$coefficients

library(ggplot2)
p<-ggplot(zdata,aes(y=smile,x=zw_beer,color=bar))
p<-p+geom_point()
p<-p+geom_smooth(method = "lm")
p
(model0 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
library(simr)
model<-model0
fixef(model)<-c(250,5)
vv<-VarCorr(model)
vv$Subject
VarCorr(model)<-matrix(c(600,9,9,55),ncol=2)
op<-simr::powerSim(model,test=simr::fixed("Days",method="sa"),nsim=1000)
op
sp<-paml::modelpower(model,method="sim",df="sat",nsim = 1000)
sp
paml::modelpower(model,method="test",df="kr")
paml::modelpower(model,method="test",df="sat")


## removing Days
(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
anova(fmLarge,fmSmall)
library(pbkrtest)
KRmodcomp(fmLarge, fmSmall)  ## 17 denominator df's
get_Lb_ddf(model2, c(1,0))  ## 17 denominator df's
summary(model1)
