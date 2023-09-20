
library(paml)


sigma2<-matrix(c(1,.2,.2,1), ncol=2)

clusters<-list(
          list(name="level3",n=10,type="nested"),
          list(name="level2",n=10,type="nested")
          )
vars<-list(list(name="y",n=30,varying="dependent"),
           list(name="x",type="numeric",varying="within"),
           list(name="f",type="nominal", n=2,varying="level2")
)
sigmaw<-matrix(c(1,.5,.5,1), ncol=2)
sample<-.asample(vars,clusters,sigmaw=sigmaw,type = "nested")
names(sample)
sample$f<-factor(sample$f)

mod<-lmer(y~x*f+(1+x|level2)+(1+f*x|level3),data=sample)
round(fixef(mod),digits=3)
mod
library(simr)
VarCorr(mod)
xsample<-sample

varcorr<-list(matrix(c(1,0,0,.2),ncol = 2),
              matrix(c(1,0,0,0,
                       0,.2,0,0,
                       0,0,.1,0,
                       0,0,0,.1),ncol = 4))
anew<-simr::makeLmer(y~x*f+(1+x|level2)+(1+f*x|level3), c(1,.5,.2,.5), VarCorr=varcorr, sigma=2,data=xsample)
anew
model<-anew
summary(model)
ses<-summary(as_lmerModLmerTest(model))
ses
summary(model)
library(paml)

#sp<-paml::modelpower(anew,method="sim",nsim = 1000)
sp
summary(model)
spt<-paml::modelpower(model)
spt

.summary<-summary(as_lmerModLmerTest(model))
.coefs<-.summary$coefficients[,1]
.se<-.summary$coefficients[,2]
df<-.summary$coefficients[,3]
ncp<-.coefs/.se
tside<-2
crit<-qt(alpha/tside,df,lower.tail = F)
power<-pt(crit, df, ncp = ncp, lower.tail = FALSE) +  pt(-crit, df, ncp = ncp, lower.tail = TRUE)



b<-fixef(anew)[3]
se<-ses$coefficients[3,2]
ncp<-b/se
df<-ses$coefficients[3,3]

sig.level<-.05
tside<-2
t0 <- qt(sig.level/tside,df=df)
(1-pt(t0,df=df,ncp=ncp)-pt(-t0,df=df,ncp=ncp))

crit<-qt(sig.level/tside,df,lower.tail = F)
pt(crit, df, ncp = ncp, lower.tail = FALSE) +  pt(-crit, df, ncp = ncp, lower.tail = TRUE)


b<-fixef(anew)[4]
se<-ses$coefficients[4,2]
ncp<-b/se
df<-ses$coefficients[4,3]

crit<-qt(sig.level/tside,df,lower.tail = F)
pt(crit, df, ncp = ncp, lower.tail = FALSE) +  pt(-crit, df, ncp = ncp, lower.tail = TRUE)


table(sample$level1,sample$f)
table(sample$level2,sample$f)
table(sample$level3,sample$f)
amat<-matrix(c(1,.5,.5,1),ncol=2)
wmat<-matrix(c(1,.1,.1,1),ncol=2)
wvar<-sqrt(2*(1-abs(amat))/12)
diag(wvar)<-1
wvar
sample<-.asample(vars,clusters,type = "nested",sigmaw = amat)
names(sample)
sample$xx<-sample$x*.787
sample$y<-sample$y+1
smodel<-lmer(y~xx+(1|bar),data=sample)
ss<-summary(smodel)

sdy<-mean(tapply(sample$y, sample$bar, sd))
sdx<-mean(tapply(sample$xx, sample$bar, sd))
ub<-ss$coefficients[2,1]
ub*sdx/sdy

zsample<-.levelstandard(sample,"y","xx","bar")
zmodel<-lmer(y~zxx+(1|bar),data=zsample)
summary(zmodel)


anew<-simr::makeLmer(y~x+(1+x|bar), fixef(amod), VarCorr=varcorr, sigma=sigma(amod),data=data)
summary(as_lmerModLmerTest(anew))
ndata<-data
ndata$y<-unlist(simulate(anew))
amod<-lmer(y~x+(1+x|bar),data=ndata)
amod
lapply(unique(ndata$bar),function(x) cor(ndata[ndata$bar==x,c("y","x")]))
var(ndata$y)
tapply(ndata$y,ndata$bar,var)

data$y<-data$y+.1
amod<-lmer(y~x+(1|bar),data=data)
summary(amod)

data$f<-factor(data$f)
fixed=c(0,.5)
VarCorr<-matrix(c(1,0.1,0.1,1),2)
sigma<-1
mm<-lmer(y~f+(1+x|bar),data=data)
mm
amod<-simr::makeLmer(y~x+(1+x|bar), fixed, VarCorr=VarCorr, sigma,data=data)
amod
amod2<-simr::extend(amod,within="bar",n=50)
mp<-modelpower(amod2,effects = "random")
mp
required_n(amod2,"bar",direction = "between",nsim=500)
