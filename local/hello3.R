library(gamlj)
library(paml)
data("beers_bars")
data<-beers_bars
data$bar<-factor(data$bar)

zzdata<-.levelstandard(data,"beer","beer","bar","case")
cdata<-.levelcenter(data,"beer","bar",level="within")
cdata<-.levelcenter(cdata,"beer","bar",level="between")
head(cdata)
model1<-lmer(smile~beer+(1+beer|bar),data=data)
model2<-lmer(smile~centered_beer+mean_beer+(1+centered_beer|bar),data=cdata)
round(betas(model1),digits=3)

library(gamlj)
library(paml)
data<-subjects_by_stimuli
data$subj<-factor(data$subj)
data$stimulus<-factor(data$stimulus)
names(data)
model<-lmer(y~cond+(1|subj)+(1|stimulus),data=data)
tapply(data$cond, data$subj, unique)
.levelstandard(data,"y","cond","subj",level = "between")
betas(model)
data2<-data
data2$cond<-factor(data2$cond)
data2$cond2<-factor(rbinom(dim(data2)[1],3,.5))
tt<-tapply(data$y, data$subj, mean)
sd(tt)
tt<-tapply(data$y, list(data$subj), mean)
sd(tt)

model2<-lmer(y~cond+(1|subj)+(1|stimulus),data=data2)
betas(model2)
model3<-lmer(y~1+(1|stimulus)+(1|subj),data=data2)
model3
betas(model3)
model3<-lmer(y~cond+(1|subj/stimulus),data=data2)
summary(model3)
betas(model3)
ff<-y~cond+(1|subj/stimulus)
af<-y~subj:stimulus

x<-findbars(ff)[[1]][[3]]
as.character(x)
head(data)
names(data2)
.levelstandard(data2,"y","cond","subj",)
a<-c(0,5.2851)
sd(a)
mean(tapply(data2$y, data2$subj, sd))
model<-model2
.data<-model@frame
.names<-names(.data)
.terms<-colnames(model.matrix(model))[-1]
.dep<-formula(model)[[2]]
 cluster<-"bar"
.zdata<-.data
.zdata$".id."<-1:dim(.zdata)[1]

for (term in .terms) {

    if (term %in% .names) {
      c0<-tapply(.data[[term]], .data[[cluster]],unique)
      if (any(sapply(c0, length)>1)) level="within"  else level="between"
      cat("\n",term,level,"\n")
      .zdata<-.levelstandard(.zdata,.dep,term,cluster,level = level,id=".id.",overwrite=T)
    }

}
tapply(.zdata$smile, .zdata$bar, sd)
tapply(.zdata$centered_beer, .zdata$bar, sd)

mod<-update(model,data=.zdata)
ss<-summary(model)
cc<-as.data.frame(ss$coefficients)
cc$beta<-fixef(mod)
cc$beta[1]<-0
cc

xmodel<-lmer(smile~centered_beer+mean_beer+(1+centered_beer|bar),data=cdata)
summary(xmodel)
(ss<-summary(model))
sdy<-mean(tapply(data$smile, data$bar, sd))
sdx<-mean(tapply(data$beer, data$bar, sd))
ub<-ss$coefficients[2,1]
ub*sdx/sdy

data$id<-1:length(data$case)
zdata<-.levelstandard(cdata,"smile","centered_beer","bar",id="case",level="within")
zdata<-.levelstandard(zdata,"smile","mean_beer","bar",id="case",level="between")
zdata
cor(zdata$zw_centered_beer,zdata$zb_mean_beer)
names(zdata)
model<-lmer(smile~zw_centered_beer+zb_mean_beer+(1+zw_centered_beer|bar),data=zdata)
summary(model)
mean(tapply(data$smile, data$bar, sd))
mean(tapply(zdata$zbeer, data$bar, sd))
sd(zdata$zb_mean_beer)

a<-tapply(data$smile, data$bar, mean)

b<-tapply(data$beer, data$bar, mean)
cor(a,b)

mr<-sapply(unique(data$bar), function(x) {
    d<-data[data$bar==x,c("smile","beer")]
    cor(d)[1,2]
})
mean(mr)

a<-scale(data$smile,scale=T)*5
sd(a)
fixed=c(5.6,.556)
VarCorr<-matrix(c(8.33,-0.1,-0.1,.30),2)
sigma<-1.431
amod<-simr::makeLmer(smile~beer+(1+beer|bar), fixed, VarCorr=VarCorr, sigma,data=data)


lmod<-lmerTest::as_lmerModLmerTest(amod)
summary(amod)

summary(lmod)
a<-modelpower(amod,effects = "random", nsim=250)
a
required_n(amod,expand = "bar",direction = "between",nsim = 250)

modelform<-y~0*1+c(.4,.2)*x+(.4*1|bar)

aa<-as.character(modelstring)
aa

    first<-as.character(modelform)
    rsh<-first[[2]]
    lsh<-first[[3]]
    lsh<-gsub(" ","",lsh)
    lsh
    params<-gsub("*","",stringr::str_extract_all(lsh,paste0(rgx_real,"\\*"),simplify = TRUE),fixed = TRUE)
    params
    termsstring<-gsub("*","",gsub("[0-9]\\*","",lsh))
    form<-as.formula(termsstring)
    terms(form)

clusters<-list(list(name="bar",n=15,type="cross"))
vars<-list(list(name="y",n=30,role="dependent"),
           list(name="x",type="numeric",role="within"),
           list(name="f",type="nominal", n=3,role="within")
           )
amat<-matrix(c(1,.5,.5,1),ncol=2)
wmat<-matrix(c(1,.1,.1,1),ncol=2)

for (w in wmat) print(rbeta(1,1,1))

hist(rbeta(1000,2,2))

a<-runif(10000,.4,.5)
var(a)

wmat<-matrix(c(1,.1,.1,1),ncol=2)

wmat<-matrix(c(1,.1,.1,1),ncol=2)
rdraw(wmat)
NC<-100
r<-.5
v<-2*(1-r)/12
v
rr<-rnorm(NC,r,sqrt(v))
rr[rr>.99999]<-.9999999
rr[rr<-.99999]<--.9999999
summary(rr)
hist(rr)
mean(rr)
ii<-rnorm(NC,0,1)

b<-rnorm(NC,1,1)
data<-list()
mus<-c(0,0)
for (i in seq_len(NC)) {
    sigmaw<-matrix(c(1,rr[i],rr[i],1),ncol=2)
    dd<-MASS::mvrnorm(100,Sigma = sigmaw ,mu=c(ii[i],0),empirical = T)
    colnames(dd)<-c("y","x")
    data[[length(data)+1]]<-cbind(dd,i)
}
orr<-lapply(data, function(d) {
    d<-as.data.frame(d)
    cor(d$y,d$x)
    })
var(unlist(orr))
var(unlist(rr))
mean(unlist(orr))
#data<-.asample(vars,clusters,type="nested",sigmaw = amat,varcorr=wmat)
ddata<-as.data.frame(do.call(rbind,data))
ddata$bar<-factor(ddata$i)
ddata$xx<-ddata$x
ddata$yy<-ddata$y*10

amod<-lmer(yy~xx+(1+xx|bar),data=ddata)
summary(amod)
amod<-lmer(y~x+(1+x|bar/car),data=ddata)
a<-y~x+(1+x|bar/car)



library(lme4)
findbars(a)
summary(amod)


varcorr<-matrix(c(1,.1,.1,1),ncol=2)

clusters<-list(
          list(name="level3",n=10,type="nested"),
          list(name="level2",n=10,type="nested"),
          )
vars<-list(list(name="y",n=30,role="dependent"),
           list(name="x",type="numeric",role="within"),
           list(name="f",type="nominal", n=2,role="level1")
)
sample<-.asample(vars,clusters,type = "nested")
names(sample)
sample$f<-factor(sample$f)

mod<-lmer(y~x*f+(1+x|level2)+(1+f*x|level3),data=sample)
mod

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
