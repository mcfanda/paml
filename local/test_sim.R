
library(paml)
library(lmerTest)
k<-100
nk<-21
cysd<-1
ysd<-runif(k,1,15)
ysd<-1
clusters<-list(list(name="cluster1",n=k,ysd=cysd))
vars<-list(list(name="y",n=nk,type="numeric",varying="dependent",ysd=ysd),
           list(name="c",type="factor",k=3,varying="within"),
           list(name="x1",type="numeric",varying="cluster1")

)
### declare the poulation model with expected parameters ###

form<-"y~[.10]*1+[.5]*c1+[.2]*c2+[.4]*x1+([1]*1+[1]*x1|cluster1)"

info<-get_formula_info(form)
info
### generate the data ###
library(paml)
sample<-make_model(vars,clusters,"nested",formula=form,empirical = T)
modelpower(sample,.05,method = "sim")
mod<-lmer(y~c1+c2+(1|cluster1),data=sample)
summary(mod)
betas(mod,verbose=T)
mod<-lmer(y~c+(1|cluster1),data=sample)
fixef(mod)[2]
50/12
betas(mod,verbose=T)

mm<-as.data.frame(model.matrix(~c,data=sample))

sd(mm$c2)
sd(mm$c3)


a<-contrasts(sample$c)
k<-nrow(a)
N<-length(sample$c)
val<-sqrt(apply(a-mean(a),2,function(x) sum(x^2))*N/(k*(N-1)))
val


mod<-lmer(y~1+(1|cluster1),data=sample)
summary(mod)

y0<-levelcenter(sample$y,sample$cluster1)
y1<-levelcenter(sample$y,sample$cluster1,level="between")
N<-length(y0)-1
var(y1)-var(y0)/nk
var(sample$cluster1_dep_)
var(sample$within_dep_)

psych::cohen.d(y0,sample$c)
cor(y0,sample$c1)

x<-levelcenter(sample$c1,sample$cluster1,"within")
y<-levelcenter(sample$y,sample$cluster1,"within")
m<-levelcenter(sample$y,sample$cluster1,"between")

cor(sample$within_dep_,sample$y)
cor(sample$y,y)
cor(x,y)
sd(y)
sd(x)
dim(sample)
x
sd(x)
a<-x[1:2]
cl<-sqrt(19/20)*x/sqrt(sum(a^2))

sd(cl)
cor(x,y)

rr<-lapply(sample$cluster1, function(x) cor(sample$y[sample$cluster1==x],sample$c1[sample$cluster1==x]))
rr
mean(unlist(rr))
y<-levelcenter(sample$y,sample$cluster1)
cor(y,sample$c1)

tapply(sample$c,sample$cluster1,sd)
tapply(y,sample$cluster1,sd)
tapply(sample$within_dep_,sample$cluster1,sd)

x<-tapply(sample$c,sample$cluster1,sd)

x<-tapply(sample$x1,sample$cluster1,mean)
y<-tapply(sample$y,sample$cluster1,mean)
cor(x,y)

info$formula
mod<-lmer(y~c+x1+(1|cluster1),data=sample)
summary(mod)
zapsmall(fixef(mod))

object<-simr::makeLmer(info$formula,fixef = c(1,.3,.3),
                       VarCorr = list(c(10)),
                       data=sample,
                       sigma = .85)

y<-.make_y(object)
sample$yy<-y*1
mod<-lmer(yy~c+x1+(1|cluster1),data=sample)
summary(mod)

mod0<-lmer(yy~1+(1|cluster1),data=sample)
summary(mod0)


one<-function() {
y<-.make_y(object)
sample$vv<-y
sample$yy<-y*1
rr<-lapply(sample$cluster1, function(x) cor(sample$yy[sample$cluster1==x],sample$c[sample$cluster1==x]))
mean(unlist(rr))
}
mean(replicate(10,one()))

one<-function() {
  y<-.make_y(object)
  sample$yy<-y
  y<-tapply(sample$yy,sample$cluster1,mean)
  x<-tapply(sample$x1,sample$cluster1,mean)
  cor(y,x)
}

mean(replicate(100,one()))

one<-function() {
  y<-.make_y(object)
  sample$yy<-y*1
  mod<-lmer(yy~c+x1+(1|cluster1),data=sample)
  #fixef(mod)[3]
  vv<-VarCorr(mod)
  as.numeric(vv$cluster1)
}
rr<-replicate(100,one())
sd(rr)
mean(rr)
33/24
exp(3+4)
exp(3)*exp(4)

n<-3
s<- matrix(c(1,0,0,1),ncol=2)
cc<-as.data.frame(MASS::mvrnorm(n,Sigma = s,mu=c(0,02),empirical = T))
cor(cc)
