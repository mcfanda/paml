library(paml)
library(lme4)

k1<-50
k2<-10
#nk<-round(runif(k,10,20))
nk<-50
cysd<-1
#ysd<-round(runif(k,10,20))
ysd<-1
clusters<-list(list(name="cluster1",n=k,ysd=cysd),
               list(name="cluster2",n=k,ysd=cysd))

vars<-list(list(name="y",n=nk,type="numeric",varying="dependent",ysd=ysd),
           list(name="c1",type="numeric",varying="cluster1"),
           list(name="c2",type="numeric",varying="cluster1"),
           list(name="x",type="numeric",varying="cluster2"),
           list(name="x2",type="numeric",varying="within")


)
### declare the poulation model with expected parameters ###

form<-"y~[.10]*1+[.4]*c1+[.3]*c2+[.4]*x2+[.3]*x+([1]*1|cluster1)+([.75]*1+[.6]*x|cluster2)"


info<-get_formula_info(form)

left<-1-sum(info$fixed[c(2,3)]^2)





info$vars<-vars
info$clusters<-clusters
### generate the data ###
one<-function() {
    summary(model)
    sample<-get_sample(model)
    mod<-lmer(info$formula,data=sample)
    m<-VarCorr(mod)$cluster1
    a<-sqrt(sum(m[lower.tri(m,T)]))
    b<-sd(levelcenter(sample$c1,sample$cluster1,level="between"))
   # b<-sd(as.numeric(sample$c1))
    c(fixef(mod),betas(mod,verbose=F),fixef(mod)[2]*b/a)
}

a<-model@frame


model<-make_model(vars,clusters,"nested",formula=form,empirical = T)
model
VarCorr(model)$cluster1
dd<-model@frame
sd(dd$c1)

v<-info$random$cluster1
sum(diag(v))
summary(model)
sample<-get_sample(model)

mod<-lmer(info$formula,data=sample)
summary(mod)
y0<-levelcenter(sample$y,sample$cluster1)
sd(y0)
sd(sample$y)

m<-VarCorr(mod)$cluster1
sqrt(sum(m[lower.tri(m,T)]))
mod2<-lmer(y~c1+x+(1|cluster1),data=sample)
m<-VarCorr(mod2)$cluster1
sqrt(sum(m[lower.tri(m,T)]))

summary(mod2)
summary(mod)

sample<-get_sample(model)
mod0<-lmer(y~(1|cluster1),data=sample)
summary(mod0)
c(fixef(mod),betas(mod,verbose=T))

doFuture::registerDoFuture()
future::plan(future::multisession)
nsim<-500
reslist  <- foreach::foreach(1:nsim) %dorng% one()
apply(do.call(rbind,reslist),2,mean)


two<-function() {
sample<-get_sample(model)
mod0<-lmer(y~(1|cluster1)+(1|cluster2),data=sample)
y0<-var(levelcenter(sample$y,sample$cluster1,"within"))
y1<-var(levelcenter(sample$y,sample$cluster1,"between"))
c(as.numeric(VarCorr(mod0)),sigma(mod0)^2,var(sample$y),y0,y1)
}

sum(info$fixed[c(4,5)]^2)+1+.6
sum(info$fixed[c(2,3)]^2)+1
VarCorr(model)$cluster1

doFuture::registerDoFuture()
future::plan(future::multisession)
nsim<-500
reslist  <- foreach::foreach(1:nsim) %dorng% two()
apply(do.call(rbind,reslist),2,mean)



mod<-lmer(info$formula,sample)
y1<-levelcenter(sample$y,sample$cluster1,level = "between")
sd(y1)
mean(tapply(sample$y, sample$cluster1, sd))
betas(mod,verbose=T)


library(lme4)

cc<-contrasts(sample$c)
sqrt(N*sum((cc-mean(cc))^2)/(2*(N-1)))
sd(as.numeric(sample$c))
summary(model)
betas(model,verbose=T)
bb<-betas(model,verbose=F)[2]
sd(as.numeric(sample$c))
N<-length(sample$c)
cc<-contrasts(sample$c)
ss<-sqrt(sum((cc-mean(cc))^2)*N/(2*(N-1)))

(d0<-bb/(ss*sqrt(1-bb^2)))
ff<-fixef(model)
cor(sample$y,sample$c1)
sample$y0<-levelcenter(sample$y,sample$cluster1)
sample$y1<-levelcenter(sample$y,sample$cluster1,level="between")

cor(sample$y1,sample$c1)

a<-tapply(sample$y,sample$cluster1,mean)
b<-tapply(sample$c1,sample$cluster1,mean)
m<-tapply(a, b, mean)
s<-tapply(a, b, sd)
dm<-m[1]-m[2]
2*(dm)/(s[1]+s[2])

psych::cohen.d(sample$y,sample$c1)

m<-tapply(sample$y, sample$c, mean)
s<-tapply(sample$y, sample$c, sd)
dm<-m[1]-m[2]
2*(dm)/(s[1]+s[2])
sp<-as.numeric(VarCorr(model)$cluster1)
se<-sigma(model)^2
dm/sqrt(sp+se)
dm/sqrt(se)

(d0<-sqrt( (N-k)/(N)) *bb/(ss*sqrt(1-bb^2)))

