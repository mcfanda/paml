k<-10
clusters<-list(list(name="cluster1",n=k,yvar=10))

nk<-round(runif(k,10,50))
#nk<-20
yk<-round(runif(k,30,50))

vars<-list(list(name="y",n=nk,type="numeric",varying="dependent",yvar=yk),
           list(name="x0",type="numeric",varying="within"),
           list(name="x1",type="numeric",varying="cluster1")
)


form<-"y~[.10]*1+[.5]*x0+[.3]*x1+([1]*1|cluster1)"

set.seed(1000)

library(paml)
sample<-make_sample(vars,clusters,"nested",formula=form)
info<-get_formula_info(form)
model<-lmer(info$formula,data=sample)
betas(model,verbose = T)
betas(model,method = "vars",verbose = T)

sample$x0[2]<-NA
sample$y[10]<-NA
tapply(sample$y,sample$cluster1,length)
tapply(sample$y,sample$cluster1,mean,na.rm=T)

model0<-lmer(y~1+(1|cluster1),data=sample)
summary(model0)
var_i<-as.numeric(VarCorr(model0)$cluster1)


sample$y<-sample$y+1
info<-get_formula_info(form)
model1<-lmer(info$formula,data=sample)
(ff<-fixef(model1))

sd_x0<-mean(tapply(sample$x0,sample$cluster1,sd,na.rm=T))
sd_x1<-sd(tapply(sample$x1,sample$cluster1,mean,na.rm=T))

sd_y0<-mean(tapply(sample$y,sample$cluster1,sd,na.rm=T))
sd_ym<-sd(tapply(sample$y,sample$cluster1,mean,na.rm=T))

ff[2]*sd_x0/sigma(model0)
ff[2]*sd_x0/sd_y0

ff[3]*sd_x1/as.numeric(sqrt(var_i))
ff[3]*sd_x1/sd_ym

cor(tapply(sample$x1,sample$cluster1,mean),tapply(sample$y,sample$cluster1,mean))

summary(model0)

sample_c<-levelcenter(sample,"x1","cluster1",level="between")
sample_c<-levelcenter(sample_c,"y","cluster1",level="between")
sample_c<-levelcenter(sample_c,"y","cluster1",level="within")

sd_x1_s<-sd(sample_c$mean_x1)
sd_ym_s<-sd(sample_c$mean_y)
summary(model0)



ff[3]*sd_x1_s/sd_ym_s
ff[3]*sd_x1_s/sqrt(as.numeric(var_i))


n<-mean(nk)
sig<-sum(sample_c$centered_y^2,na.rm=T)/(n*k-k)
sig

a<-tapply(sample$x1,sample$cluster1,mean)
b<-tapply(sample$y,sample$cluster1,mean)
cor(a,b)
bet<-var(b)
bet
bet-(sig/n)

var_i+(sigma(model0)^2)/n
var(b)
teta2<-as.numeric(var_i/(var_i+(sigma(model0)^2)/n))
summary(model0)
ff[3]*sd_x1_s/sd_ym_s
ff[3]*sd_x1_s

asample<-levelscale(sample,"y","x0","cluster1",level = "within")
asample<-levelscale(asample,"y","x1","cluster1",level = "between")
model_s<-lmer(y~sw_x0+sb_x1+(1|cluster1),data=asample)
summary(model_s)


library(paml)
clusters<-list(list(name="cluster2",n=20,ysd=5),
               list(name="cluster1",n=10,ysd=10))


vars<-list(list(name="y",n=50,type="numeric",varying="dependent"),
           list(name="w0",type="numeric",varying="within"),
           list(name="w1",type="numeric",varying="cluster1"),
           list(name="w2",type="numeric",varying="cluster2")

)


#form<-"y~[.10]*1+[.12]*w1+[.2]*w2+([1]*1+[1]*w1|cluster1)+([4]*1+[2]*w2|cluster2)"



form<-"y~[.10]*1+[.2]*w0+[.3]*w1+[.5]*w2+([1]*1|cluster2)+([1]*1|cluster1)"

#form<-"y~[.10]*1+[.212]*w1+[.3]*w2+([1]*1+[.5]*w0|cluster1)"


library(paml)
sample<-make_sample(vars,clusters,"nested",formula=form)
sample$x<-sample$w0+sample$w1+sample$w2
sample$y<-sample$y+1
info<-get_formula_info(form)
mod0<-lmer(y~(1|cluster1)+(1|cluster2),data=sample)
summary(mod0)



model1<-lmer(info$formula,data=sample)
fixef(model1)

asample<-levelcenter(sample,"w0",cluster = "cluster1",level = "within",overwrite = TRUE)
asample<-levelcenter(asample,"w1",cluster = "cluster2",level = "within",overwrite = TRUE)
asample<-levelcenter(asample,"w2",cluster = "cluster2",level = "between",overwrite = TRUE)
model2<-lmer(info$formula,data=asample)

ff<-fixef(model2)
ff
zapsmall(betas(model2,verbose = T))

x<-tapply(sample$w2,sample$cluster2,mean)
y<-tapply(sample$y,sample$cluster2,mean)
cor(x,y)
#betas(model2,method = "vars",verbose = T)


x<-tapply(asample$w1, asample$cluster1, mean)
y<-tapply(asample$y, asample$cluster1, mean)
cor(x,y)

bsample<-sample
bsample$w0<-bsample$x
bsample$w1<-bsample$x
bsample$w2<-bsample$x

bsample<-levelcenter(bsample,"w0",cluster = "cluster1",level = "within",overwrite = TRUE)
bsample<-levelcenter(bsample,"w1",cluster = "cluster1",level = "between",overwrite = TRUE)
bsample<-levelcenter(bsample,"w1",cluster = "cluster2",level = "within",overwrite = TRUE)
bsample<-levelcenter(bsample,"w2",cluster = "cluster2",level = "between",overwrite = TRUE)

cor(bsample$w0,bsample$w1)
cor(bsample$w0,bsample$w2)

model2<-lmer(info$formula,data=bsample)

ff<-fixef(model2)
zapsmall(ff)
zapsmall(betas(model2,verbose = T))
zapsmall(betas(model2,method="vars",verbose = T))

mean(tapply(bsample$y,bsample$cluster1,sd))
mean(tapply(bsample$y,bsample$cluster2,sd))




###############

form<-"y~[.10]*1+[.5]*w0+[-.5]*w1+[.0]*w2+([1]*1|cluster2)+([1]*1|cluster1)"

sample<-make_sample(vars,clusters,"nested",formula=form)
sample$x<-sample$w0+sample$w1+sample$w2
tapply(sample$x, sample$cluster2, mean)
info<-get_formula_info(form)

model1<-lmer(info$formula,data=sample)
fixef(model1)
betas(model1)
asample<-sample
asample$w0<-asample$x
asample<-levelcenter(asample,"w0",cluster = "cluster1",level = "within",overwrite = TRUE)
asample$w1<-asample$x
asample<-levelcenter(asample,"w1",cluster = "cluster1",level = "between",overwrite = TRUE)
asample<-levelcenter(asample,"w1",cluster = "cluster2",level = "within",overwrite = TRUE)
asample$w2<-asample$x
tapply(asample$w2, asample$cluster2, mean)
asample<-levelcenter(asample,"w2",cluster = "cluster2",level = "between",overwrite = TRUE)
asample$y1<-asample$y
asample<-levelcenter(asample,"y1",cluster = "cluster1",level = "between",overwrite = TRUE)
asample<-levelcenter(asample,"y1",cluster = "cluster2",level = "within",overwrite = TRUE)

var(asample$w2)

model2<-lmer(info$formula,data=asample)

ff<-fixef(model2)
zapsmall(betas(model2,verbose = T))

w<-unlist(lapply(unique(asample$cluster2), function(x) cor(asample$x[asample$cluster2==x],asample$y[asample$cluster2==x]) ))
mean(w,na.rm=T)
