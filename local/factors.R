library(paml)
ck<-9
nk<-6

clusters<-list(list(name="cluster1",n=ck,ysd=12))

vars<-list(list(name="y",n=nk,type="numeric",varying="dependent",ysd=14),
           list(name="a",type="nominal",k=3,varying="cluster1"),
           list(name="b",type="nominal",k=3,varying="within"),
           list(name="x",type="numeric",varying="within")

)
formula<-"y~[.10]*1+[.2]*x+[.2]*b1+[.2]*b2+[.2]*a1+[.3]*a2+([1]*1|cluster1)"


get_formula_info(formula)
data<-make_sample(vars,clusters,"nested",formula = formula)

model<-lmer(y~x+a+b+(1|cluster1),data)
fixef(model)

zapsmall(betas(model,verbose=T))

data<-make_sample(vars,clusters,"nested",
                  formula = formula,
                  factors="contrasts",
                  empirical = TRUE)
head(data)
info<-get_formula_info(formula)
model<-lmer(info$formula,data)
zapsmall(fixef(model))
zapsmall(betas(model))


data<-make_sample(vars,clusters,"nested",formula = formula)
head(data)
data$x<-factor(rep(c("a","b"),nrow(data)/2))
as.numeric(data$x)
model<-lmer(y~a+(1|cluster1),data)
zapsmall(fixef(model))

data$a2[1:6]
a0<-c(0,1,0,0,1,0)
b0<-a0-mean(a0)
sqrt(sum(b0^2)/5)
b0/sqrt(sum(b0^2)/5)

x<-c(0,1,0)
xx<-x-mean(x)
k=3
sum(xx^2)*n/(k*(n-1))

(a<-contrasts(data$a))
mm<-apply(a,2,mean)
b<-a-mm
s<-apply(b,2,function(x) sum(x^2))
den<-sqrt(s*(n/k)*(1/(n-1)))
den
cont<-b/den

data$a<-factor(data$a)
contrasts(data$a)<-cont
model<-lmer(y~a+(1|cluster1),data)
zapsmall(fixef(model))

tapply(data$id, data$cluster1,length)
sqrt(var(as.numeric(data$a)))

data$aa<-scale(as.numeric(data$a))
model<-lmer(y~aa+(1|cluster1),data)

ss<-lapply(data$cluster1,function(id) scale(as.numeric(data$a[data$cluster1==id])))
unlist(ss)


length(data$id)/2
dd<-lapply(data$cluster1,function(id) psych::cohen.d(data$y[data$cluster1==id],data$a[data$cluster1==id])$cohen.d[2])
unlist(dd)

dd<-lapply(data$cluster1,function(id) cor(data$y[data$cluster1==id],as.numeric(data$a[data$cluster1==id])))
unlist(dd)

d<-psych::cohen.d(data$y,data$a)
d$r
r<-.3
(2*r)/sqrt(1-r^2)

zapsmall(betas(model))
y<-sum(tapply(data$y,data$cluster1,sd))
x<-sum(tapply(data$a1,data$cluster1,sd))


### repeated measure


ck<-50
nk<-2

clusters<-list(list(name="cluster1",n=ck,ysd=12))

vars<-list(list(name="y",n=nk,type="numeric",varying="dependent",ysd=14),
           list(name="a",type="nominal",k=2,varying="within")

)
formula<-"y~[.10]*1+[.3]*a1+([1]*1|cluster1)"

get_formula_info(formula)
data<-make_sample(vars,clusters,"nested",formula = formula)

model<-lmer(y~a+(1|cluster1),data)
summary(model)
zapsmall(fixef(model))

zapsmall(betas(model,verbose=T))

n<-2
x<-c(-.994987,.994987)
sd(x)/2
X<-rep(c(-.707,.707),50)
lm(data$within_dep_~as.numeric(data$a))

beta<-rep(runif(50,0,.4),each=2)
Y<-X *  beta
cor(X,Y)
col<-1
.X<-as.matrix(X,ncol=col,nrow=n)/sqrt(n-1)
e<-scale(rnorm(100))

Y<- sqrt(beta^2/(1-beta^2))* X + e
cor(X,Y)

b<-c(-0.1414,.1414)
sum(b*a)
