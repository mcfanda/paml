library(lmerTest)
library(simr)
library(paml)
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


sigma1<-matrix(c(2,0,0,1),ncol=2)

clusters<-list(list(name="cluster",n=20),
               list(name="cluster2",n=20))

head(data)
vars<-list(list(name="y",n=5,varying="dependent"),
           list(name="w",type="numeric",varying="within"),
           list(name="x",type="numeric",varying="cluster")
)

form<-"y~[.10]*1+[.12]*w+[.2]*x+([1]*1+[1]*w|cluster)+([3]*1+[2]*x|cluster2)"

sigmab<-matrix(c(1,0,0,1),ncol=2)
n<-20
aa<-MASS::mvrnorm(n,Sigma = sigmab,mu=c(0,0),empirical = T)
apply(aa,2,scale)

data<-get_sample(vars,clusters,structure="nested",formula = form)
head(data)
info<-get_formula_info(form)
model<-lmer(info$formula,data=data)
VarCorr(model)
w<-info$random
mats<-lapply(w, function(x) diag(x))
VarCorr(model)<-mats
fixef(model)<-info$fixed

data$y<-model.matrix(~x+w,data) %*% c(3,.3,.5)
rr<-mkReTrms(findbars(y~w+(1+w|cluster)+(1+x|cluster2)),data)
dim(rr$Zt)
tapply(data$w, data$cluster, mean)
rownames(rr$Zt)


get_formula_info(form)

mmfixed<-model.matrix(nobars(res$formula),data=data)
vv<-findbars(res$formula)
l<-length(vv[[1]][[2]])



mmrandom<-mkReTrms(findbars(y~1+(1+x|cluster)+(0+x|cluster2)),data)
rownames(mmrandom$Zt)
mmrandom$cnms
rr$Zt
data$y<- as.numeric(m %*% rr$Zt )
data$y<-data$y+rnorm(length(data$y))
var(tapply(data$y,data$cluster,mean))
lmer(y~w+(1+w|cluster),data)
qq<-as.data.frame(cbind(y,data$cluster))
var(tapply(qq$y, data$cluster, mean))
subbars(y~x+(1|cluster))

data("Pixel", package="nlme")
mform <- pixel ~ day + I(day^2) + (day | Dog) + (1 | Side/Dog)
(bar.f <- findbars(mform)) # list with 3 terms
mf <- model.frame(subbars(mform),data=Pixel)
rt <- mkReTrms(bar.f,mf)
rt$Zt
rownames(rt$Zt)

data$y<-data$y+3
var(tapply(data$y, data$cluster, mean))
smodel<-lmer(y~x+(1+x|cluster),data=data)
beta<-c(1,.5)
theta<-c(1,.5,1)
varcorr<-matrix(c(1,.1,.1,1),ncol=2)
theta<-calcTheta(varcorr)
sdata<-data
check<-1

sdata$y<-unlist(simulate(~x+(1+x|cluster),newdata=data,
                newparams = list(theta = theta, beta = beta,sigma=1)))

smodel<-lmer(y~x+(1+x|cluster),data=sdata)

summary(smodel)
smodel@beta<-beta
summary(as_lmerModLmerTest(smodel))

vv<-VarCorr(smodel)
vv$cluster
smodel@theta
ssmodel<-smodel
VarCorr(ssmodel)<-VarCorr(ssmodel)
VarCorr(ssmodel)<-matrix(vec_to_mat(c(1.07,.74,1.55)),ncol=2)
ssmodel@theta<-c(1.07,.74,1.55)
smodel@theta
VarCorr(smodel)$cluster
ssmodel@u
fixef(ssmodel)<-beta
summary(ssmodel)
summary(as_lmerModLmerTest(ssmodel))

nc<-2
ncseq <- seq_along(nc)
theta<-ssmodel@theta
thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
Li <- diag(nrow = nc)
Li[lower.tri(Li, diag = TRUE)] <- thl[[1]]
Li
spli
sc=sigma(ssmodel)
a<-tcrossprod(sc * Li)

a/sc
VarCorr(smodel)$cluster
ssmodel@theta
x<-diag(1,2)
matrix(ssmodel@theta,ncol=2)
anew<-simr::makeLmer(y~x+(1+x|cluster),fixef = ssmodel@beta, VarCorr=varcorr, sigma=sigma(ssmodel),data=data)
summary(anew)
lmerTest::calcSatterth(anew,c(0,1))
lmerTest:::contest.lmerMod
lmerTest:::contest.lmerModLmerTest
lmerTest:::contestMD.lmerModLmerTest
summary(as_lmerModLmerTest(anew))
ssmodel@vcov_varpar
aanew<-as_lmerModLmerTest(anew)
aanew@vcov_varpar
as_lmerModLmerTest
lmerTest:::as_lmerModLT
ssmodel@theta
varcorr<-matrix(c(1,.09,.09,1.1),ncol=2)
VarCorr(ssmodel)<-VarCorr(smodel)
summary(ssmodel)
ssmodel@theta
smodel@theta
summary(as_lmerModLmerTest(ssmodel))
ss
ss$coefficients

k=6


calcTheta1 <- function(V, sigma=1) {
  L <- suppressWarnings(chol(V, pivot=TRUE))
  p <- order(attr(L, "pivot"))
  L <- t(L[p, p])

  L[lower.tri(L, diag=TRUE)] / sigma
}

# All the thetas
calcTheta <- function(V, sigma) {

  if(missing(sigma)) sigma <- attr(V, "sc")
  if(is.null(sigma)) sigma <- 1

  if(!is.list(V)) V <- list(V)

  theta <- plyr::llply(V, calcTheta1, sigma)

  unname(unlist(theta))
}

s<-1
value<-matrix(c(1,.5,.5,1),ncol=2)

calcTheta(value, s)

x<-"
level: within
x ~ 3*y1 + .4*y2 + .4*y3
level: between
x~~4*y1
x~4*1
"
lavaan::lavaanify(x)


