# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(future.apply)
library(gamlj)
library(paml)
data("subjects_by_stimuli")
subjects_by_stimuli$cond<-factor(subjects_by_stimuli$cond)
subjects_by_stimuli$subj<-factor(subjects_by_stimuli$subj)
subjects_by_stimuli$stimulus<-factor(subjects_by_stimuli$stimulus)

subjects_by_stimuli$x<-rnorm(nrow(subjects_by_stimuli))

formula<-y~cond*x+(1+x|subj)+(1|stimulus)
model<-lme4::lmer(formula,data=subjects_by_stimuli)
formula0<-y~cond*x+(1|subj)
model0<-lme4::lmer(formula0,data=subjects_by_stimuli)
lmerTest::ranova(model0,reduce.terms = F)

extend

nsim<-10
effects<-"fixed"
alpha<-.05
names(subjects_by_stimuli)
required_n(model,expand = "subj",direction = "between", target="cond1",nsim=500, seed=100)
obj<-modelpower(model, nsim=10)
methods(class="pamlpowermix")

.evaluatePower(obj, .80)




modelpower(model,alpha = .10, nsim=50)
library(doRNG)
library(simr)
time<-Sys.time()
a<-powerSim(model,nsim=150,test=fixed("x",method="t"))
a
Sys.time()-time
methods(modelpower)
obj<-modelpower(model,alpha=.05,nsim=150)
obj<-modelpower(model,alpha=.05,nsim=9)
obj

obj<-modelpower(ex,alpha=.05,nsim=150,effects="fixed")
obj
obj<-modelpower(ex,alpha=.05,nsim=150,effects="random")
obj
obj<-modelpower(ex,alpha=.05,nsim=150)
obj


obj<-.powerMix(model,nsim=100)
obj
VarCorr(model)
a<-lmerTest::ranova(model)
rownames(a)
simr:::doSim.merMod
table(model@frame$subj)
ex<-extend(model,along="subj", n=60)
obj<-.powerMix(ex,nsim=30)

str(model)
obj<-.powerMix(model,nsim=100)
obj
df<-.evaluatePower(obj,.80)
df
n<-clusters$levels[[cluster]]
n<-clusters$numerosity[[cluster]]
target="cond1:x"
keep=TRUE
while (keep) {
  n<-n+10
  newmodel<-extend(model,within=cluster,n=n)
  tn<-dim(attr(newmodel,"newData"))[1]
  cat("\n Evaluating for",cluster,"with ",n,"levels. Total sample N:",tn,"\n")
  obj<-.powerMix(newmodel,nsim=50)
  df<-.evaluatePower(obj,.80)
  if (is.null(target))
    keep<-any(df$power>power)
  else
    keep<-df$power[rownames(df)==target]<power
  print(obj)
  print(df)
  print(keep)
}

newmodel<-extend(model,within=cluster,n=1000)
obj<-.powerMix(newmodel,nsim=50)
obj


if (length(oks)==length(df$power))
  print("done")
clusters$levels[[cluster]]+10
newmodel<-extend(model,along=cluster,n=clusters$levels[[cluster]]+10)

obj<-.powerMix(newmodel,nsim=50)
obj


data<-subjects_by_stimuli
clusters<-list()
clusters$levels<-lme4::getME(model,"l_i")
clusters$names<-names(clusters$levels)
clusters$numerosity<-sapply(names(obj$clusters),function(x) round(mean(table(data[[x]]))))
clusters

nobj<-simr::extend(model,within = "subj",n = 200)
nd<-attr(nobj,"newData")
nd<-droplevels(nd)

update(model,data=nd)
table(nd$subj)
table(nd$stimulus)


nobj<-simr::extend(model,along = "subj",n = 60)
nd<-attr(nobj,"newData")
nd<-droplevels(nd)

update(model,data=nd)
table(nd$subj)
table(nd$stimulus)



table(nd$stimulus,nd$subj)
table(data$stimulus,data$subj)
o<-data[data$subj==1,]
table(o$stimulus)
dim(data)


data$id<-1:dim(data)[1]
cluster<-"subj"
afew<-sample(levels(data[[cluster]]),replace = TRUE,size = 200)
ids<-unlist(tapply(data$id, data$subj, function(x) sample(x,size=20)))
sdata<-data[data$id %in% ids,]
sdata<-droplevels(sdata)
table(sdata$subj)
table(sdata$stimulus)
table(sdata$stimulus,sdata$subj)

doFuture::registerDoFuture()
future::plan(future::multisession)
future::plan(future::sequential)

afff<-function() {
a=4
afun<-function(model) {
  Sys.sleep(.5)
  return(a)
}
foreach::foreach(1:5) %dorng% afun(model)
}
time<-Sys.time()
a<-afff()
Sys.time()-time

paml::required_n()

a<-settings::options_manager(foo=1,bar='a')
a()

library(simr)
fm <- lmer(y ~ x + (1|g), data=simdata)
fmx1 <- extend(fm, along="x", n=20)
nrow(getData(fmx1))
fmx2 <- extend(fm, along="x", values=c(1,2,4,8,16))
nrow(getData(fmx2))
table(simdata$g)
fmx1 <- extend(fm, along="g", n=2)
nrow(getData(fmx2))
