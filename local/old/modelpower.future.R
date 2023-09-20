#' Compute power parameters mixed models
#'
#' Compute required N, power level or minimum effect size for a
#' linear mixed model estimated.
#'
#' @param  model of class XXXX
#' @param  n
#' @return an object
#' @details Exactly one of the parameters n,
#'  power, alpha must be passed as NULL, and that parameter
#'  is determined from the others.
#'  Notice that alpha has non-NULL default so NULL must be
#'  explicitly passed if you want it computed.
#' @author Marcello Gallucci
#' @export

modelpower<- function(model,...) UseMethod(".power")

.power.default<-function(model,alpha=.05,nsim=10) {

  print(class(model))
  test  <- ((is.null(n))+(is.null(power))+(is.null(alpha)))
  if (test!=1)
       stop("exactly one of 'n', 'power', and 'alpha' must be NULL ")

}

.power.lmerMod<-function(model,alpha=.05,nsim=10,effects=c("fixed"),seed=NULL) {

  results<- .powerMix(model,nsim = nsim,seed = seed,alpha = alpha, effects=effects)
  results
}


library(gamlj)
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
modelpower(model, nsim=50)
modelpower(model,alpha = .10, nsim=50)

library(simr)

time<-Sys.time()
a<-powerSim(model,nsim=150,test=fixed("x",method="t"))
a
Sys.time()-time

modelpower(ex,alpha=.05,nsim=150)

obj<-modelpower(ex,alpha=.05,nsim=150,effects="fixed")
obj
obj<-modelpower(ex,alpha=.05,nsim=150,effects="random")
obj
obj<-modelpower(ex,alpha=.05,nsim=150)
obj


obj<-.powerMix(model,nsim=10)
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


data<-subjects_by_stimuli



data$id<-1:dim(data)[1]
cluster<-"subj"
afew<-sample(levels(data[[cluster]]),replace = TRUE,size = 200)
ids<-unlist(tapply(data$id, data$subj, function(x) sample(x,size=20)))
sdata<-data[data$id %in% ids,]
sdata<-droplevels(sdata)
table(sdata$subj)
table(sdata$stimulus)
table(sdata$stimulus,sdata$subj)

a<-modelpower(ex,n=10,alpha=.05,nsim=10)
a
print(a)
.effects<-colnames(model.matrix(model))
.model<-.model

library(progress)
pb <- progress_bar$new(total = 100)
for (i in 1:100) {
  pb$tick()
  Sys.sleep(1 / 100)
}

pb <- progress_bar$new(
  format = "  downloading [:bar] :percent time: :elapsed",
  total = 100, clear = FALSE, width= 60)
for (i in 1:100) {
  pb$tick()
  Sys.sleep(1 / 100)
}

str(pb)
pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent",
  total = 30, clear = FALSE, width = 60)
for (i in 1:30) {
  pb$tick()
  Sys.sleep(3 / 100)
}

pb <- progress_bar$new(
  format = "  downloading :what [:bar] :percent eta: :eta",
  clear = FALSE, total = 200, width = 60)
f <- function() {
  for (i in 1:100) {
    pb$tick(tokens = list(what = "foo   "))
    Sys.sleep(2 / 100)
  }
  for (i in 1:100) {
    pb$tick(tokens = list(what = "foobar"))
    Sys.sleep(2 / 100)
  }
}
f()



## From Venables and Ripley (2002) p.165.
N <- c(0,1,0,1,1,1,0,0,0,1,1,0,1,1,0,0,1,0,1,0,1,1,0,0)
P <- c(1,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,0,0,1,0,1,1,0)
K <- c(1,0,0,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,1,1,0,1,0)
yield <- c(49.5,62.8,46.8,57.0,59.8,58.5,55.5,56.0,62.8,55.8,69.5,
           55.0, 62.0,48.8,45.5,44.2,52.0,51.5,49.8,48.8,57.2,59.0,53.2,56.0)

npk <- data.frame(block = gl(6,4), N = factor(N), P = factor(P),
                  K = factor(K), yield = yield)
replications(~ . - yield, npk)

replications(~ . , data)
