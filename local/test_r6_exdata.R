library(data.table)



library(GAMLj3)
exdata<-subjects_by_stimuli
head(exdata)
table(exdata$cond)
table(exdata$subj)
table(exdata$cond,exdata$subj)
table(exdata$subj,exdata$stimulus)
names(exdata)
dim(exdata)
exdata$cond<-factor(exdata$cond)
exdata$subj<-factor(exdata$subj)
exdata$stimulus<-factor(exdata$stimulus)

library(paml)
d<-Design$new()
d$cluster<-c(name="subj", layer=2, size=20)
d$cluster<-c(name="stimulus", layer=2, size=30)
d$cluster<-c(name="within", layer=1, size=10)

d$variable<-c(name="cond",type="factor",layer=2,levels=2)
d$variable<-c(name="x",type="numeric",layer=2)
d$variable<-c(name="2",type="numeric",layer=1)

#d$variables<-c(name="x3",type="numeric",layer=3,mean=100)
#d$variables<-c(name="y",type="numeric", dependent=T)
#data<-d$size(endo_id=10,video_id=10, within=20)
form<-"y~[1]*1+[3]*cond+([3]*1+[2]*cond|subj)+([5]*1|stimulus)"
#form<-"y~[1]*1+[3]*cond+([3]*1+[2]*cond|subj)+([5]*1|subj/stimulus)"

s_data<-function() structure(Design$new(),class="sim_data")

d<-s_data()

d<-d+s_integer("x") + s_factor("z",levels=3)
d <- d + s_factor("z")
d$print()

d<-s_data()

d<- d + s_integer("x") + s_factor("z",levels=3)+(s_factor("z",levels=3)|s_cluster("q",levels=10))
d<- d + s_integer("x") * s_factor("z",levels=3)
d<- d + s_dependent("y")

form<-y ~ .p(1) * x + .p(3) * q

form
d$print()
d$print()
d
d$print()

s_integer("x")+s_integer("x")



ggplot2:::`+.gg`()
school/subj:stimuli

#form<-"y~1+ [.25]*afac2+(1+afac2|endo_id)+([.3]*1+[.25]*afac2|video_id)"
#form<-"y~[1]*1+[0,0]*afac2+[0]*x1+([.2]*1+[.2]*x1|endo_id)"
d$formula<-form
print(d$create_model())
zdata<-d$one_sample()
zdata[zdata$subj==1,]
table(zdata$stimulus,zdata$subj,zdata$cond)
zdata$layer_3<-factor(zdata$layer_3)
zdata$layer_2<-factor(zdata$layer_2)
library(lmerTest)
lmer(y~1+(1+cond|layer_3/layer_2), data=zdata)
names(zdata)
zdata$layer_3<-NULL

write.csv(zdata,"local/sub_on_stimuli.csv",row.names = F)
head(d$data)
d$info
d$size(subj=10)
dd<-d$size(cluster=10)
head(dd)

table(dd$cluster)
table(d$data$cluster)

table(d$data$afac3,d$data$cluster)
table(d$data$afac3,d$data$afac2)


u<-unique(d$data$big)
uu<-sample(u,15,replace = T)
acN<-length(u)
reN<-30
ceiling(reN/acN)
nu<-rep(u,3)
nu
K<-2
dt<-data.table(d$data)
dim(dt[on=big==x,])

fun<-function() {
  n<-length(a)
  sig<-matrix(0,ncol = K,nrow = K)
  diag(sig)<-rep(1,K)
  m<-MASS::mvrnorm(n,mu = rep(1,K),Sigma = sig,empirical = T)
  ly(1:K,function(x) m[,x])
}
q<-list(c = 1 , d = 2)
DT[, (names(q)):=fun()]
DT
cols.d0 = paste0("d0.", 1:K)
data[,cols.d0]<-matrix(0,ncol = K,nrow = dim(data)[1])
DT<-data.table(data)
DT$layer._.1<-paste(data$endo_id,data$video_id,sep=".")
names(DT)

qq<-DT[, cols.d0 := lapply(.SD, function(x) fun(x)), keyby=layer._.1,.SDcols = cols.d0][]
names(qq)
qq$d0.1<-as.numeric(qq$d0.1)
qq[, sd(d0.1),by=layer._.1]
qq[, sd(x1),by=layer._.1]
str(qq)
d$options(parallel="multisession")

#p1<-.7
#p2<-p1+.05
#odd<-(p2/(1-p2))/(p1/(1-p1))
#@cat("Odd: ",odd," slope variance:",log(odd),"\n")

#ss<-d$simulate(Nsim=5,endo_id=c(20,25),video_id=10,fixed=list("(Intercept)"=.3,afac2.1=c(.4,.5)),random=list(afac2.1=2))
#ss<-d$simulate(Nsim=5,endo_id=20,video_id=20,random=list('endo_id.(Intercept)'=.3,'endo_id.x1'=.3))

ss
#ss


DT = data.table(a = LETTERS[c(3L,1:3)], b = 4:7)
DT[, d := -8]
DT[b > 4][, b := d * 2L]

DT[2, d := 10L][]


DT_1 = data.table(x=rep(c("b","a","c"),each=3), v=c(1,1,1,2,2,1,1,2,2), y=c(1,3,6), a=1:9, b=9:1)
DT_2 = data.table(x=c("c","b"), v=8:7, foo=c(4,2))
DT_3 = data.table(x=c("c","d"), v=8:7, foo=c(4,2))

X

DT[.N]                                 # last row, only special symbol allowed in 'i'
DT[, .N]                               # total number of rows in DT
DT[, .N, by=x]                         # number of rows in each group
DT[, .SD, .SDcols=x:y]                 # select columns 'x' and 'y'
DT[, .SD[1]]                           # first row of all columns
DT[, .SD[1], by=x]                     # first row of 'y' and 'v' for each group in 'x'
DT[, c(.N, lapply(.SD, sum)), by=x]    # get rows *and* sum columns 'v' and 'y' by group
DT[, .I[1], by=x]                      # row number in DT corresponding to each group
DT[, .N, by=rleid(v)]                  # get count of consecutive runs of 'v'
DT[, c(.(y=max(y)), lapply(.SD, min)),
   by=rleid(v), .SDcols=v:b]      # compute 'j' for each consecutive runs of 'v'
DT[, grp := .GRP, by=x]                # add a group counter
DT[, grp_pct := .GRP/.NGRP, by=x]      # add a group "progress" counter
X[, DT[.BY, y, on="x"], by=x]          # join within each group
data[]

a<-DT_2[DT_1, on = c("x")]
DT_3[a, on=c("x")]

simr::extend()

