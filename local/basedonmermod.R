k<-100
#nk<-round(runif(k,10,20))
nk<-16
cysd<-1
#ysd<-round(runif(k,10,20))
ysd<-1
clusters<-list(list(name="cluster1",n=k,ysd=cysd))
vars<-list(list(name="y",n=nk,type="numeric",varying="dependent",ysd=ysd),
           list(name="c",type="factor",k=2,varying="cluster1"),
           list(name="x",type="numeric",varying="within")

)

form<-"y~[.10]*1+[.4]*c1+[.2]*x+([1]*1|cluster1)"

info<-get_formula_info(form)

### generate the data ###
sample<-make_sample(vars,clusters,"nested",formula=form,empirical = T)
model<-lmer(y~c1+x+(1|cluster1),data=sample)
simr::fixef(model)<-info$fixed
simr::VarCorr(model)<-info$random
summary(model)
sample$y<-simulate(model)[,1]
head(sample$y)
model<-lmer(y~c1+x+(1|cluster1),data=sample)

betas(model)
