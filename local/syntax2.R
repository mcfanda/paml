library(rlist)
x <- list()
library(lmerTest)
library(GAMLj3)
names(clustermanymodels)
model<-lmer(ycont~x+(1+x|cluster),data=clustermanymodels)
ww<-VarCorr(model)
ww$cluster

list.find(x, type==as.character(3), 1)
data<-data.frame(z=1:4,y=1:4,x=1,gr=c(1,1,2,2))
data$y<-factor(data$y)
levels<-1:4
nLevels<-length(levels)
dummy <- stats::contr.treatment(nLevels)
coding <- matrix(rep(1/nLevels, prod(dim(dummy))), ncol=nLevels-1)
contrast <- (dummy - coding)
colnames(contrast)<-paste0(".",1:3)
contrasts(data$y)<-contrast
model.matrix(z~y,data)
lm(z~y,data)

lme4::mkDataTemplate(y~x+(1|gr),data=data,nGrps = 10,nPerGrp = 1,dd)
a<-"afac"
b<-c("afac.1","xafac.2")

grep(paste0("^",a,"\\.[0-9]*$"),b)
form<-y~1+afac2+(1+afac2|endo_id)
formula.tools::split_terms(formula.tools::rhs(form),recursive = T)
formula.tools::rhs(form)
Reduce(form,deparse(form))
paste(deparse(form, width.cutoff = 500), collapse="")


astring<-"y~[1]*1+[2,2.5]*afac2+x1+(0+[4,4.5]*afac2|endo_id)+([5]*1|endo_id)+([6]*1+[7]*x1|video_id)"

      form0<-stringr::str_split_fixed(astring,"~",n=Inf)
      rhs<-form0[[2]]
      clean<-.make_raw(rhs)
      fixed<-stringr::str_split(stringr::str_remove_all(clean,"\\((.*?)\\)"),"\\+")[[1]]
      fixed<-fixed[unlist(sapply(fixed,function(x) x!=""))]
      lapply(fixed, function(x) {
        coef<-as.numeric(stringr::str_remove_all(stringr::str_extract_all(x,"\\[(.*?)\\]"),"[\\[\\]]"))
        terms<-gsub("\\*","",stringr::str_remove_all(x,"\\[(.*?)\\]"))
        print(terms)
        print(coef)
      })

      ### first we deal with the random component
      random<-stringr::str_extract_all(clean,"\\((.*?)\\)")[[1]]
      clusters<-list()
      params<-lapply(as.list(random), function(x) {
              str<-stringr::str_split(x,"\\|",n=Inf)[[1]]
              cl<-gsub("\\)","",str[[length(str)]])
              clname<-cl
              if (cl %in% clusters)
                    clname<-paste(cl,sum(as.numeric(clusters==cl)),sep=".")
              clusters[[length(clusters)+1]]<<-clname
              str<-stringr::str_split(x,"\\+",n=Inf)[[1]]
              suppressWarnings(
                    coefs<-sapply(str, function(z) as.numeric(stringr::str_remove_all(stringr::str_extract_all(z,"\\[(.*?)\\]"),"[\\[\\]]")))
              )
              coefs<-coefs[!is.na(coefs)]
              str<-names(coefs)
              terms<-gsub("\\[(.*?)\\]","",str)
              terms<-gsub("[\\(\\)\\*]","",terms)
              terms<-gsub(paste0("\\|",cl),"",terms)
              terms<-stringr::str_replace_all(string = terms,pattern = "^1","(Intercept)")
              names(coefs)<-paste(clname,terms,sep=".")
              coefs
      })
params
