

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
form<-"y~[1]*1+[3]*cond+(1|subj)+(1|stimulus)"
info<-get_formula_info(form)
info
data<-exdata
.clusters<-list(subj=list(name="subj",layer=2))
for (x in unlist(info$clusters)) {

  if (!is.factor(data[[x]])) {
    data[[x]]<-factor(data[[x]])
    warning("Cluster variable ",x," was coerced to factor.")
  }
  n<-nlevels(data[[x]])
  if (x %in% names(.clusters)) {
    if (hasName(.clusters[[x]],"size"))
       message("Cluster variable ",x," size is set to its size in the the data.frame")
    .clusters[[x]]$size<-n

  } else {
    .clusters[[x]]<-list(name=x, size=n)
  }
}

.expand<-function(...) {

  args<-list(...)
  .data<-data
  for (x in names(args)) {
    if (x %in% names(.clusters)) {
      newdata<-list()
      newid<-1
      obj<-.clusters[[x]]
      n<-num(args[[x]])
      print(dim(.data))
      ids<-as.character(unique(.data[[obj$name]]))
      ids<-sample(ids,n, replace = T)
      for (i in ids) {
        local<-.data[as.character(.data[[x]])==i,]
        local[[x]]<-newid
        newid<-newid+1
        newdata[[length(newdata)+1]]<-local
      }

    } else
       stop("Variable ",x," not in data.frame.")

    .data<-as.data.frame(do.call(rbind,newdata))
    .data[[x]]<-factor(.data[[x]])
  }
  .data
}

.clusters
dd<-.expand(subj=24, stimulus=13)
