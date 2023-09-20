#' Create a dataset based on a data generating model
#'
#' Create a dataset based on a data generating model passed
#' with a parametrized formula or a set of parameter matrices
#'
#' @param vars  a list defining variables characteristics
#' @param clusters a list defining clusters carateristics
#' @param structure structure of multi-level clusters.
#'                  `nested` (default) or `cross` for cross-classified
#' @param formula a string containing a string formula with terms and parameters
#' @param beta    fixed effect coefficient vector. Overrides parameters in formula.
#' @param varcorr random effect variances-covariances.
#'                 A list of matrices, one for each cluster.
#'                 Overrides parameters in formula.
#'
#' @return an object of class \code{dataframe}
#' @details more to come
#' @author Marcello Gallucci
#' @export

get_sample<-function(vars,clusters,structure,sigmaw=NULL) {

    nclusters<-lapply(clusters, function(x) 1:x$n)
    clusternames<-unlist(rlist::list.select(clusters,name))

    if (structure=="nested") {
      clustervalues<-as.data.frame(.nestedclusters(nclusters))
      names(clustervalues)<-clusternames
    }
    else {
      clustervalues<-expand.grid(nclusters)
    }

    colnames(clustervalues)<-clusternames
    clustervalues$id<-1:dim(clustervalues)[1]
    clusterdata<-list()
    dep<-rlist::list.find(vars, varying=="dependent",n=1)[[1]]
    covs<-rlist::list.find(vars, type=="numeric",n=Inf)
    factors<-rlist::list.find(vars, type=="nominal",n=Inf)

    for (cluster in clusters) {
      ## first, we prepare a frame with the cluster values
      cv<-data.frame(unique(clustervalues[,cluster$name]))
      names(cv)<-cluster$name
      ## we select the continuous between variables
      bvars<-rlist::list.find(covs, grep(cluster$name,varying)>0,n=Inf)
      if (is.something(bvars)) {
        bnames<-unlist(rlist::list.select(bvars,name))
        bn<-length(bnames)
        if (hasName(cluster,"sigma"))
          sigmab<-cluster$sigma
        else {
          sigmab<-matrix(0,ncol=bn,nrow=bn)
          diag(sigmab)<-1
        }
        mus<-rep(0,bn)
        n<-dim(cv)[1]
        ccovs<-as.data.frame(MASS::mvrnorm(n,Sigma = sigmab,mu=mus,empirical = T))
        names(ccovs)<-bnames
        cv<-cbind(cv,ccovs)
      }
      ## here we deal with between factors
      fvars<-rlist::list.find(factors, grep(cluster$name,varying)>0,n=Inf)
      if (is.something(fvars)) {
        nfactors<-lapply(fvars, function(x) 1:x$n)
        fnames<-lapply(fvars, function(x) x$name)
        eg<-.expandto(nfactors,cluster$n)
        colnames(eg)<-paste0(cluster$name,".",fnames)
        cv<-cbind(cv,eg)
      }
      clusterdata[[cluster$name]]<-cv
    }

    # now the within variables
    NC<-dim(clustervalues)[1]
    ## continuous variables

    cvars <- rlist::list.find(covs, varying=="within",n=Inf)
    if (!is.something(cvars))
      stop("No within clusters variable has been defined.")
    wnames<-unlist(rlist::list.select(cvars,name))
    wn<-length(wnames)

    if (is.null(sigmaw))
      sigmaw<-diag(1,ncol=wn,nrow=wn)
    mus<-rep(0,wn)
    ldata<-list()
    for (i in 1:NC) {
      cv<-as.data.frame(MASS::mvrnorm(dep$n,Sigma = sigmaw ,mu=mus,empirical = T))
      names(cv)<-wnames
      cls<-subset(clustervalues,id==i)
      cls[["id"]]<-NULL
      for (cl in names(cls)) {
        cv[,cl]<-cls[[cl]]
        cv<-merge(cv,clusterdata[[cl]],by=cl,all.x = T)
      }
      ldata[[length(ldata)+1]]<-cv


    }

    data<-as.data.frame(do.call(rbind,ldata))
    data$id<-1:dim(data)[1]

    fvars <- rlist::list.find(factors, varying=="within",n=Inf)


    if (is.something(fvars)) {
      nfactors<-lapply(fvars, function(x) 1:x$n)
      fnames<-lapply(fvars, function(x) x$name)
      eg<-.expandto(nfactors,dim(data)[1])
      colnames(eg)<-fnames
      data<-cbind(data,eg)
    }

    finaldata<-data

    for (n in clusternames)
      finaldata[,n]<-factor(finaldata[,n])

    finaldata$id<-1:dim(data)[1]
    finaldata[[dep$name]]<-rnorm(dim(data)[1])
    as.data.frame(finaldata)
  }


#' Create a model with required parameters
#'
#' Create a model with fixed and random effects defined
#' by the user.
#'
#' @param formula model formula. It can be a character vector containing
#'                a string formula with terms and parameters (see details). Alternatively,
#                  can be passed as formula. When passed as \code{\link{formula}},, \code{beta} and \code{varcorr} are required.
#' @param data    a data.frame containing the variables used in the formula.
#' @param beta    fixed effects coefficients vector. Overrides parameters in formula.
#' @param varcorr random effect variances-covariances.
#'                 A list of matrices, one for each cluster.
#'                 Overrides parameters in formula.


make_mer_model<-function(formula,data,beta=NULL,varcorr=NULL) {

  if (inherits(formula,"character")) {
       info<-get_formula_info(formula)
       .formula<-info$formula
        w<-rev(info$random)
       .varcorr<-mats<-lapply(w, function(x) diag(x))
       .beta <- info$fixed
  } else {
       .formula<-formula
       .varcorr<-varcorr
       .beta<-beta
  }
  suppressMessages(model<-lme4::lmer(.formula,data=data))
  simr::sigma(model)<-1
  simr::VarCorr(model)<-.varcorr
  simr::fixef(model)<-.beta
  model@optinfo$conv$lme4$messages<-NULL
  model


}




### helper functions


.expandto<-function(alist,n) {

  eg<-expand.grid(alist)
  res<-list()
  l<-dim(eg)[1]
  j<-0
  for (i in 1:n) {
    j<-j+1
    ladd(res)<-eg[j,]
    if (j==l) j=0
  }
  res<-do.call(rbind,res)
  rownames(res)<-1:dim(res)[1]
  res
}


.levelstandard<-function(data,y,var,cluster,level="within", overwrite=FALSE) {

  if (is.factor(data[[var]])) {

     xdata<-data
     if (level=="between") {
        s0<-tapply(xdata[[y]], list(xdata[[cluster]]), mean)
        mm<-sd(s0)
     } else {
       s0<-tapply(xdata[[y]], list(xdata[[cluster]]), sd)
       mm<-mean(s0)
     }
     codes<-contrasts(xdata[[var]])
     codes<-mm*codes/apply(codes,2,sd)
     print(codes)
     contrasts(xdata[[var]])<-codes


  } else {

  zname<-var
  if (level=="between") {

    if (!overwrite) zname<-paste0("zb_",var)
    xdata<-data
    c0<-round(tapply(data[[var]], list(data[[cluster]]), sd),digits = 2)

    if (any(c0>.001))
        warning("Variable ",var," does not seem to be a between levels variable")

    s0<-tapply(xdata[[y]], list(xdata[[cluster]]), mean)
    mm<-sd(s0)
    xdata[[zname]]<-as.numeric(scale(data[[var]]))*mm
  } else {
     data$id....<-1:dim(data)[1]
     zname<-paste0("zw_",var)
     s0<-tapply(data[[y]], list(data[[cluster]]), sd)
     mm<-mean(s0)
     vdata<-c()
     cls<-unique(data[[cluster]])
     tot<-length(cls)
     i<-0
     for (cl in cls) {
      i<-i+1
      ldata<-data[data[[cluster]]==cl,c("id....",var)]
      ex<-scale(ldata[[var]],scale=mm)
      ldata[[zname]]<-as.numeric(scale(ldata[[var]]))*s0[i]
      ldata[[2]]<-NULL
      vdata<-rbind(ldata,vdata)
    }
    xdata<-merge(data,vdata,by="id....")
    if (overwrite) {
        xdata[[var]]<-xdata[[zname]]
        xdata[[zname]]<-NULL
      }
    xdata$id....<-NULL
   }
  }
  xdata

}

.levelcenter<-function(data,var,cluster,level="within") {


   zname<-paste0("centered_",var)
   m0<-tapply(data[[var]], list(data[[cluster]]), mean,na.rm=TRUE)
   mm<-as.data.frame(cbind(names(m0),m0))
   names(mm)<-c(cluster,".mm.")
   zdata<-merge(data,mm,by=cluster)
   zdata[[".mm."]]<-as.numeric(zdata[[".mm."]])
   if (level=="within") {
      zdata[[zname]]<-zdata[[var]]-zdata[[".mm."]]
    } else {
      zname<-paste0("mean_",var)
      zdata[[zname]]<-zdata[[".mm."]]
    }
    zdata[[".mm."]]<-NULL
    zdata
}


## this produce a frame of codes for the nested clusters
## the mother cluster is 100,200,300 ..., the child
## will be 1001 1002 ... 2001, 2002 ...

.nestedclusters=function(clusters) {

  eg<-rev(expand.grid(rev(clusters)))
  alist<-lapply(1:(ncol(eg)), function(i) {
    a<-eg[,i]*1000^(ncol(eg)-i)
    a
  }  )
  alist
  eg<-do.call(cbind,alist)
  res<-eg
  for (i in 1:nrow(eg))
    for (j in 1:ncol(eg))
      res[i,j]<-sum(eg[i,1:j])/(10^(ncol(eg)-j))

  res
}



.rdraw<-function(sigma,varcorr) {

  ncol<-ncol(varcorr)
  res<-sigma
  for (i in 1:ncol)
    for (j in 1:ncol) {
      if (j==i) next
      s<-sqrt(varcorr[i,j])
      res[i,j]<-rnorm(1,sigma[i,j],s)
    }
  res
}
