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

get_sample<-function(model) {

  data<-model@frame
  dep<-formula(model)[[2]]
  data[[dep]]<-as.numeric(simulate(model)[,1])
  data
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
#' @export

make_model<-function(vars,clusters,
                      structure,
                      formula=NULL,
                      fixed=NULL,
                      random=NULL,
                      factors="variable",
                      empirical=TRUE) {

  if (is.null(vars)) stop("vars argument is required")
  dep<-FALSE
  within<-FALSE
  for (var in vars) {
    if (!hasName(var,"name")) stop("Variable definition requires a `name` option")
    if (!hasName(var,"type")) stop("Variable ",var$name," requires a `type` option")
    if (!(var$type %in% c("numeric","factor"))) stop("Variables can be of type `numeric` or `factor`. Check variable ",var$name)
    if (var$type=="factor")  {
        if (is.null(var$k)) stop("Variable ",var$name," requires to specify number of levels with options `k=`")
        vary<-var$varying
    }
    if (!hasName(var,"varying")) stop("Variable ",var$name," requires a `varying` option")
    if  (var$varying=="dependent") dep<-TRUE
    if  (var$varying=="within") within<-TRUE

  }
  if (isFALSE(dep)) stop("Plese specify a dependent variable with option `varying='dependent'` in variables list")
  if (isFALSE(within)) stop("Plese specify at least one independent variable varying within-clusters")

  if (inherits(formula,"character")) {
       info<-get_formula_info(formula)
       info$vars<-vars
       info$clusters<-clusters
       info$varcorr<-lapply(info$random, function(x) diag(x))
       info$factors<-factors
       info$empirical<-empirical

  } else {
       .formula<-formula
       .varcorr<-varcorr
       .beta<-beta
  }
  info$structure<-structure

  clusters<-info$clusters
  clusters[[length(clusters)+1]]<-list(name="within")
  for (i in 1:length(clusters)) {
    cluster<-clusters[[i]]
    vars<-rlist::list.find(info$vars, varying==cluster$name,n=Inf)
    varsnames<-unlist(rlist::list.select(vars,name))
    betas<-sum(as.numeric(info$fixed[varsnames])^2)
    int<-as.numeric(info$random[[cluster$name]][1])
    slopes<-unlist(lapply(info$random,function(x) x[varsnames]))
    slopes<-sum(slopes[!is.na(slopes)])
    clusters[[i]]$variances<-list(betas=betas,int=int,slopes=slopes)
  }
  print(clusters)
  sample<-.make_sample(info)
  model<-lme4::lmer(info$formula,data=sample)
  simr::fixef(model)<-(info$fixed)
  for (cluster in names(info$random)) {
    info$random[[cluster]]<-diag(info$random[[cluster]])
  }
  simr::VarCorr(model)<-(info$random)
  ### compute sigma ###
  sigma<-1
  ###
  model<-simr::makeLmer(info$formula,info$fixed,info$random,sigma,sample)
  return(model)

}

## make the sample

.make_sample<-function(info) {

  ## gather info for clusters

  clusters<-info$clusters
  nclusters<-lapply(clusters, function(x) 1:x$n)
  clusternames<-unlist(rlist::list.select(clusters,name))

  if (info$structure=="nested") {
    clustervalues<-as.data.frame(.nestedclusters(nclusters))
    names(clustervalues)<-clusternames
  } else {
    clustervalues<-expand.grid(nclusters)
  }

  colnames(clustervalues)<-clusternames
  clustervalues$id<-1:dim(clustervalues)[1]
  clusterdata<-clustervalues

  ## gather info for variables
  vars<-info$vars
  dep<-rlist::list.find(vars, varying=="dependent",n=1)[[1]]
  covs<-rlist::list.find(vars, type=="numeric",n=Inf)
  factors<-rlist::list.find(vars, type=="factor",n=Inf)

  ## create cluster-level data
  for (cluster in clusters) {

    values<-clustervalues[[cluster$name]]
    n<-length(unique(values))
    onedata<-.make_clusterdata(cluster$name,info,n)
    onedata[[cluster$name]]<-unique(values)
    clusterdata<-merge(clusterdata,onedata,by=cluster$name,all.x=T)
  }

  info$clusters<-list(list(name="within",n=dep$n,ysd=dep$ysd))
  wdata<-list()
  nk<-length(unique(clustervalues$id))
  n<-dep$n
  if (length(n)==1) n<-rep(n,nk)
  ucluster<-unique(clustervalues$id)
  ysd<-1
  if (is.something(dep$ysd)) ysd<-dep$ysd
  if (length(ysd)==1) ysd<-rep(ysd,nk)
  for (i in seq_along(ucluster)) {
    onedata<-.make_clusterdata("within",info,n[i],ysd[i])
    onedata[["id"]]<-ucluster[i]
    wdata[[i]]<-onedata
  }
  wdata<-do.call(rbind,wdata)
  finaldata<-merge(clusterdata,wdata,by="id")
  finaldata[[dep$name]]<-runif(dim(finaldata)[1])
  for (cluster in clusternames) finaldata[[cluster]]<-factor(finaldata[[cluster]])
  for (factor in factors)
    if (factor$name %in% names(finaldata))
      finaldata[[factor$name]]<-factor(finaldata[[factor$name]])

  for (name in names(info$variances)) {
#    finaldata[[name]]<-as.numeric(scale(finaldata[[name]]))*info$variances[name]
     finaldata[[name]]<-as.numeric(scale(finaldata[[name]]))

  }


  finaldata
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


#' @export

.levelscale<-function(data,y,var,cluster,level="within", overwrite=FALSE) {

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
     contrasts(xdata[[var]])<-codes


  } else {

  zname<-var
  if (level=="between") {

    if (!overwrite) zname<-paste0("sb_",var)
    xdata<-data
    c0<-round(tapply(data[[var]], list(data[[cluster]]), sd,na.rm=T),digits = 2)
    c0<-c0[!is.na(c0)]


    if (any(c0>.001))
        warning("Variable ",var," does not seem to be a purely between levels variable")

    s0<-tapply(xdata[[y]], list(xdata[[cluster]]), mean,na.rm=T)
    mm<-sd(s0)
    ldata<-aggregate(data[[var]],list(data[[cluster]]),mean,na.rm=T)
    names(ldata)<-c(cluster,zname)
    ldata[[zname]]<-as.numeric(scale(ldata[zname]))*mm
    message("scaling ",var," with sd(y)=",mm)
#    ldata[[zname]]<-ldata[[zname]]-mean(ldata[[zname]],na.rm=T)
    if (zname %in% names(data)) data[[zname]]<-NULL
    xdata<-merge(data,ldata,by=cluster)

  } else {
     data$id....<-1:dim(data)[1]
     zname<-paste0("sw_",var)
     vdata<-c()
     cls<-unique(data[[cluster]])
     tot<-length(cls)
     i<-0
     for (cl in cls) {
      i<-i+1
      .vars<-c("id....",as.character(var),as.character(y))
      ldata<-data[data[[cluster]]==cl,.vars]
      s<-sd(ldata[[y]],na.rm=T)
      ldata[[zname]]<-as.numeric(scale(ldata[[var]]))*s
      ldata[[var]]<-NULL
      ldata[[y]]<-NULL
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


#' @export

levelstandardize<-function(data,var,cluster,level="within", overwrite=FALSE) {

  if (is.factor(data[[var]])) {

    xdata<-data
    if (level=="between") {
      s0<-tapply(xdata[[y]], list(xdata[[cluster]]), mean, na.rm=T)
      mm<-sd(s0)
    } else {
      s0<-tapply(xdata[[y]], list(xdata[[cluster]]), sd, na.rm=T)
      mm<-mean(s0,na.rm=T)
    }
    codes<-contrasts(xdata[[var]])
    codes<-mm*codes/apply(codes,2,sd)
    contrasts(xdata[[var]])<-codes


  } else {

    zname<-var
    if (level=="between") {

      if (!overwrite) zname<-paste0("zb_",var)
      xdata<-data
      c0<-round(tapply(data[[var]], list(data[[cluster]]), sd, na.rm=T),digits = 2)

      if (any(c0>.001))
        warning("Variable ",var," does not seem to be a purely between levels variable")

      ldata<-aggregate(data[[var]],list(data[[cluster]]),mean,na.rm=T)
      names(ldata)<-c(cluster,zname)
      ldata[[zname]]<-as.numeric(scale(ldata[zname]))
      if (zname %in% names(data)) data[[zname]]<-NULL
      xdata<-merge(data,ldata,by=cluster)

    } else {
      data$id....<-1:dim(data)[1]
      zname<-paste0("zw_",var)
      vdata<-c()
      cls<-unique(data[[cluster]])
      tot<-length(cls)
      i<-0
      for (cl in cls) {
        i<-i+1
        ldata<-data[data[[cluster]]==cl,c("id....",var)]
        ldata[[zname]]<-as.numeric(scale(ldata[[var]]))
        ldata[[var]]<-NULL
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


.levelcenter<-function(data,var,cluster,level="within",overwrite=FALSE) {


   zname<-paste0("centered_",var)
   if (overwrite) zname<-var
   m0<-tapply(data[[var]], list(data[[cluster]]), mean,na.rm=TRUE)
   mm<-as.data.frame(cbind(names(m0),m0))
   names(mm)<-c(cluster,".mm.")
   zdata<-merge(data,mm,by=cluster)
   zdata[[".mm."]]<-as.numeric(zdata[[".mm."]])
   if (level=="within") {
      zdata[[zname]]<-zdata[[var]]-zdata[[".mm."]]
    } else {
      zname<-paste0("mean_",var)
      if (overwrite) zname<-var
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


### This produced standardized y with given betas and X
### X must be standardized
### the code is adapted from

### @MISC {508107,
### TITLE = {Simulate regression with specified standardized coefficients},
### AUTHOR = {whuber (https://stats.stackexchange.com/users/919/whuber)},
### HOWPUBLISHED = {Cross Validated},
### NOTE = {URL:https://stats.stackexchange.com/q/508107 (version: 2021-02-05)},
### EPRINT = {https://stats.stackexchange.com/q/508107},
### URL = {https://stats.stackexchange.com/q/508107}
### }

.yield_y<-function(X,beta) {

  n<-dim(X)[1]
  col<-dim(X)[2]
  .X<-as.matrix(X,ncol=col,nrow=n)/(sqrt(n-1))
  dual <- function(aX, tol=1e-16) {
    n <- dim(aX)[1]
    p <- dim(aX)[2]
    qr.Q(qr(cbind(aX, rep(1, n), diag(1, n, n)), tol=tol))[, -seq_len(p+1)]
  }
  y.0 <- .X %*% beta
  y2 <- sum(y.0^2)
  if (y2 > 1) stop("variances and betas not compatible")
  eps.norm <- sqrt(1 - y2) # This will be the length of any error term
  Phi <- with(svd(.X), dual(u[, which(d > 1e-8), drop=FALSE]))
  ##############################################
  # This is the entire algorithm.              #
  # It uses precomputed y.0, Phi, and eps.norm.#
  ##############################################
  eps <- Phi %*% rnorm(dim(Phi)[2])            # Unnormalized errors
  y <- y.0 + eps * eps.norm / sqrt(sum(eps^2)) # Normalize and add
  #
  # Check the generated `y` meets all constraints.
  # maybe we do something with it later on
  #  ok=(1 - sqrt(sum(y^2)) > tol*1e2)
  #
  y*sqrt(n-1)

}


.make_factors<-function(vars,N) {

  .names<-unlist(rlist::list.select(vars,name))
  .values<-list()
  for (var in vars) {
    ladd(.values)<-1:var$k
  }
  .data<-as.data.frame(.expandto(.values,N))
  names(.data)<-.names
  .data
}
.make_dummies<-function(vars,data) {

  .names<-unlist(rlist::list.select(vars,name))
  for (var in vars) {
     data[[var$name]]<-factor(data[[var$name]])
     levs<-nlevels(data[[var$name]])
     cont<-contr.sum(levs)
     colnames(cont)<-1:ncol(cont)
     contrasts(data[[var$name]])<-cont
  }
  .formula<-as.formula(paste0("~",paste0(.names,collapse = "+")))
   data<-as.data.frame(model.matrix(.formula,data))
#  data[[1]]<-NULL
#  for (name in names(data)) data[[name]]<-scale(data[[name]])
  data[[1]]<-NULL

  data
}


.make_clusterdata<-function(cluster,info,n,ysd) {

  bnames<-NULL
  fnames<-NULL
  ysd<-1
  cluster<-rlist::list.find(info$clusters, name==cluster,n=Inf)[[1]]
  .dep<-paste0(cluster$name,"_dep_")
  ### X will be the model matrix
  X<-NULL
  ### finaldata will be the actual data (factors not dummies)
  finaldata<-data.frame(y=rep(0,n))
  names(finaldata)<-.dep
  ## we select the continuous between variables
  vars<-rlist::list.find(info$vars, grep(cluster$name,varying)>0,n=Inf)
  covs<-rlist::list.find(vars, type=="numeric",n=Inf)
  factors<-rlist::list.find(vars, type=="factor",n=Inf)

  if (is.something(covs)) {
    bnames<-unlist(rlist::list.select(covs,name))
    bn<-length(bnames)
    if (hasName(cluster,"sigma") && n>3) {
      sigmab<-cluster$sigma
    } else {
      sigmab<-matrix(0,ncol=bn,nrow=bn)
      diag(sigmab)<-1
    }
    mus<-rep(0,bn)
    X<-as.data.frame(MASS::mvrnorm(n,Sigma = sigmab,mu=mus,empirical = info$empirical))
    names(X)<-bnames
    finaldata<-apply(X,2,scale)
  }
  ## now we check for factors
  if (is.something(factors)) {
    fnames<-unlist(rlist::list.select(factors,name))
    data<-.make_factors(factors,n)
    dummies<-.make_dummies(factors,data)
    fnames<-colnames(dummies)
    if (is.null(X)) X<-dummies else X<-cbind(X,dummies)
    data<-cbind(data,dummies)
    if (is.null(finaldata)) finaldata<-data else finaldata<-cbind(finaldata,data)
  }
  finaldata<-as.data.frame(finaldata)
  vnames<-c(bnames,fnames)
  finaldata
}







