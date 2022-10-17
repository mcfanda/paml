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

get_sample<-function(info) {

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
    clusterdata<-list()

    ## gather info for variables
    vars<-info$vars
    dep<-rlist::list.find(vars, varying=="dependent",n=1)[[1]]
    covs<-rlist::list.find(vars, type=="numeric",n=Inf)
    factors<-rlist::list.find(vars, type=="nominal",n=Inf)

    ## create cluster-level data
    for (cluster in clusters) {

      bnames<-NULL
      fnames<-NULL
      ysd<-1
      if (is.something(cluster$ysd)) ysd<-cluster$ysd
      .dep<-paste0(cluster$name,"_dep_")
      ## first, we prepare a frame with the cluster values
      cv<-data.frame(unique(clustervalues[,cluster$name]))
      N<-dim(cv)[1]
      names(cv)<-cluster$name
      ### here we need an empty dataset to be filled later
      X<-NULL

      ## we select the continuous between variables
      bvars<-rlist::list.find(covs, grep(cluster$name,varying)>0,n=Inf)
      if (is.something(bvars)) {
        bnames<-unlist(rlist::list.select(bvars,name))
        bn<-length(bnames)
        if (hasName(cluster,"sigma")) {
          sigmab<-cluster$sigma
        } else {
          sigmab<-matrix(0,ncol=bn,nrow=bn)
          diag(sigmab)<-1
        }
        mus<-rep(0,bn)
        .ccovs<-MASS::mvrnorm(N,Sigma = sigmab,mu=mus,empirical = T)
        X<-as.data.frame(.ccovs)
        names(X)<-bnames
      }
      ## now we check for factors
      fvars<-rlist::list.find(factors, grep(cluster$name,varying)>0,n=Inf)
      if (is.something(fvars)) {
        data<-.make_factors(fvars,n)
        fnames<-colnames(data)
        cv<-cbind(cv,data)
      }

      vnames<-c(bnames,fnames)
      if (is.something(vnames)) {
      ccovs[[.dep]]<-.yield_y(.ccovs,as.numeric(info$fixed[as.character(names)]))
      ccovs[[.dep]]<-ccovs[[.dep]]*ysd
      cv<-cbind(cv,ccovs)
      }
      clusterdata[[cluster$name]]<-cv
    }

    # now the within variables
    within=FALSE
    NC<-dim(clustervalues)[1]
    ## continuous variables

    cvars <- rlist::list.find(covs, varying=="within",n=Inf)
    if (is.something(cvars)) {
        within=TRUE
        wnames<-unlist(rlist::list.select(cvars,name))
        wn<-length(wnames)
        sigmaw<-diag(1,ncol=wn,nrow=wn)
        mus<-rep(0,wn)
        ldata<-list()
        .dep<-"within_dep_"
        ysd<-dep$ysd
        if (is.null(ysd))
             ysd<-1
        if (length(ysd)==1)
             ysd<-rep(ysd,NC)

        n<-dep$n
        if (length(n)==1)
           n<-rep(dep$n,NC)

    for (i in 1:NC) {
     .cv<-MASS::mvrnorm(n[i],Sigma = sigmaw ,mu=mus,empirical = T)
      cv<-as.data.frame(.cv)
      names(cv)<-wnames
      cv[[.dep]]<-.yield_y(.cv,as.numeric(info$fixed[as.character(wnames)]))
      cv[[.dep]]<-cv[[.dep]]*ysd[i]
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

    ### compose the dependent variable
     depnames<-paste0(c("within",clusternames),"_dep_")
     depnames<-depnames[depnames %in% names(finaldata)]
     finaldata[[dep$name]]<-apply(finaldata[,depnames],1,sum)

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
#' @export

make_sample<-function(vars,clusters,structure,formula=NULL,fixed=NULL,random=NULL) {

  if (inherits(formula,"character")) {
       info<-get_formula_info(formula)
       info$vars<-vars
       info$clusters<-clusters
       info$varcorr<-lapply(info$random, function(x) diag(x))
  } else {
       .formula<-formula
       .varcorr<-varcorr
       .beta<-beta
  }
  info$structure<-structure

  sample<-get_sample(info)

  return(sample)

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

levelscale<-function(data,y,var,cluster,level="within", overwrite=FALSE) {

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

#' @export

levelcenter<-function(data,var,cluster,level="within",overwrite=FALSE) {


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
  .X<-as.matrix(X)/sqrt(n-1)

  dual <- function(aX, tol=1e-16) {
    n <- dim(aX)[1]
    p <- dim(aX)[2]
    qr.Q(qr(cbind(aX, rep(1, n), diag(1, n, n)), tol=tol))[, -seq_len(p+1)]
  }
  print(.X)
  y.0 <- .X %*% beta
  y2 <- sum(y.0^2)

  if (y2 > 1) stop("No solutions are possible: betas are too large.")
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

  for (var in vars) {
    .data[[var$name]]<-factor(.data[[var$name]])
    contrasts(.data[[var$name]])<-contr.sum(var$k)
  }
  .formula<-as.formula(paste0("~",paste0(.names,collapse = "+")))
  data<-as.data.frame(model.matrix(.formula,.data))
  data[[1]]<-NULL
  data
}



.make_clusterdata<-function(cluster,info) {

  bnames<-NULL
  fnames<-NULL
  ysd<-1
  cluster<-rlist::list.find(info$clusters, name==cluster,n=Inf)[[1]]

  if (is.something(cluster$ysd)) ysd<-cluster$ysd
  .dep<-paste0(cluster$name,"_dep_")
  ## first, we prepare a frame with the cluster values
  N<-cluster$n
  ### here we need an empty dataset to be filled later
  X<-NULL

  ## we select the continuous between variables
  vars<-rlist::list.find(info$vars, grep(cluster$name,varying)>0,n=Inf)
  covs<-rlist::list.find(vars, type=="numeric",n=Inf)
  factors<-rlist::list.find(vars, type=="nominal",n=Inf)

  if (is.something(covs)) {
    bnames<-unlist(rlist::list.select(covs,name))
    bn<-length(bnames)
    if (hasName(cluster,"sigma")) {
      sigmab<-cluster$sigma
    } else {
      sigmab<-matrix(0,ncol=bn,nrow=bn)
      diag(sigmab)<-1
    }
    mus<-rep(0,bn)
    X<-as.data.frame(MASS::mvrnorm(N,Sigma = sigmab,mu=mus,empirical = T))
    names(X)<-bnames
  }
  ## now we check for factors
  if (is.something(factors)) {
    fnames<-unlist(rlist::list.select(factors,name))
    data<-.make_factors(factors,N)
    fnames<-colnames(data)
    X<-cbind(X,data)
  }

  vnames<-c(bnames,fnames)
  test<-!all(vnames %in% names(info$fixed))
  if (test) stop("Not all fixed effect have an effect size. Fixed effects are: ",paste(vnames,collapse = ", "))
  if (is.something(vnames)) {
    print(as.numeric(info$fixed[as.character(vnames)]))
    X[[.dep]]<-.yield_y(X,as.numeric(info$fixed[as.character(vnames)]))
    X[[.dep]]<-X[[.dep]]*ysd
  }
  X
}







