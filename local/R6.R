source("R/functions.R")
source("R/constants.R")
source("R/methods.R")
source("R/syntax.R")

library(R6)
Design <- R6Class("Design",
                  public = list(
                    vars = NULL,
                    initialize = function() {
                    },
                    print=function() {
                      cat("\nCluster variables\n")
                      cat(paste(private$.clusters),"\n\n")
                      cat("Variables\n")
                      garb<-lapply(private$.vars,function(x) {
                        print(unlist(x))
                        cat("\n")

                      } )

                    },
                    check=function() {

                        if (length(private$.clusters)==0) stop("Please specify at least one clustering variable with the option: `obj$clusters<-c(name='aname', layer='integer')`")
                        if (length(private$.vars)==0) stop("Please specify at least independent variable with the option: `obj$variables<-c(name='aname',type='atype', layer='integer')`")

                        v1<-rlist::list.find(private$.vars,layer=="1",n=Inf)
                        test<-(length(v1)==0)
                        if (test) stop("There is no variable varying at level 1. Mixed models are useless in this design.")
                        v1f<-rlist::list.find(v1,type=="factor",n=Inf)
                        l1 <-rlist::list.find(private$.clusters,name=="within",n=1)
                        if (is.something(l1)) l1<-l1[[1]]
                        if (!is.something(v1f) && !is.something(l1))
                              stop("Within clusters size cannot be determined. Please explicitly declare it as `$clusters<-c(name='within', layer=1, size='integer').")
                        if (!is.something(v1f) && !hasName(l1,"size"))
                             stop("Within clusters size cannot be determined. Please explicitly declare it with argument `size='integer'`.")

                        if (is.something(l1)) {
                              if (l1$layer!="1") stop("`within` layer should have `layer=1`")
                        }

                    },
                    size=function(...) {

                       private$.data<-private$.size(...)

                    },


                    create_model=function(family=NULL) {

                       if (is.null(private$.data)) {
                         cls<-unlist(rlist::list.select(private$.clusters,hasName(.,"size")))
                         opt<-sapply(private$.clusters,function(x) {
                           if (hasName(x,"size"))
                             return(x$size)
                           else
                             10
                         },simplify = F)
                         do.call(self$size,opt)
                       }
                       data<-private$.data
                        data[[private$.dep$name]]<-rnorm(nrow(data))
#                       data[[private$.dep$name]]<-rbinom(nrow(data),1,.5)

                       suppressMessages({
                           suppressWarnings({
                             model<-lme4::lmer(private$.info$formula,data=data)
                             oldvc<-lme4::VarCorr(model)
                       })})
                       info<-private$.info$random
                       vc<-sapply(names(oldvc),function(x) {
                           o<-oldvc[[x]]
                           m<-matrix(0,ncol=ncol(o),nrow=nrow(o))
                           colnames(m)<-colnames(o)
                           rownames(m)<-rownames(o)
                           diag(m)<-info[[x]]
                           m
                         },simplify = F)

                       .data<-private$.data
                       suppressMessages(suppressWarnings({
                       if (is.null(family)) {
                            private$.root_model <-  simr::makeLmer(
                            private$.info$formula,
                            fixef = private$.info$fixed,
                            VarCorr = vc,
                            sigma = self$sigma,
                            data = .data
                            )
                       } else {
                         private$.root_model <-  simr::makeGlmer(
                           private$.info$formula,
                           family=family,
                           fixef = private$.info$fixed,
                           VarCorr = vc,
                           data = .data
                         )
                       }
                       }))
                       message("Model created")

                      invisible(private$.root_model)

                     },
                   update_model=function(data=NULL,fixed=NULL,varcorr=NULL,sigma=NULL,...) {

                       private$.root_model<-private$.update_model(data=data,fixed=fixed,varcorr=varcorr,sigma=sigma,...)

                     },

                    one_sample=function(...) {
                      opt<-list(...)
                      if (is.null(private$.root_model))
                         stop("Please create a model before drawing a sample.")
                      if (hasName(opt,"data")) {
                          .data<-opt$data
                          opt$data<-NULL
                      } else
                          .data<-private$.data

                      .model<-private$.root_model

                      if (is.something(opt)) {
                        .data<-do.call(private$.size,opt)
                        .model<-private$.update_model(data=.data)
                      }
                      y<-simr::doSim(.model)
                      .data[[private$.dep$name]]<-y
                      return(.data)
                    },
                    simulate=function(...) {

                      opt<-list(...)
                      if (!hasName(opt,"Nsim"))
                        stop("Please specify the number of simulations with the option `Nsim='integer'.")
                      Nsim<-opt$Nsim
                      opt$Nsim<-NULL

                      clusters<-opt[names(opt) %in% private$.clustersname]
                      params<- opt[! names(opt) %in% private$.clustersname]

                      lapply(private$.clusternames, function(x) {
                          if (! x %in% names(clusters)) {
                            opt[[x]]<-cl$size
                            message("Number of levels for cluster variable ",cl$name," is set to the default ",cl$size)
                          }
                      })
                      test<-names(params)[!names(params) %in% SIMPARAMS]
                      if (length(test)>0)
                           stop("Argument(s)", paste(test,collapse = ", ")," not allowed here.")

                      ### prepare experimental factors
                      alist<-list(clusters=clusters)
                      if (hasName(opt,"fixed")) alist[["fixed"]]<-opt[["fixed"]]
                      if (hasName(opt,"random")) alist[["random"]]<-opt[["random"]]
                      args<-do.call(expand.list,alist)
                      message("Total number of simulations: ",attr(args,"nitems"), " (combinations) X ", Nsim," (Nsim) = ",(Nsim*attr(args,"nitems")))

                      model<-private$.root_model

                      ####


                      doFuture::registerDoFuture()
                      .options<-self$options()
                      if (.options$parallel=="multicore")
                        future::plan(future::multicore)
                      else
                        future::plan(future::multisession)

                      msg<-"Parallel computation started"
                      if (Nsim<11 | .options$parallel==FALSE) {
                        future::plan(future::sequential)
                        msg<-"Sequential computation started"
                      }

                      progressr::handlers(global = TRUE)
                      progressr::handlers("progress")


                      two<-function(model,p) {

                        y<-simr::doSim(model)
                        suppressMessages(suppressWarnings({
                        mod<-simr::doFit(y,model)
                        }))
                        ss<-summary(mod)
                        res<-c()
                        for (i in seq_len(nrow(ss$coefficients))) {
                          tab<-ss$coefficients[i,]
                          opt<-self$options()
                          tab[["power"]]<-as.numeric(tab[[length(tab)]]<opt$alpha)
                          res<-c(res,tab)
                        }
                        rnames<-rownames(ss$coefficients)
                        cnames<-c(colnames(ss$coefficients),"power")
                        .names<-paste(rep(rnames,each=length(cnames)),cnames,sep="_")
                        names(res)<-.names
                        res[["conv"]]<-as.numeric(is.null(mod@optinfo$conv$lme4$messages))
                        p()
                        return(res)
                      }
                      one<-function(i) {
                        params<-lists.atIndex(i,args)
                        clusters<-params[["clusters"]]
                        params[["clusters"]]<-NULL
                        message("Simulating  ",Nsim," times for ",paste(names(clusters),clusters,sep=" = ", collapse = ", "), paste(names(params),params,sep="=",collapse = ", "))

                        data<-do.call(private$.size,as.list(clusters))
                        params[["data"]]<-data
                        .model<-do.call(private$.update_model,params)
                        #                        res<- lapply(1:Nsim, function(j) two(model))
                        p <- progressr::progressor(along = 1:Nsim)
                        res  <- foreach::foreach(1:Nsim) %dorng% two(.model,p)
                        res<-do.call(rbind,res)
                        res<-apply(res,2,mean)
                        res<-c(Nsim=Nsim,params,res)
                        res
                      }
                      time<-Sys.time()
                      results<-lapply(1:attr(args,"nitems"), function(i) one(i))
                      results<-do.call(rbind,results)
                      message("Elapsed time: ",Sys.time()-time," secs")
                      results

                    },

                    options=function(...) {
                      a<-list(...)
                      if (length(a)>0) {
                        private$.options[[names(a)]]<-a[[1]]
                      } else
                         private$.options

                    }
                  ),
                  active= list(
                    N=function(val) {
                      if (missing(val))
                        private$.N
                      else
                        private$.N<-val

                    },
                    clusters = function(cluster) {
                      if (missing(cluster))
                           private$.clusters
                      else {

                          if (!("name" %in% names(cluster)))
                              stop("Please define a cluster as named vector containing the `name=` option")
                        if (!("layer" %in% names(cluster)))
                          stop("Please define a cluster layer (1,2,etc.) usining the `layer=` option")

                          if (!inherits(cluster,"list")) cluster<-as.list(cluster)

                          private$.clusters[[cluster[["name"]]]]<-cluster
                          private$.clustersname[[cluster[["name"]]]]<-cluster[["name"]]

                      }
                    },
                    variables = function(obj) {
                      if (missing(obj))
                        private$.vars
                      else {
                        if (!is.list(obj)) obj<-as.list(obj)
                        if ( (!"name" %in% names(obj)) | (!"type" %in% names(obj)))
                           stop("Please specify a variable as c(name='varname',type='vartype',layer='integer') where 'vartype' can be `numeric`  or `factor`, and varying is a cluster variable across which the variable levels varies. Use `layer=X` for variables varying within level X+1")

                        if ((!"layer" %in% names(obj)) & (!"dependent" %in% names(obj)))
                           warning("No layer specified for variable ",obj$name," use `layer=X` to specify a variable varying within layer X+1 or `dependent=TRUE` for the dependent variable")

                        if (obj$type=="factor" & (!"levels" %in% names(obj)))
                            stop("Factors should have `levels=K` defined, where K is the number of levels")
                        if ((obj$type=="factor") & num(obj$levels)<2)
                          stop("Factor ",obj$name," should have at leat 2 levels")

                        if (obj$type=="numeric" & (!"mean" %in% names(obj))) {
                          obj$mean<-0
                          message("Variable ",obj$name," mean has been set to 0. Use `mean=xxx` to set a different value")
                        }
                        if (obj$type=="numeric" & (!"sd" %in% names(obj))) {
                          obj$sd<-1
                          message("Variable ",obj$name," sd has been set to 1. Use `sd=xxx` to set a different value")
                        }
                        if (("dependent" %in% names(obj)) && (obj["dependent"]==TRUE)) {
                          private$.dep<-obj
                          return()
                        }

                        ladd(private$.vars, key=obj$name)<-obj
                      }
                    },
                    info=function() {
                      private$.info
                    },
                    formula=function(formula) {
                      if (missing(formula))
                          return(private$.formula)
                      else {
                        private$.formula<-formula
                        private$.info<-get_formula_info(formula)
                        private$.checkformula()
                      }
                    },
                    data=function() {
                          return(private$.data)
                      },
                    sigma=function(s) {
                      if (missing(s))
                        return(private$.sigma)
                      else
                        private$.sigma<-s
                    }

                  ), ## end of active
                  private = list(
                    .type="lmer",
                    .N=NULL,
                    .clusters=list(),
                    .clustersname=list(),
                    .vars=list(),
                    .names=NULL,
                    .info=NULL,
                    .formula=NULL,
                    .data=NULL,
                    .dep=list(name="y"),
                    .sigma=1,
                    .root_model=NULL,
                    .options=list(messages=TRUE,
                                  parallel="multisession",
                                  alpha=.05),
                    .checkformula=function() {

                      .var_checks<-function(x) {

                        x<-gsub(" ","",x)

                        if (x=="(Intercept)")
                          return()

                        if (x %in% vars) {
                          if (x %in% facs) {
                            f<-rlist::list.find(private$.vars,type == "factor", n=1)[[1]]
                            ncontr<-as.numeric(as.character(f$levels))-1
                            if (ncontr==1)
                               names(private$.info$fixed)[names(private$.info$fixed)==x] <-paste0(x,".1")
                            else
                              stop("Variable ",x," has ",ncontr+1," levels and requires ",ncontr," contrast variables. Please define its coefficients as ",paste0("`[",paste(1:ncontr,collapse = ","),"]*",x,"`"))
                          }
                          return()
                        }
                        root<-gsub("\\.[0-9]*$","",x)
                        var<-private$.vars[[root]]

                        if (! root %in% vars)
                          stop("Variable ",x," is not in the design.")
                        if (root %in% covs)
                          stop("`var.digit` names are reserved for categorical contrasts but variable ",test," is numeric")


                        test<-grep(paste0("^",root,"\\.[0-9]*$"),names(private$.info$fixed))
                        if (length(test)!=(as.numeric(var$levels)-1))
                          stop("Variable ",
                               var$name,
                               " is defined with ",
                               var$levels,
                               " levels but ",
                               length(test),
                               " contrast variables coefficients are provided. ",
                               "Variable ",var$name," requires ",length(test)-1,
                               " contrast coefficient(s).")
                      }

                      self$check()
                      facs<-names(rlist::list.find(private$.vars,type == "factor", n=Inf))
                      covs<-names(rlist::list.find(private$.vars,type == "numeric", n=Inf))
                      vars<-c(facs,covs)
                      garb<-lapply(names(private$.info$fixed), .var_checks)
                      garb<-lapply(private$.info$clusters,function(x){
                        if (x %in% names(private$.clustersname))
                          return()
                        stop("Cluster Variable ",x," is not in the design")
                      })

                      lapply(private$.info$random, function(x) {
                        lapply(names(x), .var_checks)
                      })

                    },
                    .make_layer=function(layer,match=NULL) {

                      wfac<-list()
                      wcon<-list()
                      l<-as.character(layer)
                      clusters<-private$.clusters[unlist(rlist::list.map(private$.clusters,layer==as.character(l)))]
                      if (length(clusters)==0) {
                        cdata<-NULL
                        Nc<-0
                      } else {
                            clusters<-rlist::list.find(private$.clusters,layer==as.character(l),n = Inf)
                            clist<-lapply(clusters, function(x) {
                                 if (!is.null(match))
                                   paste(match,1:x$size,sep="_")
                                 else
                                    1:x$size
                            })
                            cdata<-as.data.frame(expand.grid(clist))
                            Nc<-nrow(cdata)
                      }
                      vars<-private$.vars[unlist(rlist::list.map(private$.vars,layer==l))]
                      for (x in vars) {
                          if (x$type=="factor")
                            ladd(wfac)<-x
                          else
                            ladd(wcon)<-x
                      }
                      if (is.something(wfac)) {
                          varsnames<-unlist(rlist::list.select(wfac,name))
                          varslevels<-as.numeric(unlist(rlist::list.select(wfac,levels)))
                          wdata<-as.data.frame(expand.grid(lapply(varslevels,function(x) 1:x)))
                          names(wdata)<-varsnames
                          Nw<-nrow(wdata)

                          if (Nc==0) Nc<-Nw
                          dif<-Nc-Nw
                          if (dif<0) stop("The levels in layer ",layer,"
                                          are less than the required levels given by factor(s): ",paste(varsnames,collapse = ", "),
                                          ". Layer levels:",Nc," Factors levels:",Nw)
                          if (dif>0) {
                            t<-ceiling(Nc/Nw)
                            wdata<-wdata[rep(1:nrow(wdata),t),]
                            if (!inherits(wdata,"data.frame")) {
                                wdata<-as.data.frame(wdata)
                            }
                          }
                          .names<-c(names(cdata),varsnames)
                          cdata<-as.data.frame(cbind(cdata,wdata[1:Nc,]))
                          names(cdata)<-.names
                      }
                      varsnames<-unlist(rlist::list.select(wcon,name))
                      if (length(varsnames)>0) {
                            varsmean<-as.numeric(unlist(rlist::list.select(wcon,mean)))
                            varssd<-as.numeric(unlist(rlist::list.select(wcon,sd)))
                            k<-length(varsmean)
                            sigma<-matrix(0,nrow = k,ncol=k)
                            diag(sigma)<-varssd^2
                            ndata<-as.data.frame(MASS::mvrnorm(Nc,mu = varsmean,Sigma = sigma,empirical = T))
                            names(ndata)<-varsnames
                            cdata<-cbind(cdata,ndata)
                      }
                      cdata[[paste0("layer_",layer)]]<-1:nrow(cdata)
                      if (is.something(match))
                        cdata[[paste0("layer._.",layer)]]<-paste(match,1:nrow(cdata),sep="_")

                      return(cdata)
                    },
                    .make_contrast=function(nlevels) {

                      dummy <- stats::contr.treatment(nlevels)
                      coding <- matrix(rep(1/nlevels, prod(dim(dummy))), ncol=nlevels-1)
                      contrast <- (dummy - coding)
                      colnames(contrast)<-paste0(".",1:(nlevels-1))
                      contrast
                    },
                    .size=function(...) {
                      alist<-list(...)
                      sizes<-list()
                      cnames<-private$.clustersname

                      test<-setdiff(names(private$.clusters),names(alist))
                      if (length(test)>0) {
                        cls<-private$.clusters[test]
                        for (x in cls) {
                          if (hasName(x,"size"))
                            alist[[x$name]]<-x$size
                        }
                      }
                      test<-setdiff(names(private$.clusters),names(alist))
                      if (length(test)>0)
                        stop("Cluster variable ",paste(test,collapse = ", ")," requires a size (# of groups)")

                      if (hasName(alist,"within") & !hasName(private$.clusters,"within"))
                        private$.clusters[["within"]]<-list(name="within",layer=1,size=as.numeric(alist$within))

                      self$check()
                      lapply(names(alist),function(x) {
                        if (!x %in% names(private$.clusters))
                          stop("Cluster variable ",x," not defined.")
                        private$.clusters[[x]]$size<-alist[[x]]
                      })
                      layers<-as.numeric(unlist(rlist::list.select(private$.clusters,layer)))
                      layers<-unique(sort(c(layers,1)))
                      data<-private$.make_layer(max(layers))
                      for (j in rev(layers)[-1]) {
                        rdata<-list()
                        for (i in 1:nrow(data)) {
                          suppressWarnings({
                            ladd(rdata)<-cbind(data[i,],private$.make_layer(j,i))
                          })
                        }
                        data<-do.call(rbind,rdata)
                      }
                      facs<-rlist::list.find(private$.vars,type == "factor", n=Inf)
                      for (v in facs) {
                        data[[v$name]]<-factor(data[[v$name]])
                        contrasts(data[[v$name]])<-private$.make_contrast(as.numeric(v$levels))
                      }
                      clusters<-names(private$.clusters)
                      if ("within" %in% clusters)  clusters<-clusters[which(clusters!="within")]
                      for (v in clusters) data[[v]]<-factor(data[[v]])
                      return(data)
                    },
                    .update_model=function(data=NULL,fixed=NULL,random=NULL,sigma=NULL,...) {

                      if (is.null(private$.root_model))
                        stop("No model to update")

                      if (is.null(data)) .data<-private$.data else .data<-data

                      .fixed<-private$.info$fixed
                      if (!is.null(fixed)) {
                        .fixed[names(fixed)]<-fixed
                      }
                      .random<-private$.info$random
                      if (!is.null(random)) {
                        for (param in names(random)) {
                          for (acluster in names(.random ) ) {
                                  cl<-grep(paste0("^",acluster,"\\."),param)
                                  clean<-gsub(paste0("^",acluster,"."),"",param)
                                  if (length(cl)>0) {
                                     .random[[acluster]][[clean]]<-random[param]
                                  }
                        }
                        }
                      }
                      .sigma<-private$.sigma

                      oldvc<-lme4::VarCorr(private$.root_model)
                      vc<-sapply(names(oldvc),function(x) {
                        o<-oldvc[[x]]
                        m<-matrix(0,ncol=ncol(o),nrow=nrow(o))
                        colnames(m)<-colnames(o)
                        rownames(m)<-rownames(o)
                        diag(m)<-.random[[x]]
                        m
                      },simplify = F)
                      .family<-family(private$.root_model)


                      suppressMessages(suppressWarnings({
                        if (.family$family=="gaussian") {
                          model <-  simr::makeLmer(
                            private$.info$formula,
                            fixef = .fixed,
                            VarCorr = vc,
                            sigma = .sigma,
                            data = .data
                          )
                        } else {
                          model <-  simr::makeGlmer(
                            private$.info$formula,
                            family=.family,
                            fixef = .fixed,
                            VarCorr = vc,
                            data = .data
                          )
                        }
                      }))

                      model

                    }

                  ) # end of private
)

d<-Design$new()
d$clusters<-c(name="endo_id", layer=2)
d$clusters<-c(name="video_id", layer=2)
#d$clusters<-c(name="within", layer=1, size=50)

#d$variables<-c(name="afac1",type="factor",layer=1,levels=2)
d$variables<-c(name="afac2",type="factor",layer=1,levels=2)
#d$variables<-c(name="x1",type="numeric",layer=1)
#d$variables<-c(name="x3",type="numeric",layer=3,mean=100)
#d$variables<-c(name="y",type="numeric", dependent=T)
#data<-d$size(endo_id=10,video_id=10, within=20)
#form<-"y~[1]*1+[2,3]*afac2+[4]*x1+(0+[5]*x1|endo_id)+([6]*1|endo_id)+([7]*1+[8]*x1|video_id)"

form<-"y~1+ [.25]*afac2+(1+afac2|endo_id)+([.3]*1+[.25]*afac2|video_id)"
#form<-"y~[1]*1+[0,0]*afac2+[0]*x1+([.2]*1+[.2]*x1|endo_id)"
#form<-"y~[1]*1+[.3]*x1+([.3]*1+[.3]*x1|endo_id)"
d$formula<-form
d$info
q<-d$create_model(family=binomial())
#simulate(q,newdata=as,allow.new.levels=T)

d$options(parallel="multisession")

p1<-.7
p2<-p1+.05
odd<-(p2/(1-p2))/(p1/(1-p1))
cat("Odd: ",odd," slope variance:",log(odd),"\n")

#ss<-d$simulate(Nsim=5,endo_id=c(20,25),video_id=10,fixed=list("(Intercept)"=.3,afac2.1=c(.4,.5)),random=list(afac2.1=2))
#ss<-d$simulate(Nsim=50,endo_id=20,video_id=405)
#ss
