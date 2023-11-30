#' @title  Design R6 Class
#' @description
#' The R6 object that defines the sample, the design matrix and runs the simulations for power calculations
#'
#' @export

RDesign <- R6::R6Class("PDesign",
                  public = list(
                    vars = NULL,
                    #' @description
                    #' Initialize the Design object. No arguments are required.
                    #' @md
                    initialize = function() {
                    },
                    #' @description
                    #' Print some information regarding the design being set up
                    print=function() {
                      cat("\nCluster variables\n")
                      cat(paste(private$.clusters),"\n\n")
                      cat("Variables\n")
                      garb<-lapply(private$.vars,function(x) {
                        print(unlist(x))
                        cat("\n")
                      } )
                      if (!is.null(private$.info)) {
                        print(private$.info)
                      } else
                        message("No formula has been defined")
                    },
                    #' @description
                    #' Checks some crucial parameters of the `Design` object being defined.
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
                    #' @description
                    #' Define new sample size parameters. It accepts clusters name as arguments in the form
                    #' `obj$size(cluster1='integer', cluster2='integer', ...)`. Not necessary to use when `obj$simulate()` is used.
                    size=function(...) {
                       private$.data<-private$.size(...)
                    },
                    #' @description
                    #' Create a root model as a basis of the power simulations.
                    #' @param family a \code{\link[stats]{stats::family}}  that defines the type of model to be estimated. Any family accepted by
                    #' \code{\link[lme4]{lme4::glmer}} is accepted. If not defined (default) a linear model is assumed. If random coefficients
                    #' are defined (such as `(x|cluster)` in the formula), a mixed model is assumed, eastimated with \code{\link[lme4]{lme4::lmer}}. If no random effects are defined, a general(ized) linear
                    #' model is estimated with \code{\link[stats]{lme}}
                    create_model=function(family=NULL) {

                      if (is.null(private$.info))
                         stop("No formula specified, please enter a formula with the command `obj$formula<-'astring'`")

                       if (is.null(private$.data)) {
                         cls<-unlist(rlist::list.select(private$.clusters,hasName(.,"size")))
                         opt<-sapply(private$.clusters,function(x) {
                           if (hasName(x,"size"))
                             return(x$size)
                           else {
                             message("Number of clusters for variable ",x$name, " was not defined. It is set to 10 as provisional size.")
                             10
                           }
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
                      #' @description
                      #' Update a previously defined model with new coefficients or new sample sizes.
                      #' @param data a new dataframe to use to estimate the model
                      #' @param fixed a new set of fixed coefficient defined as `fixed=list(varname1=number,varname2=number)`
                      #' @param random a new set of random coefficient defined as `random=list(acluster=list(varname1=number,varname2=number))`


                   update_model=function(data=NULL,fixed=NULL,random=NULL,sigma=NULL,...) {

                       private$.root_model<-private$.update_model(data=data,fixed=fixed,random=random,sigma=sigma,...)

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
                      #' @description
                      #' Define a simulation experiment. It accepts several arguments, depending to the type of experiment the user wishes to use
                      #' @param Nsim the number of simulations (repetitions) for each cell of the simulation experimental design.
                      #' @param fixed a named list of of expected parameters of the form `(fixed=list('(Intercept)'='number',aname='number',...)`.
                      #' Parameters not defined here are set accordingly to the formula passed with `obj$formula<-astring` method.
                      #' @param random a named list of of expected parameters of the form `(random=list(acluster=list('(Intercept)'='number',aname='number',...))`.
                      #' Parameters not defined here are set accordingly to the formula passed with `obj$formula<-astring` method.
                      #' @param ... any cluster variable defined with `obj$clusters` with a required size, passed as `obj$simulate(aclustername='integer'` where
                      #' `integer` is the number of clusters required. Clusters not sized here are set accordingly to the size passed by the user with
                      #' `obj$cluster(...,size='integer'`) or to the default provisional size=10
                      #'

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
                           stop("Argument(s) ", paste(test,collapse = ", ")," not allowed here.")

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
                    cluster = function(cluster) {
                      if (missing(cluster))
                           private$.clusters
                      else {
                          if (!inherits(cluster,"list")) cluster<-as.list(cluster)

                          if (!("name" %in% names(cluster)))
                              stop("Please define a cluster as named vector containing the `name=` option")

                        if (cluster$name!="within" && !("layer" %in% names(cluster)))
                          stop("Please define a cluster layer (2,3,etc.) usining the `layer=` option")

                        if (cluster$name!="within" && num(cluster$layer)<2)
                           stop("Clusters layer should be larger than 1. Layer 1 is the within cluster layer.")


                          private$.clusters[[cluster[["name"]]]]<-cluster
                          private$.clustersname[[cluster[["name"]]]]<-cluster[["name"]]

                      }
                    },
                    variable = function(obj) {

                      if (missing(obj))
                        private$.vars
                      else {
                        if (!is.list(obj)) obj<-as.list(obj)

                        if ( (!"name" %in% names(obj)) | (!"type" %in% names(obj)))
                           stop("Please specify a variable as c(name='varname',type='vartype',layer='integer') where 'vartype' can be `numeric`  or `factor`, and varying is a cluster variable across which the variable levels varies. Use `layer=X` for variables varying within level X+1")

                        if ((!"layer" %in% names(obj)) & (!"dependent" %in% names(obj)))
                           warning("No layer specified for variable ",obj$name," use `layer=X` to specify a variable varying within layer X+1 or `dependent=TRUE` for the dependent variable")

                        if (obj$type=="factor" && (!"levels" %in% names(obj)))
                            stop("Factors should have `levels=K` defined, where K is the number of levels")

                        if ((obj$type=="factor") && num(obj$levels)<2)
                           stop("Factor ",obj$name," should have at leat 2 levels")

                        if (obj$type=="numeric" && (!"mean" %in% names(obj))) {
                          obj$mean<-0
                          message("Variable ",obj$name," mean has been set to 0. Use `mean='number` to set a different value")
                        }
                        if (obj$type=="numeric" && (!"sd" %in% names(obj))) {
                          obj$sd<-1
                          message("Variable ",obj$name," sd has been set to 1. Use `sd='number` to set a different value")
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
                    data=function(data) {
                      if (missing(data))
                          return(private$.data)
                      else {
                         for (v in names(data)) {

                           if (is.factor(data[[v]])) {
                             cnt<-contrasts(data[[v]])
                             colnames(cnt)<-paste0(".",1:length(colnames(cnt)))
                             contrasts(data[[v]])<-cnt
                           }
                         }

                        for (cluster in names(self$info$random)) {
                          clustervar<-data[[cluster]]
                          n<-length(unique(clustervar))

                          if (!hasName(private$.clusters, cluster)) {
                               private$.clusters[[cluster]]<-list(name=cluster,size=n)
                          } else {
                            if (hasName(private$.clusters[[cluster]], "size")) {
                              if (num(private$.clusters[[cluster]]$size)!=n)
                                  message("Cluster variable ",cluster," size defined by the users (N=",private$.clusters[[cluster]]$size,") but in the data it has size N=",n,". The data size is used.")
                            }
                          }
                        }

                         private$.data<-data

                         private$.external_data=TRUE
                      }
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
                    .external_data=FALSE,
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
                      cat("layer ", l ,"\n")
                      clusters<-private$.clusters[unlist(rlist::list.map(private$.clusters,layer==l))]
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
                      ### here we create data for each layer. It may seem a bit slow
                      ### but as compared with running simulations this procedure
                      ### time does not really impact the computational time of the
                      ### powe analysis. A data.frame is create once for all simulation repetitions
                      ### within one sample size
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
                      ### data are ready, now simply check variables types

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
                        stop("No model to update. Please Specify a model with the command `obj$create_model(...)` ")

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
                      .sigma<-sigma
                      if (is.null(sigma))
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
