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

                        v1<-private$.vars[unlist(rlist::list.map(private$.vars,as.numeric(layer)==1))]
                        test<-(length(v1)==0)
                        if (test) stop("There is no variable varying at level 1. Mixed models are useless in this design.")

                        test1<-!is.something(v1[unlist(rlist::list.map(v1,type=="factor"))])
                        test2<-!is.something(private$.clusters[unlist(rlist::list.map(private$.clusters,name=="within"))])
                        if (test1 && test2) stop("Within clusters size cannot be determined. Please specify it with obj$size(within=xx)")

                        if (is.null(private$.dep))
                               warning("No dependent variable has been defined. The dependent variable distribution will be determined by the type of model being estimated.")



                    },
                    size=function(...) {

                       alist<-list(...)
                       sizes<-list()
                       cnames<-private$.clustersname

                       test<-setdiff(names(private$.clusters),names(alist))
                       if (length(test)>0) stop("Cluster variable ",paste(test,collapse = ", ")," requires a size (# of groups)")

                       if (hasName(alist,"within"))
                         private$.clusters[["within"]]<-list(name="within",layer=1)

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
                       private$.data<-data

                    },
                    create=function(type) {
                      if (missing(type))
                          type<-private$.type

                      if (type=="lmer") {

                       info<-private$.info$random
                       vc<-sapply(info,function(x) {
                           m<-matrix(0,ncol=length(x),nrow=length(x))
                           colnames(m)<-names(x)
                           rownames(m)<-names(x)
                           diag(m)<-x
                           m
                         })
                      sim_model <-  simr::makeLmer(
                        private$.info$formula,
                        fixef = private$.info$fixed,
                        VarCorr = vc,
                        sigma = 1,
                        data = private$.data
                      )
                      simdata$dep<-simr::doSim(sim_model)


                     }

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
                      print(private$.info)
                    },
                    formula=function(formula) {
                      if (missing(formula))
                          return(private.formula)
                      else {
                        private$.formula<-formula
                        private$.info<-get_formula_info(formula)
                        private$.checkformula()
                      }
                    },
                    data=function() {
                          return(private$.data)
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
                    .dep=NULL,
                    .checkformula=function() {

                      .var_checks<-function(x) {

                        if (x=="(Intercept)")
                          return()
                        if (x %in% vars)
                          return()
                        test<-gsub("\\.[0-9]*$","",x)

                        var<-private$.vars[[unique(test)]]

                        if (! unique(test) %in% vars)
                          stop("Variable ",x," is not in the design.")
                        if (unique(test) %in% covs)
                          stop("`var.digit` names are reserved for categorical contrasts but variable ",test," is numeric")

                        if (length(test)!=(as.numeric(var$levels)-1))
                          stop("Variable ",var$name," is defined with ",var$levels," levels but ",length(test)," contrasts variables coefficients are provided")

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
                      print(private$.info)

                      lapply(private$.info$random, function(x) {
                        lapply(names(x), .var_checks)
                      })

                    },
                    .make_layer=function(layer,match=NULL) {

                      wfac<-list()
                      wcon<-list()
                      l<-as.character(layer)
                      clusters<-private$.clusters[unlist(rlist::list.map(private$.clusters,layer==l))]

                      if (length(clusters)==0) {
                        cdata<-NULL
                        Nc<-0
                      } else {
                            clusters<-private$.clusters[unlist(rlist::list.map(private$.clusters,layer==l))]
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
                          lapply(varsnames,function(x) {
                            cdata[[x]]<-factor(cdata[[x]])
                            contrasts(cdata[[x]])<-contr.sum(as.numeric(vars[[x]]$levels))
                          })
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
                    }
                  ) # end of private
)

d<-Design$new()
d$clusters<-c(name="endo_id", layer=2)
d$clusters<-c(name="video_id", layer=2)
#d$variables<-c(name="afac1",type="factor",layer=1,levels=2)
d$variables<-c(name="afac2",type="factor",layer=1,levels=2)
d$variables<-c(name="x1",type="numeric",layer=1,mean=100)
d$variables<-c(name="x3",type="numeric",layer=3,mean=100)
d$variables<-c(name="y",type="numeric", dependent=T)
data<-d$size(endo_id=16,video_id=10)
data
form<-"y~[1]*1+[2,3]*afac2+[4]*x1+(0+[5]*x1|endo_id)+([6]*1|endo_id)+([7]*1|video_id)"
d$formula<-form
nrow(data)
length(table(data$endo_id))
length(table(data$video_id))


