### This is inspired by  simr::powerSim, but extracts all p-values rather
### than one in each simulations
nsim=10
alpha=.05
effects=c("fixed", "random")

.powerMix<-function (model, test = "t", nsim=10, seed=NULL,alpha=.05, effects=c("fixed","random")) {

  pb <- progress::progress_bar$new(
    format = "  simulating [:bar] :percent time: :elapsed",
    total = nsim, clear = FALSE, width= 80)

  if (is.null(seed))
      set.seed(seed)
  .warnings<-0
  .errors<-0

  testfixed<-("fixed" %in% effects)
  testrandom<-("random" %in% effects)

  one<-function(id) {
    pb$tick()
    res<-try_hard({
      values<-list()
      y       <-  do.call(doSim, list(model))
      amodel  <-  do.call(doFit, c(list(y, model)))
      amodel  <-  lmerTest::as_lmerModLmerTest(amodel)
      if (testfixed) {
          fixtest <-  summary(amodel)$coefficients
          values<-fixtest[,5]
      }
      if (testrandom) {

          retest  <-  lmerTest::ranova(amodel)
          retest  <-  retest[-1,]
          rval    <-  t(subset(retest,select=6))
          names(rval)<-colnames(rval)
          values  <-  c(values,rval)

      }
      values
    })
    if (!isFALSE(res$error)) .errors<<-.errors+1
    if (!isFALSE(res$warning)) .warnings<<-.warnings+1
    res$obj
  }

  numCores <- parallel::detectCores()
  reslist<-parallel::mclapply(1:10,one, mc.cores = numCores)
  results<-do.call(rbind,reslist)
  power<-apply(apply(results,2,function(x) x < alpha),2,sum)

  data<-simr::getData(model)
  N<-dim(data)[1]
  clusters<-lme4::getME(model,"l_i")
  clusters<-sapply(names(clusters),function(cluster) length(unique(data[[cluster]])))
  power<-list(nsim=nsim,
              power=power,
              clusters=clusters,
              N=N,
              warnings=.warnings,
              errors=.errors)
  class(power)<-"pamlpowermix"
    #conf<-binom::binom.confint(power,n = nsim)
  power
}


print.pamlpowermix<-function(obj) {

  if (obj$errors>0)
    warning("There were ",.errors," replications with errors")
  if (obj$warnings>0)
    message("There were ",.warnings," replications with warnings")
   df<-as.data.frame(obj)
   cat("Power parameters based on",obj$nsim,"simulations\n\n")
   print(df)
   cat("\nCluster size:",paste(names(obj$clusters),obj$clusters),"\n")
   cat("\nNumber of units (sample size):",obj$N,"\n")

}

as.data.frame.pamlpowermix<-function(obj) {

  df<-data.frame(power=c(obj$fixed,obj$random))
  rownames(df)<-c(names(obj$fixed),names(obj$random))
  suppressWarnings(conf<-binom::binom.confint(df$power,n = obj$nsim,methods = "exact"))
  df$ci_lower<-round(conf$lower,5)
  df$ci_upper<-round(conf$upper,5)
  df$power<-df$power/obj$nsim
  df$type<-c(rep("fixed",length(obj$fixed)),rep("random",length(obj$random)))
  df

}


.evaluatePower<-function(pamlobj,power) {


  df<-as.data.frame(pamlobj)
  req<-matrix(0,ncol=length(levels),nrow=length(df$power))
  levels<-c(units=obj$N,obj$clusters)
  colnames(req)<-names(levels)
  df<-cbind(df,req)
  oks<-which((obj$fixed/obj$nsim)> power)
  for (ok in oks)
    for (lev in names(levels))
      if (df[ok,lev]==0) df[ok,lev]<-levels[lev]
  df

}
