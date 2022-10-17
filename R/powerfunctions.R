alpha=.05
.powerMixTest<-function (model,
                         test = "t",
                         df   = "sat",
                         alpha=.05,
                         effects=c("fixed","random")) {
         test<-paste0(".",test,"_",df)
         res<-do.call(eval(parse(text=test)),list(model))

         ncp<-res$coefs/res$se
         tside<-2
         crit<-qt(alpha/tside,res$df,lower.tail = F)
         pwr<-pt(crit, res$df, ncp = ncp, lower.tail = FALSE) +  pt(-crit, res$df, ncp = ncp, lower.tail = TRUE)
         power<-list(coefs=res$coefs,df=res$df,power=pwr)
         data<-model@frame
         N<-dim(data)[1]
         clusters<-lme4::getME(model,"l_i")
         clusters<-sapply(names(clusters),function(cluster) length(unique(data[[cluster]])))
         clustersn<-sapply(names(clusters),function(cluster) round(mean(table(data[[cluster]]))))

         struct<-list(power=power,
                      clusters=clusters,
                      clustersn=clustersn,
                      N=N,
                      dfmethod=df,
                      messages=NULL,
                      warnings=NULL,
                      errors=NULL,
                      time=NULL)
         class(struct)<-"pamlpowermix"
         struct
}



### This is inspired by  simr::powerSim, but extracts all p-values rather
### than one in each simulations

.powerMixSims<-function (model,
                     test = "t",
                     nsim=10,
                     seed=NULL,
                     alpha=.05,
                     effects=c("fixed","random"),
                     parallel=TRUE) {

  time<-Sys.time()

  if (is.null(seed))
      set.seed(seed)


  testfixed<-("fixed" %in% effects)
  testrandom<-("random" %in% effects)

  one<-function(model) {

    .error<-NULL
    .warning<-NULL
    .message<-NULL

      .exp<-function() {
      values<-list()
      y       <-  do.call(simr::doSim, list(model))
      amodel  <-  do.call(simr::doFit, c(list(y, model)))
      amodel  <-  lmerTest::as_lmerModLmerTest(amodel)
      if (testfixed) {
          fixtest <-  summary(amodel)$coefficients
          values<-fixtest[,5]
      }
      if (testrandom) {
          newData <-  amodel@frame
          retest  <-  lmerTest::ranova(amodel)
          retest  <-  retest[-1,]
          rval    <-  t(subset(retest,select=6))
          names(rval)<-colnames(rval)
          values  <-  c(values,rval)
      }
      values
      }
      results<-withCallingHandlers(
        tryCatch(.exp(), error=function(e) {
          .error<<-conditionMessage(e)
          NULL
        }), warning=function(w) {
          .warning<<-conditionMessage(w)
          invokeRestart("muffleWarning")
        }, message = function(m) {
          .message<<-conditionMessage(m)
          invokeRestart("muffleMessage")
        })

      p()
      if (is.null(results)) results<-list()
      if (!is.null(.error))   attr(results,"error")   <- .error
      if (!is.null(.warning)) attr(results,"warning") <- .warning
      if (!is.null(.message)) attr(results,"message") <- .message

      results
  }

  doFuture::registerDoFuture()
  future::plan(future::multisession)
  msg<-"Parallel computation started"
  if (nsim<11 || !parallel) {
     future::plan(future::sequential)
     msg<-"Sequential computation started"
  }
  message(msg)
  progressr::handlers(global = TRUE)
  progressr::handlers("progress")
  if (is.null(attr(model,"newData")))
    attr(model,"newData")<-model@frame

  if (is.something(seed)) set.seed(seed)

  p <- progressr::progressor(along = 1:nsim)
  reslist  <- foreach::foreach(1:nsim) %dorng% one(model)
  .errors<-unlist(lapply(reslist, function(x) attr(x,"error")))
  .warnings<-unlist(lapply(reslist, function(x) attr(x,"warning")))
  .messages<-unlist(lapply(reslist, function(x) attr(x,"message")))
  results<-do.call(rbind,reslist)
  message("Computation ended")
  ## now we compute the power
  pwr<-apply(apply(results,2,function(x) x < alpha),2,sum)
  power=list(coefs=fixef(model),power=pwr)
  ## then we get some info
  data<-simr::getData(model)
  N<-dim(data)[1]
  clusters<-lme4::getME(model,"l_i")
  clusters<-sapply(names(clusters),function(cluster) length(unique(data[[cluster]])))
  clustersn<-sapply(names(clusters),function(cluster) round(mean(table(data[[cluster]]))))

  struct<-list(nsim=nsim,
              power=power,
              clusters=clusters,
              clustersn=clustersn,
              N=N,
              messages=.messages,
              warnings=.warnings,
              errors=.errors,
              time=as.numeric(Sys.time()-time))
  class(struct)<-"pamlpowermix"
  struct
}




.evaluatePower<-function(obj,power) {

  df<-as.data.frame(obj)
  levels<-c(units=obj$N,obj$clusters)
  req<-matrix(0,ncol=length(levels),nrow=length(df$power))
  colnames(req)<-names(levels)
  df<-cbind(df,req)
  oks<-which(df$power> power)
  for (ok in oks)
    for (lev in names(levels))
      if (df[ok,lev]==0) df[ok,lev]<-levels[lev]
  df

}
