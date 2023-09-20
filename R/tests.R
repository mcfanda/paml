
.t_sat<- function(model) {

  doFuture::registerDoFuture()
  future::plan(future::multisession)

  one<-function() {
    y<-do.call(simr::doSim, list(model))
    amodel  <-  do.call(simr::doFit, c(list(y, mod)))
    .summary<-summary(as_lmerModLmerTest(amodel))
    .coefs<-.summary$coefficients[,1]
    .se<-.summary$coefficients[,2]
    .df<-.summary$coefficients[,3]
    res<-list(coefs = .coefs,
              se    = .se,
              df    = .df)
    return(res)
  }

  results  <- foreach::foreach(x=1:10) %dorng%  one()
  qq<-do.call(cbind,results)
  .coefs<-apply(do.call(rbind,qq[1,]),2,mean)
  .se<-apply(do.call(rbind,qq[2,]),2,mean)
  .df<-apply(do.call(rbind,qq[3,]),2,mean)
  res<-list(coefs = .coefs,
                  se    = .se,
                  df    = .df)
  return(res)
}


.t_kr<- function(model) {

  .summary<-summary(model)
  .coefs<-.summary$coefficients[,1]
  .se<-.summary$coefficients[,2]
   nv<-length(.coefs)
   lmat<-diag(1,nv,nv)
   .df<-apply(lmat,1, function (l) {
     unlist(pbkrtest::get_Lb_ddf(model,l))
   })
   names(.df)<-names(.coefs)
  res<-list(coefs = .coefs,
            se    = .se,
            df    = .df)

  return(res)
}
