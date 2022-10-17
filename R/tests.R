
.t_sat<- function(model) {
        .summary<-summary(as_lmerModLmerTest(model))
        .coefs<-.summary$coefficients[,1]
        .se<-.summary$coefficients[,2]
        .df<-.summary$coefficients[,3]
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
