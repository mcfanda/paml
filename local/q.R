.make_y<-function(object, seed = NULL) {



  if(!is.null(seed)) set.seed(seed)
  if(!exists(".Random.seed", envir = .GlobalEnv))
    runif(1) # initialize the RNG if necessary
  RNGstate <- .Random.seed

  sigma <- sigma(object)
  etapred <- predict(object, newdata=NULL, re.form=NA, type="link")
  n <- length(etapred)


  f<-formula(object)
  compReForm<-reformulate(paste0("(", vapply(findbars(f), deparse1, ""), ")"), env = environment(f))

  ## (1) random effect(s)
    newRE <- mkNewReTrms(object, NULL, compReForm)
    U <- t(newRE$Lambdat %*% newRE$Zt) # == Z Lambda
    u <- as.numeric(scale(rnorm(ncol(U))))
    ## UNSCALED random-effects contribution:
    reval<-as(U %*% matrix(u, ncol = 1), "matrix")

    ## result will be matrix  n x nsim :
    err<-rnorm(n)
    err<-as.numeric(scale(resid(lm(err~reval+etapred))))
    sigma<-sqrt(1-var(etapred))
    val<-etapred + reval + sigma * err
    return(val)
}
