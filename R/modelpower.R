#' Compute the power of a linear model
#'
#' Compute power of parameters tests linear models
#'
#' @param  model A model of of class \code{\link[lme4]{lmerMod}}
#' @param  alpha Significance level (Type I error probability)
#' @param ...    additional options that depends on the class of model
#'               passed to the function.
#' @return an object of class \code{pamlpower*},
#'         depending on the class of the model that is evaluated
#' @details Exactly one of the parameters n,
#'  power, alpha must be passed as NULL, and that parameter
#'  is determined from the others.
#'  Notice that alpha has non-NULL default so NULL must be
#'  explicitly passed if you want it computed.
#' @author Marcello Gallucci
#' @export

modelpower<- function(model,alpha,method="test",nsim=100,parallel=TRUE,seed=NULL,...) UseMethod("modelpower")

#' @export
modelpower.default<-function(model,alpha=.05,nsim=100) {

  stop("class",class(model),"not implemented")
}

#' Compute power mixed models
#'
#' Compute the power of a mixed model of the class
#' \code{\link[lme4]{lmerMod}} or \code{\link[lmerTest]{lmerModLmerTest-class}}

#' @param  model of class \code{\link[lme4]{lmerMod}}
#' @param  alpha significance level (Type I error probability)
#' @param  method `test` (default) computes power based on the t-test and
#'                model degrees fo freedom. `sim` computes power based on simulations.
#' @param  nsim   Number of simulations for `method="sim"`.Default is 100 to allow trying
#'                out the function. For accurate results, at least 1000 is recommended.
#' @param  effects a character string or character vector, describing for which
#'                 type of coefficients the power is required.
#'                 It could be \code{fixed} (default),\code{random} or both (\code{c("fixed","random")}).
#'
#' @param  parallel TRUE (default) of FALSE. Whether the simulations are run
#'                  in parallel or in sequence. If the required number of
#'                  simulations (`nsim`) is equal or less than 10, FALSE is
#'                  enforced. Parallel computation is implemented based on \code{\link[future]} package.

#' @param  seed     optional argument to \code{\link[base]{set.seed}}
#' @return an object of class \code{\link[paml]{pamlpowermix}}
#' @details Exactly one of the parameters n,
#'  power, alpha must be passed as NULL, and that parameter
#'  is determined from the others.
#'  Notice that alpha has non-NULL default so NULL must be
#'  explicitly passed if you want it computed.
#' @author Marcello Gallucci
#' @export

modelpower.lmerMod<-function(model,
                             alpha=.05,
                             method="test",
                             nsim=100,
                             test="t",
                             df="sat",
                             effects=c("fixed"),
                             seed=NULL,
                             parallel=TRUE) {

  if (method=="sim")
      results<- .powerMixSims(model,
                            nsim = nsim,
                            alpha = alpha,
                            effects=effects,
                            seed = seed,
                            parallel=parallel)
  else
      results<-.powerMixTest(model=model,
                                   test=test,
                                   df=df,
                                   alpha=alpha,
                                   effects=effects)

  results
}


