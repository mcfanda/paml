
#' Print results of paml power analysis
#'
#' Prints in a nice tabular shape the results of
#' \code{\link[paml]{modelpower}} applied to a mixed model.
#'
#' @param  obj of class \code{\link[paml]{pamlpowermix}}
#' @return an object of class \code{\link[paml]{pamlpowermix}}
#' @author Marcello Gallucci
#' @examples
#'  methods(class="pamlpowermix")
#' @export

print.pamlpowermix<-function(obj,verbose=TRUE) {

  res<-as.data.frame(obj)
  if (is.something(obj$nsim))
        cat("Power parameters based on",obj$nsim,"simulations\n")
  else
        cat("Power parameters based on model standard errors and t-tests \n")

  if ("df" %in% names(obj))
    cat("Degrees of freedom method:",obj$dfmethod," \n")

  if ("target" %in% names(obj))
    cat("Target effect:",obj$target,"\n\n")
  print(res)
  cat("\nClusters levels:",paste(names(obj$clusters),obj$clusters))
  cat("\nWithin clusters avg size:",paste(names(obj$clusters),obj$clustersn))

  cat("\nNumber of units (sample size):",obj$N,"\n")
  if (verbose) {

    if (is.something(obj$time))
       cat("\nElapsed time:",obj$time,"s\n")

    if (length(obj$errors)>0) message(length(obj$errors)," runs returned an error")
    if (length(obj$messages)>0) message(length(obj$messages)," runs returned a message")
    if (length(obj$warnings)>0) message(length(obj$warnings)," runs returned a warning")
  }
  if (!is.something(obj$nsim))
    message("Power estimation based on the t-distribution may be missleading. Please consider `method='sim'` for simulations based power estimates.\n")

}


#' Print results of paml power analysis
#'
#' Prints in a nice tabular shape the results of
#' \code{\link[paml]{modelpower}}) applied to a mixed model.
#'
#' @param  obj of class \code{\link[paml]{pamlpowermix}}
#' @return an object of class \code{\link[base]{data.frame}}
#' @author Marcello Gallucci
#' @examples
#'  methods(class="pamlpowermix")
#' @export


as.data.frame.pamlpowermix<-function(obj) {

  res<-data.frame(obj$power)
  if (is.something(obj$nsim)) {
        suppressWarnings(conf<-binom::binom.confint(res$power,n = obj$nsim,methods = "exact"))
        res$ci_lower<-round(conf$lower,5)
        res$ci_upper<-round(conf$upper,5)
        res$power<-res$power/obj$nsim
  }
  res

}


