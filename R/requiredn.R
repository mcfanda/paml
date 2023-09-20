#' Find required N for linear models
#'
#' Compute required N for several classes of linear models.
#'
#' @param  model of class \code{\link[lme4]{lmerMod}}
#' @param  alpha significance level (Type I error probability)
#' @param ...    additional options that depends on the class of model
#'               passed to the function
#' @return an object of class \code{\link[paml]{pamlpowermix}}
#' @details Exactly one of the parameters n,
#'  power, alpha must be passed as NULL, and that parameter
#'  is determined from the others.
#'  Notice that alpha has non-NULL default so NULL must be
#'  explicitly passed if you want it computed.
#' @author Marcello Gallucci
#' @export


required_n<-function(model,
                     expand,
                     direction="within",
                     nsim=100,
                     effects="fixed",
                     target=NULL,
                     alpha=.05,
                     power=.80,
                     tol=.05,
                     method="sim",
                     seed=NULL) {

   time<-Sys.time()

   data<-model@frame

   if (!(expand %in% names(data)))
       stop("The variable",expand,"to expand the sample is not present in the model data.frame")
   if (direction=="within") {
        original_n<-round(mean(table(data[[expand]])))
        fun<-function(n) simr::extend(model,within=expand,n=n)
        msg0<-"Original "
        msg1<-paste(expand,"clusters mean size ")
   }
   else {
        original_n<-length(table(data[[expand]]))
        fun<-function(n) simr::extend(model,along=expand,n=n)
        msg0<-"Original "
        msg1<-paste(expand,"number of clusters ")
   }
   .model<-model
   n<-original_n
   keep<-TRUE
   i<-1



   cat("\nPower estimates...\n")
   while (keep) {
       .model<-fun(n)
        message(msg0, " ",msg1," ",n)
        msg0<-"Trying: "
        estimate<-modelpower(.model,nsim=nsim,effects = effects,alpha=alpha,seed=seed,method=method)
        levels<-c(units=estimate$N,estimate$clusters)
        bpower<-as.data.frame(estimate)
        if (!is.null(target) && !(target %in% (rownames(bpower))) )
           stop("Target coefficient",target,"is not present in the model")


       .bpower<-bpower

        if (!is.null(target))
           .bpower<- .bpower[rownames(.bpower)==target,]
        else {
           .bpower<- .bpower[row.names(.bpower)!="(Intercept)",]
           .bpower<- .bpower[which.min(.bpower$power),]
         }

        tn<-which.min(.bpower$power)
       .target<-row.names(.bpower)[tn]
        mp<-max(min(.bpower$power),.05)
        message("Found power ",mp," with ",expand," at ",n,"\n")


        n<-ceiling(power*n/(mp))
       if (.bpower$power>=power) {
         if ((.bpower$power-power)<tol)
                 keep=FALSE
         else
            if (n<original_n) {
                 keep=FALSE
                 message(paste("The input sample is larger",
                       "than the sample required for power ",
                        power,"of effect ",row.names(.bpower),
                       ". Consider using a smaller input sample for minimum required sample size."))
        }
       }
       estimate$target<-.target
       estimate$time<-as.numeric(Sys.time()-time)
       print(estimate,verbose=FALSE)
       cat("\n")
   }
  return(estimate)
}
