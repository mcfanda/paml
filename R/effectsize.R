#' Compute standardized coefficients
#'
#' Compute the standardized coefficients, betas, for a mixed model.
#'
#' @param  obj of class \code{\link[lme4]{merMod-class}}
#' @return a data.frame with the original coefficients and the beta coefficients
#' @author Marcello Gallucci
#' @examples
#'  methods(class="pamlpowermix")
#' @export
betas<- function(model,...) UseMethod("betas")

#' @export
betas.default<-function(model,...) {
  stop("No standardized coefficients for models of class ",class(model))
  }

#' @export
betas.lm<-function(model) {

  .data<-model$model
  vars<-as.list(attr(terms(mod),"predvars"))[-1]
  for (x in vars) .data[[x]]<-as.numeric(scale(.data[[x]]))
  coef(update(mod,data=.data))

}

#' @export
betas.lmerMod <- function(model, method="coef" ,verbose=FALSE) {

  data <- model@frame
  dep <- formula(model)[[2]]
  effects <- fixef(model)
  terms<- terms(model)
  preds <- attr(terms,"term.labels")
  clusters <- names(lme4::getME(model, "cnms"))


    for (pred in preds) {
        found = FALSE
        type  <- class(data[[pred]])
        if (verbose) message("Scaling ", pred, " as ", type)

        for (cluster_i in 1:length(clusters)) {

          cluster<-clusters[cluster_i]
          nested_test<-FALSE
          if (cluster_i>1) {nested_test<-lme4::isNested(data[[clusters[cluster_i-1]]],data[[cluster]])}
          if (length(grep(":", cluster))>0 & !(cluster %in% names(data)) ) {
            cls <- stringr::str_split(cluster, ":")[[1]]
            data[[cluster]]<-0
            n<-dim(data)[1]
            e<-floor(log10(n))
            e<-10^e
            k<-length(cls)
            for (i in 1:k)
              data[cluster]<-data[cluster]+e^(k-i+1)*as.numeric(data[[cls[i]]])
          }
          .data<-data

            if (nested_test && type!="factor") {
                 if (verbose) message("variable ", pred," aggregated over ",clusters[cluster_i-1]," nested within ",clusters[cluster_i])
                .data[[pred]]<-levelcenter(data[[pred]],data[[clusters[cluster_i-1]]],level = "between")
            }
            numerator <- .numerator(.data[[pred]],.data[[cluster]])

            if (!isFALSE(numerator)) {
                  if (!found) {
                    if (is.factor(.data[[pred]]))
                      pred<-paste0(pred,names(numerator))

                      if (nested_test) {
                        .data[[dep]]<-levelcenter(data[[dep]],data[[clusters[cluster_i-1]]],level = "between")
                      }

                    found = TRUE

                    y<-levelcenter(.data[[dep]],data[[clusters[cluster_i]]],level = "within")
                    if (method=="coef") {
                          denominator <- sd(y)
                          effects[pred] <- effects[pred] * numerator / denominator
                    }

                  if (verbose) message("Term ",paste(pred,callapse=", ")," scaled within cluster ", cluster)
                  if (verbose) message("with sd(y)=",paste(denominator,collapse = ", ")," and sd(x)=",paste(numerator,collapse =", " ))

                  }
            }
          }

        if (!found) {
           numerator <- .between_numerator(data[[pred]], data[[cluster]])
           if (is.factor(.data[[pred]]))
                pred<-paste0(pred,names(numerator))

           denominator <- sd(levelcenter(data[[dep]],data[[clusters[cluster_i]]],level = "between"))
           if (verbose) message("Term ",paste(pred,callapse=", ")," scaled between cluster ", cluster)
           if (verbose) message("with sd(y)=",paste(denominator,collapse = ", ")," and sd(x)=",paste(numerator,collapse =", " ))

           effects[pred] <- effects[pred] * numerator / denominator
        }
    }
    effects
  }




######### helper functions

.numerator<-function(var,cluster) {

    if (is.factor(var)) {
          # test that var varies within cluster
          test<-tapply(var, cluster, function(x) length(unique(x)))
          test2<-sum(test<2)/length(test)
          if (test2>.95)
                return(FALSE)
          ## if it varies, compute the numerator from the contrasts codes
          a<-contrasts(var)
          k<-nrow(a)
          N<-length(var)
          val<-sqrt(apply(a-mean(a),2,function(x) sum(x^2))*N/(k*(N-1)))
          } else {
          val<-sd(levelcenter(var,cluster))
          if (val<.00001) val<-FALSE
    }

    return(val)

}

.between_numerator<-function(var,cluster) {

  if (is.factor(var)) {

    n<-length(unique(cluster))
    a<-contrasts(var)
    k<-nlevels(var)
    mm<-apply(a,2,mean)
    b<-a-mm
    s<-apply(b,2,function(x) sum(x^2))
    val<-sqrt(s*(n/k)*(1/(n-1)))
  } else {
    val<-sd(levelcenter(var,cluster,level="between"))
  }

  return(val)

}
