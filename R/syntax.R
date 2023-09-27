#
#  Functions to extract model formula and input parameters from
#  a string of type `y~[1]*1+[2]*x+([10]*1|cluster)`
#
#
#' @export


get_formula_info<-function(astring) {

    form0<-stringr::str_split_fixed(astring,"~",n=Inf)
    dep<-form0[[1]]
    rhs<-form0[[2]]
    rawformula<-.make_raw(rhs)
    rawterms<-stringr::str_split(rawformula,"\\+",n=Inf)[[1]]
    rawcoefs<-stringr::str_remove_all(stringr::str_extract_all(rawformula,"\\[(.*?)\\]")[[1]],"[\\[\\]]")
    test<-grep("\\(0$",rawterms)
    if (length(test)>0)
       rawterms<-rawterms[-test]

    if (length(rawterms)!=length(rawcoefs)) stop("All terms must have a coefficient value. Only 0 random intercepts are defined without a value.")

    merformula<-stringr::str_remove_all(rhs,"\\[(.*?)\\]\\*")
    rawformula<-stringr::str_remove_all(rawformula,"\\[(.*?)\\]\\*")

    formula<-as.formula(paste0(dep,"~",rawformula))
    fixed<-lme4::nobars(formula)
    terms<-attr(terms(fixed),"term.labels")
    nfixed<-stringr::str_count(as.character(fixed)[[3]],"\\+")+1
    check<-length(grep("\\*",fixed))>0
    if (check) stop("Interaction should be explicitely defined with the ':' operator")
    check<-length(grep("^1",unlist(stringr::str_split(as.character(fixed),"\\+"))))==0
    if (check) stop("Please explicitly specify the fixed intercept value using `~[value]*1+..`")
    check<-length(grep("^0",unlist(stringr::str_split(as.character(fixed),"\\+"))))>0
    if (check) stop("Zero fixed intercept models are not allowed, but one can specify a data generating model with intercept equal to zero with the syntax `~[0]*1+..`")
    fixedcoefs<-as.numeric(rawcoefs[1:(nfixed)])
    names(fixedcoefs)<-c("(Intercept)",terms)
    randcoefs<-rawcoefs[(nfixed+1):length(rawcoefs)]
    rands<-lme4::findbars(formula)
    randoms<-list()
    clustersname<-list()
    k<-1

    for (r in rands) {
             cluster<-as.character(r[[3]])
             terms<-unlist(stringr::str_split(as.character(r[[2]]),"\\+"))
             terms<-terms[terms!=""]
             terms<-stringr::str_replace(terms,"^1","(Intercept)")
             check<-grep("^0",unlist(stringr::str_split(as.character(terms),"\\+")))
             if (length(check)>0) terms<-terms[-check]
             clustersname<-unique(c(clustersname,cluster))
             .coefs<-as.numeric(randcoefs[k:(k+length(terms)-1)])
             names(.coefs)<-terms
             if (cluster %in% names(randoms)) {
               test<-grep(cluster,names(randoms))
               cluster<-paste(cluster,length(test),sep = ".")
             }
             randoms[[cluster]]<-.coefs
             k<-k+length(terms)
     }

    structure(list(
              formula=  as.formula(paste0(dep,"~",merformula)),
              rawformula=as.formula(paste0(dep,"~",rawformula)),
              fixed=fixedcoefs,
              random=randoms,
              clusters=clustersname,
              dep=dep
      ))
}

.make_raw<-function(rhs) {

  pre<-stringr::str_split(rhs,"[\\+\\|~]",n=Inf)[[1]]
  g<-lapply(pre,function(x) {
    v<-stringr::str_extract_all(x,"\\[(.*?)\\,(.*?)\\]")[[1]]
    if (length(v)>0) {
      n<-stringr::str_split(x,"\\*")[[1]][[2]]
      q<-as.numeric(stringr::str_split(gsub("\\[|\\]","",v),"\\,")[[1]])
      x<-paste(paste0("[",q,"]"),paste0(n,".",seq_along(q)),sep = "*",collapse = "+")
    } else {
      if (length(grep("\\)",x))>0) x<-paste0("|",x)
    }

    x

  })
  gsub("+|","|",paste0(g,collapse = "+"), fixed = T)

}

.flat<-function(alist) {

  y<-list()
  for (i in seq_along(alist)) {
    if (length(unlist(alist[[i]]))>1) {

       for (j in seq_along(unlist(alist[[i]])) )
           y[[length(y)+1]]<-unlist(alist[[i]])[[j]]
    }
    else
       y[[length(y)+1]]<-alist[[i]]

  }
 y

}
