#
#  Functions to extract model formula and input parameters from
#  a string of type `y~[1]*1+[2]*x+([10]*1|cluster)`
#
#
#' @export

astring<-"y~[.10]*1+[.5,.5]*afac2+([1]*1+[1,.5,.3]*x1+[2]*x2|acluster)+([1]*1|bcluster)"

get_formula_info<-function(astring) {

    form0<-stringr::str_split_fixed(astring,"~",n=Inf)
    dep<-form0[[1]]
    rhs<-form0[[2]]
    coefs<-stringr::str_remove_all(stringr::str_extract_all(rhs,"\\[(.*?)\\]")[[1]],"[\\[\\]]")
    strformula<-stringr::str_remove_all(rhs,"\\[(.*?)\\]\\*")
    pre<-gsub("\\(","",gsub("\\|(.*?)\\)","",strformula))
    terms<-stringr::str_split(pre,"\\+",n=Inf)[[1]]

    formula<-as.formula(paste0(dep,"~",strformula))
    fixed<-lme4::nobars(formula)
    stringr::str_split(as.character(fixed[[3]]),"+")
    nfixed<-length(fixed)
    random<-lme4::findbars(formula)
    rterms<-.flat(lapply(random,function(x) stringr::str_split_fixed(as.character(x[[2]]),"\\+",n = Inf)))
    rterms<-rterms[sapply(rterms,function(x) x!="")]
    terms<-stringr::str_split(strformula,"\\+")[[1]]
#    pos<-str_split(rhs,"\\+")[[1]]
    if (length(terms)!=length(coefs)) stop("All terms must have a coefficient value")


    check<-stringr::str_split(stringr::str_extract_all(rhs,"\\+\\(.+\\)")[[1]],"\\+")[[1]]
    check<-check[check!=""]
    if (!all(stringr::str_detect(check,"[\\(\\)]"))) stop("Fixed effects should be specified before the random effects")
    test<-grep(",",coefs)
    print(terms)
    for (w in test) {
      new<-stringr::str_split(coefs[[w]],",")
      coefs[[w]]<-new
      coefs<-.flat(coefs)
      terms[[w]]<-list(paste0(unlist(terms[[w]]),".",1:length(unlist(new[[1]]))))
      terms<-.flat(terms)
    }
    print(terms)
    rawformula<-as.formula(paste0(dep,"~",paste(terms,collapse = " + ")))
    formula<-as.formula(paste0(dep,"~",strformula))
    fixed<-lme4::nobars(formula)
    terms<-attr(terms(fixed),"term.labels")
    nfixed<-length(terms)

    check<-length(grep("\\*",fixed))>0
    if (check) stop("Interaction should be explicitely defined with the ':' operator")
    check<-length(grep("^1",unlist(stringr::str_split(as.character(fixed),"\\+"))))==0
    if (check) stop("Please explicitly specify the fixed intercept value using `~[value]*1+..`")
    check<-length(grep("^0",unlist(stringr::str_split(as.character(fixed),"\\+"))))>0
    if (check) stop("Zero intercepts models are not allowed, but one can specify a data generating model with intercept equal to zero with the syntax `~[0]*1+..`")
    fixedcoefs<-as.numeric(coefs[1:(nfixed+1)])
    names(fixedcoefs)<-c("(Intercept)",terms)

    randcoefs<-coefs[(nfixed+2):length(coefs)]
    rands<-lme4::findbars(formula)
    randoms<-list()
    k<-1
    for (r in rands) {
             cluster<-r[[3]]
             terms<-unlist(stringr::str_split(as.character(r[[2]]),"\\+"))
             terms<-terms[terms!=""]
             terms<-stringr::str_replace(terms,"^1","(Intercept)")
            .coefs<-as.numeric(randcoefs[k:(k+length(terms)-1)])
             names(.coefs)<-terms
             randoms[[cluster]]<-.coefs
             k<-k+length(terms)
     }

    structure(list(
              formula=formula,
              rawformula=rawformula,
              fixed=fixedcoefs,
              random=rev(randoms),
              dep=dep
      ))
}
get_formula_info(astring)

.make_raw<-function(rhs) {

  pre<-stringr::str_split(rhs,"[\\+\\|~]",n=Inf)[[1]]
  g<-lapply(pre,function(x) {
    v<-stringr::str_extract_all(x,"\\[(.*?)\\,(.*?)\\]")[[1]]
    if (length(v)>0) {
      n<-stringr::str_split(x,"\\*")[[1]][[2]]
      q<-as.numeric(stringr::str_split(gsub("\\[|\\]","",v),"\\,")[[1]])
      x<-paste(paste0("[",q,"]"),paste0(n,".",seq_along(q)),sep = "*",collapse = "+")
    }
    x

  })
  paste0(g,collapse = "+")

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
