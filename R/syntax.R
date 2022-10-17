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

    coefs<-stringr::str_remove_all(stringr::str_extract_all(rhs,"\\[(.*?)\\]")[[1]],"[\\[\\]]")
    strformula<-stringr::str_remove_all(rhs,"\\[(.*?)\\]\\*")
    terms<-stringr::str_split(strformula,"\\+")[[1]]
#    pos<-str_split(rhs,"\\+")[[1]]
    if (length(terms)!=length(coefs)) stop("All terms must have a coefficient value")

    check<-stringr::str_split(stringr::str_extract_all(rhs,"\\+\\(.+\\)")[[1]],"\\+")[[1]]
    check<-check[check!=""]
    if (!all(stringr::str_detect(check,"[\\(\\)]"))) stop("Fixed effects should be specified before the random effects")

    formula<-as.formula(paste0(dep,"~",strformula))
    fixed<-nobars(formula)
    terms<-attr(terms(fixed),"term.labels")
    nfixed<-length(terms)

    check<-length(grep("\\*",fixed))>0
    if (check) stop("Interaction should be explicitely defined with the ':' operator")
    check<-length(grep("^1",unlist(stringr::str_split(fixed,"\\+"))))==0
    if (check) stop("Please explicitly specify the fixed intercept value using `~[value]*1+..`")
    check<-length(grep("^0",unlist(stringr::str_split(fixed,"\\+"))))>0
    if (check) stop("Zero intercepts models are not allowed, but one can specify a data generating model with intercept equal to zero with the syntax `~[0]*1+..`")

    fixedcoefs<-as.numeric(coefs[1:(nfixed+1)])
    names(fixedcoefs)<-c("(Intercept)",terms)

    randcoefs<-coefs[(nfixed+2):length(coefs)]
    rands<-findbars(formula)
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
              fixed=fixedcoefs,
              random=rev(randoms),
              dep=dep
      ))
}
