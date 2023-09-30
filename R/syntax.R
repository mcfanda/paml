#
#  Functions to extract model formula and input parameters from
#  a string of type `y~[1]*1+[2]*x+([10]*1|cluster)`
#
#
#' @export

get_formula_info<-function(astring) {

   if (inherits(astring,"formula"))
      astring<-paste(deparse(form, width.cutoff = 500), collapse="")

   form0<-stringr::str_split_fixed(astring,"~",n=Inf)
   rhs<-form0[[2]]
   dep<-form0[[1]]
   clean<-.make_raw(rhs)
   fixed<-stringr::str_split(stringr::str_remove_all(clean,"\\((.*?)\\)"),"\\+")[[1]]
   fixed<-fixed[unlist(sapply(fixed,function(x) x!=""))]
   fixed<-sapply(fixed, function(x) {
     coef<-as.numeric(stringr::str_remove_all(stringr::str_extract_all(x,"\\[(.*?)\\]"),"[\\[\\]]"))
     terms<-gsub("\\*","",stringr::str_remove_all(x,"\\[(.*?)\\]"))
     terms<-stringr::str_replace_all(terms,"^1","(Intercept)")
     names(coef)<-terms
     coef
   },USE.NAMES = FALSE)

   ### first we deal with the random component
   .random<-stringr::str_extract_all(clean,"\\((.*?)\\)")[[1]]
   clusters<-list()
   clusterstag<-list()
   random<-lapply(as.list(.random), function(x) {
     str<-stringr::str_split(x,"\\|",n=Inf)[[1]]
     cl<-gsub("\\)","",str[[length(str)]])
     clname<-cl
     if (cl %in% clusters)
       clname<-paste(cl,sum(as.numeric(clusters==cl)),sep=".")
     else
       clusters[[length(clusters)+1]]<<-clname
     ladd(clusterstag)<<-clname
     str<-stringr::str_split(x,"\\+",n=Inf)[[1]]
     suppressWarnings(
       coefs<-sapply(str, function(z) as.numeric(stringr::str_remove_all(stringr::str_extract_all(z,"\\[(.*?)\\]"),"[\\[\\]]")))
     )
     coefs<-coefs[!is.na(coefs)]
     str<-names(coefs)
     terms<-gsub("\\[(.*?)\\]","",str)
     terms<-gsub("[\\(\\)\\*]","",terms)
     terms<-gsub(paste0("\\|",cl),"",terms)
     terms<-stringr::str_replace_all(string = terms,pattern = "^1","(Intercept)")
     names(coefs)<-paste(clname,terms,sep=".")
     names(coefs)<-terms

     coefs
   })
   names(random)<-clusterstag

   formula<-stringr::str_remove_all(rhs,"\\[(.*?)\\]")
   formula<-stringr::str_remove_all(formula,"\\*")
   formula<-as.formula(paste(dep,formula,sep="~"))
   structure(list(
              formula=  formula,
              fixed=fixed,
              random=random,
              clusters=clusters,
              dep=dep
      ))
}

.make_raw<-function(rhs) {
  msg<-NULL
  pre<-stringr::str_split(rhs,"[\\+\\|~]",n=Inf)[[1]]
  g<-lapply(pre,function(x) {
    ## we search for factors coefficients passed as [1,2]*factor
    v<-stringr::str_extract_all(x,"\\[(.*?)\\,(.*?)\\]")[[1]]
    if (length(v)>0) {
      ## if found, we split in dummies
      n<-stringr::str_split(x,"\\*")[[1]][[2]]
      q<-as.numeric(stringr::str_split(gsub("\\[|\\]","",v),"\\,")[[1]])
      x<-paste(paste0("[",q,"]"),paste0(n,".",seq_along(q)),sep = "*",collapse = "+")
    } else {
      ## here we remove the father/offspring notation and set it to [FIXIT]
      if (length(grep("\\)",x))>0) x<-paste0("|",x)
      ## then we check if something has no coefficient, and set to zero
      nocof<-length(grep("\\[(.*?)\\]",x,invert = T))>0
      nocl<-length(grep("\\|(.*?)\\)",x))==0
      nozero<-(length(grep("^\\(0",x))==0 && length(grep("^0",x))==0)
      if (nocof)
         if (nozero & nocl) {
              msg<<-"The expected value of terms without coefficients is set to 0."
              test<-grep("^\\(",x)
                if (length(test)>0)
                         x<-paste0("([0]*",gsub("\\(","",x))
                else
                         x<-paste0("[0]*",x)
          }
     }
    x

  })
  result<-gsub("+|","|",paste0(g,collapse = "+"), fixed = T)
  if (!is.null(msg)) attr(result,"warning")<-msg
  if (!is.null(msg)) attr(result,"zeros")<-TRUE

  result
}

