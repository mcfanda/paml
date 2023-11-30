
s_data<-function() return(structure(Design$new(),class=c("sim_data","sim_object")))

s_integer<-function(x,mean=0,sd=1) structure(list(name=x,type="integer",mean=mean,sd=sd),class=c("sim_variable", "sim_integer","sim_object","singleton"))

s_numeric<-function(x,mean=0,sd=1) structure(list(name=x,type="numeric",mean=mean,sd=sd),class=c("sim_variable", "sim_numeric","sim_object","singleton"))

s_factor<-function(x,levels) {
   if (missing(levels))
      stop("factor variables should have number of levels defined with argument `levels`")
   fanthom_layer<-stringi::stri_rand_strings(1,length = FANTHOM)
   structure(list(name=x,type="factor",levels=levels,layer=-99,fanthom_layer=fanthom_layer,nested=NULL),class=c("sim_variable", "sim_factor","sim_object","singleton"))
}
s_cluster<-function(x,levels=NULL,layer=99) {
  fanthom_layer<-stringi::stri_rand_strings(1,length = FANTHOM)
  structure(list(name=x,levels=levels,layer=layer,fanthom_layer=fanthom_layer),class=c("sim_cluster","sim_object","singleton"))
}
s_dependent<-function(x,family=gaussian(),levels=NULL) structure(list(name=x,family=family,levels=levels),class=c("sim_dependent","sim_object","singleton"))

s_params<-function(x) structure(x,class=c("sim_params","sim_object","singleton"))

s_model<-function(x) structure(x,class=c("sim_model","sim_object","singleton"))


`+.sim_object`<-function(x,y) {

  if (!inherits(x,"sim_object"))
    stop("object ",y," is not of class sim_object" )
  if (!inherits(y,"sim_object"))
    stop("object ",y," is not of class sim_object" )

  if (inherits(x,"sim_variable")) {
    if (inherits(x,"singleton"))   x<-list(x)
    if (inherits(y,"singleton")) y<-list(y)
    vars<-c(x,y)
   return(structure(vars,class=c("sim_variables","sim_variable","sim_object")))
  }
  if (inherits(x,"sim_data")) {
    if (inherits(y,"singleton")) y<-list(y)
    lapply(y, function(v)  x$add_object<-v)
    return(x)

  }
}



`|.sim_variable`<-function(x,y) {

  if (!inherits(x,"sim_variable") || !inherits(y,"sim_cluster"))
    stop("For 'sim_object' object, the `|` operator requires 'sim_variable|sim_cluster' operations")


  if (!inherits(y,"sim_object"))
    stop("object ",y," is not of class sim_object" )
  if (!inherits(x,"sim_variables"))
    x<-list(x)
  if (!inherits(y,"sim_clusters"))
    y<-list(y)

  layer<-min(unlist(sapply(y,function(v) v$layer)))
  for (i in seq_along(x)) x[[i]]$layer<-layer
  structure(c(x,y),class=c("sim_data","sim_cluster","sim_object"))

}


`*.sim_object`<-function(x,y) {

  classes<-class(x)

  if (inherits(x,"singleton"))
    x<-list(x)
  if (inherits(y,"singleton"))
    y<-list(y)
  vars<-c(x,y)

  check<-unlist(sapply(vars, function(v) v$name))
  if (length(check)!=length(unique(check)))
    stop("Crossing the same variable variable is not meaningful.")

  fanthom_layer<-stringi::stri_rand_strings(1,length = FANTHOM)
  for (i in seq_along(vars))
            vars[[i]][['fanthom_layer']]<-fanthom_layer

  if (inherits(x,"sim_cluster"))
    return(structure(vars,class=c("sim_clusters","sim_cluster", "sim_object")))

  if ("sim_cluster" %in% classes)
    return(structure(vars,class=c("sim_cluster", "sim_clusters","sim_object")))
  if ("sim_variable" %in% classes)
    return(structure(vars,class=c("sim_variable", "sim_variables","sim_object")))

}

`/.sim_object`<-function(x,y) {

  classes<-class(x)

  if (inherits(x,"singleton"))
    x<-list(x)
  if (inherits(y,"singleton"))
    y<-list(y)
  vars<-c(x,y)

  vars<-rev(vars)

  check<-unlist(sapply(vars, function(v) v$name))
  if (length(check)!=length(unique(check)))
    stop("Nesting the same variable is not meaningful.")

  ufl<-unique(unlist(sapply(vars,function(v) v$fanthom_layer)))
  layers<-1:length(ufl)

  for (i in layers) {
     for (j in seq_along(vars))
        if (vars[[j]]$fanthom_layer==ufl[[i]]) {
              vars[[j]]$layer<-i
              }
  }


  for (i in layers[-length(layers)])
    for (j in seq_along(vars)) {
      .vars<-rlist::list.filter(vars,layer==(i+1))
      if (vars[[j]]$layer==i)
           vars[[j]]$nested<-list(rlist::list.select(.vars,name))
    }

  vars<-rev(vars)
  if ("sim_cluster" %in% classes)
        return(structure(vars,class=c("sim_cluster", "sim_clusters","sim_object")))
  if ("sim_variable" %in% classes)
       return(structure(vars,class=c("sim_variable", "sim_variables","sim_object")))
}



