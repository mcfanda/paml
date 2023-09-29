is.something <- function(x, ...) UseMethod(".is.something")

.is.something.default <- function(obj) (!is.null(obj))

.is.something.list <- function(obj) (length(obj) > 0)

.is.something.numeric <- function(obj) (length(obj) > 0)

.is.something.character <- function(obj) (length(obj) > 0 && obj!="")

.is.something.logical <- function(obj) !is.na(obj)

is.there<-function(pattern,string) length(grep(pattern,string,fixed=T))>0


`ladd<-`<-function(x,value,key=NULL) {
  if (is.something(key))
    x[[key]]<-value
  else
    x[[length(x)+1]]<-value
  return(x)
}

## we need these to have a future
`%dorng%`<-doRNG::`%dorng%`
`%dopar%`<- foreach::`%dopar%`

## some helper

num<-function(x) as.numeric(as.character(x))

vec_to_mat<-function(vec) {
   n<-.revcomb(length(vec))
   x<-diag(1,ncol=n,nrow = n)
   x[lower.tri(x,diag = T)]<-vec
   x[upper.tri(x)]<-x[lower.tri(x)]
   x
}

.revcomb<-function(k) {
  r=0
  n<-1
  while(k!=r) {
    n<-n+1
    r<-n+(n*(n-1)/2)
    if (n>k) {
      n<-NaN
      stop("No matrix dimensions can be found")
    }
  }
  n
}


expand.list<-function(...) {
  opt<-c(...)
  optl<-list(...)
  items<-names(optl)
  eg<-expand.grid(opt)
  nitems<-dim(eg)[1]
  result<-list()
  results<-sapply(items, function(x) {
    .names<-grep((paste0("^",x,"\\.")),names(eg))
    .sel<-as.list(eg[,.names])
    if (length(names(.sel))==0) {
      .sel<-list(unlist(.sel))
      names(.sel)<-names(optl[[x]])
    }
    names(.sel)<-gsub(paste0("^",x,"\\."),"",names(.sel))
    .sel
  },simplify = F)
  attr(results,"nitems")<-nitems
  results
}

lists.atIndex<-function(i,alist) {
  lapply(alist, function(x) sapply(x,function(z) z[i]))
}

