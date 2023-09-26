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

