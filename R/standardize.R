#' @export

levelcenter<-function(var,cluster,level="within") {

  data<-data.frame(x=var)
  data$id<-1:dim(data)[1]
  data$cluster<-as.character(cluster)
  data<-data[order(data$cluster),]
  m0<-tapply(data$x, list(data$cluster), mean,na.rm=TRUE)
  mm<-as.data.frame(cbind(names(m0),m0))
  names(mm)<-c("cluster","mm.")
  zdata<-merge(data,mm,by="cluster")
  zdata$mm.<-as.numeric(zdata$mm.)
  if (level=="within") {
    zdata$x<-zdata$x-zdata$mm.
  } else {
    zdata$x<-zdata$mm
  }
  zdata<-zdata[order(zdata$id),]
  zdata$x
}


## this produce a frame of codes for the nested clusters
## the mother cluster is 100,200,300 ..., the child
## will be 1001 1002 ... 2001, 2002 ...

.nestedclusters=function(clusters) {

  eg<-rev(expand.grid(rev(clusters)))
  alist<-lapply(1:(ncol(eg)), function(i) {
    a<-eg[,i]*1000^(ncol(eg)-i)
    a
  }  )
  eg<-do.call(cbind,alist)
  res<-eg
  for (i in 1:nrow(eg))
    for (j in 1:ncol(eg))
      res[i,j]<-sum(eg[i,1:j])/(10^(ncol(eg)-j))

  res
}



.rdraw<-function(sigma,varcorr) {

  ncol<-ncol(varcorr)
  res<-sigma
  for (i in 1:ncol)
    for (j in 1:ncol) {
      if (j==i) next
      s<-sqrt(varcorr[i,j])
      res[i,j]<-rnorm(1,sigma[i,j],s)
    }
  res
}
