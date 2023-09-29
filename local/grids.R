
clusters<-list(ciao1=c(3,4),ciao2=c(3,4,5))
fixed<-list(ad1=c(3,4),af2=c(3,4))
random<-list(af2=c(1,2))

expand.list<-function(...) {
 opt<-c(...)
 optl<-list(...)
 items<-names(optl)
 eg<-expand.grid(opt)
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
 results
}
random<-NULL
q<-expand.list(clusters=clusters,fixed=fixed,random=random)
q
lists.atIndex<-function(i,alist) {
  lapply(alist,funtion(x) sapply(x,function(z) z[i]))
}
lists.atIndex()
eg<-expand.grid(c(clusters,fixed,random))
q1<-as.list(eg[,names(fixed)])
q2<-as.list(eg[,names(clusters)])
q2
q3<-as.list(eg[,names(random)])
if (length(names(q3))==0) {
         q3<-list(unlist(q3))
         names(q3)<-names(random)
}

