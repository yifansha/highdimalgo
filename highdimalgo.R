library(randomForestSRC)
load("/Users/yifansha/Downloads/expression.RData")
expression1<-expression[1:2000,]
data<-expression1
data<-t(data)
nrep<-10
n<-nrow(data)
p<-ncol(data)
var.columns<-2:(p+1)
mvars<-ceiling(ncol(data)/5)
t<-vector("list",nrep)
resample<-function(x,size,...){
  if(length(x)<=1) {
    if (!missing(size) && size == 0) x[FALSE] else x
  }
  else {
    sample(x,size,...)
  }
}
colnames(data)<-paste("x",1:p,sep="")
rownames(data)<-paste("id",1:n,sep="")
set.seed(34)
f.active<-sample(seq(p),4,replace=FALSE)
result<-0.0000008*data[,f.active[1]]*data[,f.active[2]]*data[,f.active[3]]*data[,f.active[4]]
data<-data.frame(result,data)
names(data)[1]<-c("y")
formula<-y~.
options(rf.cores=1)
options(mc.cores=1)
sim.out<-rfsrc(formula,data)
sim.int<-find.interaction(sim.out,method="maxsubtree")
sim.out<-rfsrc(formula,data,importance=TRUE)
if (is.factor(data[,1])){
  wts<-pmax(sim.out$importance[,1],0)
} else {
  wts<-pmax(sim.out$importance,0)
}
sort(wts,decreasing=TRUE)[1:20]

for(m in 1:nrep){
  if (any(wts>0)) {
    var.pt<-unique(resample(var.columns,mvars,replace=FALSE,prob=wts))
  }
  else {
    var.pt<-unique(resample(var.columns,mvars,replace=FALSE))
  }
  data1<-data[,c(1,var.pt)]
  wts1<-wts[var.pt]
  t[[m]]<-match(names(data1[,md.varsel(formula,data1,wts1)+1]),names(data[,-1]))
}

md.varsel<-function(formula,data,wts){
  sim.out<-rfsrc(formula,data,predictorWt=wts)
  sim.int<-find.interaction(sim.out,method="maxsubtree")
  sim.int1<-sim.int
  for (i in 1:5){
    sim.intcopy<-sim.int1
    diag(sim.intcopy)<-NA
    o.max<-ncol(data)-1
    order<-c()
    for (i in 1:o.max){
      splinefit<-smooth.spline(seq(ncol(data)-2),sim.int1[i,-i])$y
      splinegap<-abs(sim.int1[i,-i]-splinefit)
      optk<-which.max(splinegap[1:o.max])
      order<-c(order,optk)
    }
    k<-quantile(order,probs=0.5)
    wt<-apply(sim.intcopy,1,function(x) mean(sort(x)[1:k],na.rm=TRUE))
    wt<-wt*diag(sim.int1)
    idx<-match(colnames(data[,-1]),names(wt))
    wt<-wt[idx]
    sim.out<-rfsrc(formula,data,predictorWt=log(1/wt),statistics=TRUE)
    sim.int1<-find.interaction(sim.out,method="maxsubtree")
  }
  
  ##var selection
  sim.intcopy<-sim.int1
  diag(sim.intcopy)<-NA
  sim.v<-as.vector(sim.intcopy)
  sim.sd<-sd(sim.v,na.rm=TRUE)
  thsd<-mean(sim.v,na.rm=TRUE)-0.15*sim.sd
  perrow<-apply(sim.intcopy,1,function(x) mean(x,na.rm=TRUE))
  print(sort(perrow,decreasing=FALSE))
  var.imp<-which(perrow<=thsd)
  idx<-match(names(var.imp),colnames(data[,-1]))
}

tt<-unique(unlist(t))
sim.out<-rfsrc(formula,data[,c(1,tt+1)])
sim.int<-find.interaction(sim.out,method="maxsubtree")

#####iRF try
idx<-sample(seq(p),500-4,replace=FALSE)
datasub<-data[,c(1,idx+1,f.active+1)]
