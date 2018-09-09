library(randomForestSRC)

datapr<-function(type,seed,order,beta){
  if (type=="expression"){
    load("/Users/yifansha/Downloads/expression.RData")
    expression1<-expression[1:2000,]
    data<-expression1
    data<-t(data)
    n<-nrow(data)
    p<-ncol(data)
    colnames(data)<-paste("x",1:p,sep="")
    rownames(data)<-paste("id",1:n,sep="")
    set.seed(seed)
    f.active<-sample(seq(p),order,replace=FALSE)
    result<-exp(beta*apply(data[,f.active], 1, prod))
    data<-data.frame(result,data)
    names(data)[1]<-c("y")
    formula<<-y~.
    data<-data
  }
}

w.initial<-function(formula,data){
  sim.out<-rfsrc(formula,data,importance=TRUE)
  if (is.factor(data[,1])){
    pmax(sim.out$importance[,1],0)
  } else {
    pmax(sim.out$importance,0)
  }
}

resample<-function(x,size,...){
  if(length(x)<=1) {
    if (!missing(size) && size == 0) x[FALSE] else x
  }
  else {
    sample(x,size,...)
  }
}

bd.kink<-function(data,sim.int){
  o.max<-ncol(data)-1
  order<-sapply(seq(o.max), FUN=function(x){ 
    which.max(abs(sim.int[x,-x]-smooth.spline(seq(o.max-1),sim.int[x,-x])$y)[1:o.max])})
  quantile(order,probs=0.5)
}

wt.itr<-function(i,data,m){
  if (i==1){
    sim.int1<-m
    sim.int1<<-sim.int1
  }
  else {
    sim.intcopy<<-sim.int1
    k<-bd.kink(data,sim.intcopy)
    diag(sim.intcopy)<-NA
    wt<-apply(sim.intcopy,1,function(x) mean(sort(x)[1:k],na.rm=TRUE))
    wt<<-wt*diag(sim.int1)
    idx<-match(colnames(data[,-1]),names(wt))
    wt<-wt[idx]
    sim.out<-rfsrc(formula,data,predictorWt=log(1/wt),statistics=TRUE)
    sim.int1<<-find.interaction(sim.out,method="maxsubtree")
  }
  list(sim.int1)
}

ett<-function(t,data){
  sapply(seq(5),wt.itr,data=data,m=t)
}

init.gen<-function(x,wts,data){
  if (any(wts>0)) {
    var.pt<-unique(resample(var.columns,mvars,replace=FALSE,prob=wts))
  }
  else {
    var.pt<-unique(resample(var.columns,mvars,replace=FALSE))
  }
  data1<-data[,c(1,var.pt)]
  wts1<-wts[var.pt]
  sim.out<-rfsrc(formula,data1,predictorWt=wts1)
  sim.int<-find.interaction(sim.out,method="maxsubtree")
  sim.int1<-sim.int
  list(data1,sim.int1)
}

md.varsel<-function(mt,data){
  sim.intcopy<-mt
  diag(sim.intcopy)<-NA
  perrow<-apply(sim.intcopy,1,function(x) mean(x,na.rm=TRUE))
  thsd<-quantile(perrow,probs=0.1)
  print(sort(perrow,decreasing=FALSE))
  var.imp<-which(perrow<=thsd)
  idx<-match(names(var.imp),colnames(data[,-1]))
}

var.sel<-function(data,wts){
  init<-lapply(seq(nrep),init.gen,wts=wts,data=data)
  init.m<-lapply(init,"[[",2)
  init.dta<-lapply(init,"[[",1)
  mt.list<-mapply(ett,init.m,init.dta,SIMPLIFY=FALSE)
  mt.final<-lapply(mt.list, '[[', 5)
  list.idx<-mapply(md.varsel,mt.final,init.dta,SIMPLIFY=FALSE)
  mapply(function(x,data1){match(names(data1[,x+1]),
                                               names(data[,-1]))}
         ,list.idx,data1=init.dta,SIMPLIFY=FALSE)
}

data<-datapr("expression",62,4,0.00003)
wts<-w.initial(formula,data)
nrep<-2
var.columns<-2:(ncol(data))
mvars<-ceiling((ncol(data)-1)/5)

t<-var.sel(data,wts)

tt<-unique(unlist(t))
sim.out<-rfsrc(formula,data[,c(1,tt+1)])
sim.int<-find.interaction(sim.out,method="maxsubtree")


