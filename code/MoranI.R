MoranI <- function(sgrid,data,dclasses){

# Guillaume Larocque 2014
  
  D <- dist(sgrid)
  Dmat<-as.matrix(D)
  if (length(dclasses)==1){ #Based on equal frequencies
    sd<-sort(D)
    ngoal<-floor(length(sd)/dclasses)
    dcl<-0
    for (i in c(1:dclasses)){
      if (i==dclasses){
        tt<-tail(sd,n=1);
      }else{
        tt<-sd[ngoal];
      }
      dcl<-cbind(dcl,sd[tail(which(sd==tt),1)]);
      sd<-sd[-(1:ngoal)];
    }
    dclasses<-dcl;
  }
  cls<-t(dclasses);
  nd<-length(dclasses);
  if (!is.matrix(data)){
    data=as.matrix(data)
  }
  n<-nrow(data);
  Nh<-I<-cl<-VarI<-Z<-prob<-rep(0,length=nd-1)
  
  # Significance testing
  for (i in (1:(nd-1))){
    iHH<-(Dmat>dclasses[i]) & (Dmat<=dclasses[i+1])
    md<-mean(data[,1])
    prodd<-(t(data[,1]-md) %*% iHH %*% (data[,1]-md))/2
    Nh[i]<-sum(iHH>0)/2
    I[i]<-(1 / Nh[i])*prodd / (var(data[,1])*(n-1)/n)
    cl[i]<-mean(Dmat[iHH>0])
    Nh[i]<-Nh[i]*2
    S1<-Nh[i]*2
    S2<-sum((2*apply(iHH,1,sum))^2)
    b2<-n*sum((data[,1]-mean(data[,1]))^4)/(sum((data[,1]-mean(data[,1]))^2)^2)
    t1<-n*((n^2-3*n+3)*S1-n*S2+3*Nh[i]^2)
    t2<-b2*((n^2-n)*S1-2*n*S2+6*Nh[i]^2)
    t3<-(n-1)*(n-2)*(n-3)*(Nh[i]^2)
    VarI[i]<-((t1-t2)/t3)-solve((n-1)^2)
    if (Nh[i]>(4*(n-sqrt(n))) & Nh[i] <= (4*(2*n-3*sqrt(n)+1))) {
      if (I[i]<0){
        probs<-0.5
        tt=1
        while (tt>0.000001){
          Z<-(I[i]+(sqrt(10*probs)*solve(n-1)))/sqrt(VarI[i])
          pr<-pnorm(Z)
          tt<-abs(probs-pr)
          probs<-pr
        }
        prob[i]<-probs
      }else{
        probs<-0.5
        tt<-1
        while (tt>0.000001){
          Z<-(I[i]+(sqrt(10*probs)*solve(n-1)))/sqrt(VarI[i])
          pr<-pnorm(-Z)
          tt<-abs(probs-pr)
          probs<-pr
        }
        prob[i]<-probs;
      }
    } else {
      Z[i]<-(I[i]+solve(n-1))/sqrt(VarI[i])
      if (I[i]<0){
        prob[i]<-pnorm(Z[i]);
      }else{
        prob[i]<-pnorm(-Z[i]);
      }
    }
  }
  setClass("moranI",slots=c(distances="vector", 
                            mI="vector", 
                            number_pairs="vector",
                            prob.adjust="vector",
                            varI="vector"
  )
  )
  error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
  
  out=new("moranI",distances=cl,mI=I,number_pairs=Nh/2,prob.adjust=p.adjust(prob,'holm'),varI=VarI)
  setMethod("plot", "moranI",  function (x, xlab ="" , ylab ="" , axes = FALSE , asp =1 ){
    maxd=max(x@distances)
    maxy=max(x@mI)
    plot(x@distances,x@mI,xlab='Distance',ylab='Moran\'s I',xlim=c(0, maxd+maxd/20), ylim=c(-0.15, 1.15))
    title('Moran\'s I')
    lines(c(0,maxd+maxd/10),c(0,0))
    error.bar(x@distances,x@mI, 1.96*sqrt(x@varI))
  })
  
  setMethod("show", "moranI",  function (object) {
    print(data.frame(distances=cl,mI=I,number_pairs=Nh/2,prob.adjust=p.adjust(prob,'holm'),varI=VarI))
  })
  return(out)
}