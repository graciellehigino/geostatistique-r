GearyC <- function(sgrid,data,dclasses){
  
  # Guillaume Larocque 2014
  
  D <- dist(sgrid)
  Dmat<-as.matrix(D)
  if (length(dclasses)==1){ #Based on equal frequencies
    sd<-sort(D)
    ngoal<-floor(length(sd)/dclasses)
    dcl<-0
    for (i in c(1:dclasses)){
      if (i==dclasses){
        tt<-tail(sd,n=1)
      }else{
        tt<-sd[ngoal]
      }
      dcl<-cbind(dcl,sd[tail(which(sd==tt),1)])
      sd<-sd[-(1:ngoal)]
    }
    dclasses<-dcl
  }
  cls<-t(dclasses)
  nd<-length(dclasses)
  if (!is.matrix(data)){
     data=as.matrix(data)
  }
  n<-nrow(data);
  VarC<-Z<-prob<-Nh<-cl<-rep(0,length=nd-1)
  
  data=data.frame(data)
  names(data)='data'
  SPDF <- SpatialPointsDataFrame(coords=sgrid, data=data)
  vario <- variogram(data~1,SPDF, boundaries=cls)
  gC=vario$gamma/var(data);
  # Significance testing
  for (i in (1:(nd-1))){
    iHH<-(Dmat>dclasses[i]) & (Dmat<=dclasses[i+1])
    cl[i]<-mean(Dmat[iHH>0])
    Nh[i]=vario$np[i]*2
    S1=Nh[i]*2
    S2<-sum((2*apply(iHH,1,sum))^2)
    b2<-n*sum((data[,1]-mean(data[,1]))^4)/(sum((data[,1]-mean(data[,1]))^2)^2)
    t1=(n-1)*S1*((n^2-3*n+3)-(n-1)*b2)
    t2=n*(n-2)*(n-3)*(Nh[i]^2)
    t3=-0.25*(n-1)*S2*(n^2+3*n-6-(n^2-n+2)*b2)+Nh[i]^2*(n^2-3+(-(n-1)^2)*b2)
    VarC[i]=t1/t2+t3/t2;
    if (gC[i]>1){
      if (Nh[i]>(4*(n-sqrt(n))) & Nh[i] <= (4*(2*n-3*sqrt(n)+1))) {
        probs<-0.5
        tt<-1
        while (tt>0.000001){
          Z<-(gC[i]+(sqrt(10*probs)*solve(n-1)))/sqrt(VarC[i])
          pr<-pnorm(Z)
          tt<-abs(probs-pr)
          probs<-pr
        }
        prob[i]<-probs
      }else{
        probs<-0.5
        tt<-1
        while (tt>0.000001){
          Z<-(gC[i]-1+solve(n-1))/sqrt(VarC[i])
          pr<-pnorm(-Z)
          tt<-abs(probs-pr)
          probs<-pr
        }
        prob[i]<-probs;
      }
    } else {
      Z[i]<-(gC[i]-1)/sqrt(VarC[i])
      prob[i]<-pnorm(Z[i]);
    }
  }
  setClass("gearyC",slots=c(distances="vector", 
                          gC="vector", 
                          number_pairs="vector",
                          prob.adjust="vector",
                          varC="vector"
           )
  )
  error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
  
  out=new("gearyC",distances=cl,gC=gC,number_pairs=Nh/2,prob.adjust=p.adjust(prob,'holm'),varC=VarC)
  setMethod("plot", "gearyC",  function (x, xlab ="" , ylab ="" , axes = FALSE , asp =1 ){
    maxd=max(x@distances)
    plot(x@distances,x@gC,xlab='Distance',ylab='Geary\'s C',xlim=c(0, maxd+maxd/20), ylim=c(0, 1.15))
    lines(c(0,maxd+maxd/10),c(1,1))
    title('Geary\'s C')
    error.bar(x@distances,x@gC, 1.96*sqrt(x@varC))
  })
  
  setMethod("show", "gearyC",  function (object) {
    print(data.frame(distances=cl,gC=gC,number_pairs=Nh/2,prob.adjust=p.adjust(prob,'holm'),varC=VarC))
  })
  return(out)
}