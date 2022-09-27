# estimation
est.sim <- function(k=1,alpha=0.05,beta=0.2,rho=2,sigma=1,theta=1,K=1,f="normal",method="sar",des.f="normal",m=NULL,r,a,eps=0.001,bootstrap=0,nonpar=TRUE,B=1e4,ties="average",out.tab="b"){
  delta <- K/sigma
  # pi1 <- c(alpha*(1/k)^rho,diff(alpha*((1:k)/k)^rho))
  # pi2 <- c(beta*(1/k)^rho,diff(beta*((1:k)/k)^rho))
  if(f=="normal"){
    f <- function(n,mu=0,sig){rnorm(n=n,mean=mu,sd=sig)}
    ppsi <- function(x,mu=0,sig=1){pnorm(x,mean=mu,sd=sig)}
    dpsi <- function(x,mu=0,sig=1,lg=FALSE){dnorm(x,mean=mu,sd=sig,log=lg)}
    qpsi <- function(q,mu=0,sig=1) qnorm(q)*sig + mu
  } else if(f=="logistic"){
    f <- function(n,mu=0,sig){rlogis(n=n,location=mu,sig*sqrt(3)/pi)}
    ppsi <- function(x,mu=0,sig=1){plogis(x,location=mu,scale=sig*sqrt(3)/pi)}
    dpsi <- function(x,mu=0,sig=1,lg=FALSE){dlogis(x,location=mu,scale=sig*sqrt(3)/pi,log=lg)}
    qpsi <- function(q,mu=0,sig=1){ qlogis(q,location=mu,scale=sig*sqrt(3)/pi) }
  } else if(f=="laplace"){
    # uses package rmutil
    f <- function(n,mu=0,sig){rmutil::rlaplace(n=n,m=mu,s=sig/sqrt(2))}
    ppsi <- function(x,mu=0,sig=1){rmutil::plaplace(x,m=mu,s=sig/sqrt(2))}
    dpsi <- function(x,mu=0,sig=1,lg=FALSE){rmutil::dlaplace(x,m=mu,s=sig/sqrt(2),log=lg)}
    qpsi <- function(q,mu=0,sig=1){ rmutil::qlaplace(q,m=mu,s=sig/sqrt(2)) }
  } else if(f=="t"){
    f <- function(n,mu=0,sig){mu + sig*rt(n=n,df=3,ncp=0)/sqrt(3)}
    ppsi <- function(x,mu=0,sig=1){pt((x-mu)*sqrt(3)/sig,df=3,ncp=0)}
    dpsi <- function(x,mu=0,sig=1,lg=FALSE){
      if(lg==FALSE){ return(dt((x-mu)*sqrt(3)/sig,df=3,ncp=0)*sqrt(3)/sig)
      } else { dt((x-mu)*sqrt(3)/sig,df=3,ncp=0,log=T)+log(sqrt(3))-log(sig) }
    }
    qpsi <- function(q,mu=0,sig=1) qt(q,df=3)/sqrt(3)*sig + mu
  }
  # mixture dist functions
  rpsi.mix <- function(n,mu=0,sig=1,theta,delta) {
    ind <- sample(1:2,prob=c(1-theta,theta),size=10000,replace = TRUE)
    mus <- c(mu,mu+delta)
    f(n=n,mu=mus[ind],sig=sig)
  }
  dpsi.mix <- function(x,mu=0,sig=1,theta,delta) { (1-theta)*dpsi(x,mu=mu,sig=sig) + theta*dpsi(x,mu=mu+delta,sig=sig) }
  ppsi.mix <- function(x,mu=0,sig=1,theta,delta) { (1-theta)*ppsi(x,mu=mu,sig=sig) + theta*ppsi(x,mu=mu+delta,sig=sig) }
  qpsi.mix <- function(q,mu=0,sig=1,theta,delta) {
    temp <- function(x,q,mu=mu,sig=sig,theta=theta,delta=delta) (ppsi.mix(x,mu=mu,sig=sig,theta=theta,delta=delta) - q)^2
    optimize(f=temp,interval=c((mu+delta)+10*sig,mu-10*sig),q=q,mu=mu,sig=sig,theta=theta,delta=delta)$minimum
  }
  if(des.f=="normal"){
    des.f <- function(n,mu=0,sig){rnorm(n=n,mean=mu,sd=sig)}
    des.ppsi <- function(x){pnorm(x,mean=0,sd=1)}
    des.dpsi <- function(x,mu=0,sig=1,lg=FALSE){dnorm(x,mean=mu,sd=sig,log=lg)}
  } else if(des.f=="logistic"){
    des.f <- function(n,mu=0,sig){rlogis(n=n,location=mu,sig*sqrt(3)/pi)}
    des.ppsi <- function(x){plogis(x,location=0,scale=1*sqrt(3)/pi)}
    des.dpsi <- function(x,mu=0,sig=1,lg=FALSE){dlogis(x,location=mu,scale=sig*sqrt(3)/pi,log=lg)}
  } else if(des.f=="laplace"){
    # uses package rmutil
    des.f <- function(n,mu=0,sig){rmutil::rlaplace(n=n,m=mu,s=sig/sqrt(2))}
    des.ppsi <- function(x){rmutil::plaplace(x,m=0,s=1/sqrt(2))}
    des.dpsi <- function(x,mu=0,sig=1,lg=FALSE){rmutil::dlaplace(x,m=mu,s=sig/sqrt(2),log=lg)}
  } else if(des.f=="t"){
    des.f <- function(n,mu=0,sig){mu + sig*rt(n=n,df=3,ncp=0)/sqrt(3)}
    des.ppsi <- function(x,mu=0,sig=1){pt((x-mu)*sqrt(3)/sig,df=3,ncp=0)}
    des.dpsi <- function(x,mu=0,sig=1,lg=FALSE){
      if(lg==FALSE){ return(dt((x-mu)*sqrt(3)/sig,df=3,ncp=0)*sqrt(3)/sig)
      } else { dt((x-mu)*sqrt(3)/sig,df=3,ncp=0,log=T)+log(sqrt(3))-log(sig) }
    }
  }
  h1<-function(y,theta,K) {
    ppsi(y)*((1-theta)*dpsi(y) + theta*dpsi(y-K))
  }
  h2<-function(y,theta,K) {
    ((1-(1-theta)*ppsi(y)-theta*ppsi(y-K))^2)*dpsi(y)
  }
  h3<-function(y,theta,K) {
    (ppsi(y))^2*((1-theta)*dpsi(y) + theta*dpsi(y-K))
  }
  gamval<-integrate(h1,-Inf,Inf,theta,K)$value
  xi1<-integrate(h2,-Inf,Inf,theta,K)$value-gamval^2
  xi2<-integrate(h3,-Inf,Inf,theta,K)$value-gamval^2
  lambda <- 1/2
  m.approx <- ((qnorm(1-alpha)/sqrt(6) + qnorm(1-beta)*sqrt(xi1+xi2))/(gamval-1/2))^2
  ifelse(is.null(m), m <- ceiling(m.approx), m <- m)
  n <- m
  
  mu <- m*(m+m+1)/2
  v <- m*m*(m+m+1)/12
  
  if(method=="rerank"){w <- matrix(0,ncol = k,nrow=B); w.alt <- matrix(0,ncol = k,nrow=B)
  }else{srs <- matrix(0,ncol = k,nrow=B); srs.alt <- matrix(0,ncol = k,nrow=B)}
  stage <- vector("numeric", length=B)
  theta.hat <- vector("numeric",length=B)
  delta.hat <- vector("numeric",length=B)
  theta.hat.check <- delta.hat.check <- vector("numeric",length=B)
  theta.hat.deps <- delta.hat.deps <- vector("numeric",length=B)
  theta.hat.noeps <- delta.hat.noeps <- vector("numeric",length=B)
  theta.hat.adj <- delta.hat.adj <- vector("numeric",length=B)
  # theta.hat.mad <- delta.hat.mad <- vector("numeric",length=B)
  # theta.hat.mad.adj <- delta.hat.mad.adj <- vector("numeric",length=B)
  theta.hat.mle <- delta.hat.mle <- vector("numeric",length=B)
  theta.hat.km <- delta.hat.km <- vector("numeric",length=B)
  theta.hat.ckm <- delta.hat.ckm <- vector("numeric",length=B)
  theta.hat.ckmx <- delta.hat.ckmx <- vector("numeric",length=B)
  # theta.hat.pamxy <- delta.hat.pamxy <- vector("numeric",length=B)
  # theta.hat.mle.n <- delta.hat.mle.n <- vector("numeric",length=B)
  bc.theta.hat.p <- bc.delta.hat.p <- vector("numeric",length=B)
  bc.theta.hat.np <- bc.delta.hat.np <- vector("numeric",length=B)
  bc.theta.hat.mle <- bc.delta.hat.mle <- vector("numeric",length=B)
  # bc.theta.hat.mle.n1 <- bc.delta.hat.mle.n1 <- vector("numeric",length=B)
  # bc.theta.hat.mle.n2 <- bc.delta.hat.mle.n2 <- vector("numeric",length=B)
  out.r <- out.a <- NULL
  out <- vector("numeric",length=B)
  
  for(i in 1:B){
    ind <- sample(1:2,prob=c(1-theta,theta),size=k*m,replace = TRUE)
    mus <- c(0,delta)
    y <- f(n=k*m,mu=mus[ind],sig=sigma)
    obs.alt <- matrix(c(f(n=k*m,sig=sigma),y),ncol=2,nrow=(k*m),byrow=F)
    ts.alt <- vector("numeric",length=k)
    for(i1 in 1:k){
      if(method=="rerank"){
        ts.alt[i1] <- sum(rank(obs.alt[1:(i1*m),])[(i1*m+1):(i1*m+i1*m)])
      } else if(method=="sar"){
        ts.alt[i1] <- sum(rank(obs.alt[((i1-1)*m+1):(i1*m),])[(n+1):(n+m)])
      } else if(method=="average"){
        ts.alt[i1] <- (1/sqrt(2*m*i1*sigma^2))*(diff(colSums(obs.alt[1:(i1*m),])))
      }
    }
    if(method=="rerank"){
      mu <- (1:k)*m*((1:k)*(m+m)+1)/2
      sig <- sqrt((1:k)^2*m*m*((1:k)*(m+m)+1)/12)
      ts.alt <- (ts.alt-mu)/sig
    } else if(method=="sar"){
      mu <- m*((m+m)+1)/2
      sig <- sqrt(m*m*((m+m)+1)/12)/sqrt(1:k)
      ts.alt <- cumsum(ts.alt)/seq_along(ts.alt)
      ts.alt <- (ts.alt-mu)/sig
    }
    
    if(is.na(which(ts.alt[1:(k-1)]>=r[1:(k-1)])[1]) & is.na(which(ts.alt[1:(k-1)]<=a[1:(k-1)])[1])){
      term <- k
      if(ts.alt[k]>=r[k]){out.r <- append(out.r,k)}
      if(ts.alt[k]<=a[k]){out.a <- append(out.a,k)}
    } else if(is.na(which(ts.alt[1:(k-1)]>=r[1:(k-1)])[1])){
      term <- which(ts.alt[1:(k-1)]<=a[1:(k-1)])[1]
      out.a <- append(out.a,term)
    } else {
      term <- which(ts.alt[1:(k-1)]>=r[1:(k-1)])[1]
      out.r <- append(out.r,term)
    }
    
    stage[i] <- term
    if(eps==0){eps1 <- var(obs.alt[1:(term*m),1])*log((2*term*m)^2)/(2*term*m)}
    iqr.adj <- qpsi(0.75,mu=delta,sig=1)-qpsi(0.25,mu=delta,sig=1)
    if(eps==0){
      sx <- (IQR(obs.alt[1:(term*m),1])/(iqr.adj))
      eps.adj <- sx^2*log((2*term*m)^2)/(2*term*m)
    }
    means <- colMeans(obs.alt[1:(term*m),])
    vars <- apply(obs.alt[1:(term*m),],2,var)
    meds <- apply(obs.alt[1:(term*m),],2,median)
    iqrs <- apply(obs.alt[1:(term*m),],2,IQR)
    km.res <- kmeans(obs.alt[1:(term*m),2], centers=2, nstart=50)
    theta.hat.km[i] <- km.res$size[which(km.res$centers==max(km.res$centers))]/sum(km.res$size)
    delta.hat.km[i] <- max(km.res$centers)-min(km.res$centers)
    
    sy <- sqrt(sx^2 + theta.hat.km[i]*(1-theta.hat.km[i])*(delta.hat.km[i])^2)
    iqrs.adj <- c(sx,sy)
    temp <- 1 + max(diff(vars),0)/(max(diff(means),0)^2 + eps1)
    temp.adj <- 1 + max(diff(iqrs.adj^2),0)/(max(diff(meds),0)^2 + eps.adj)
    # mom uncorrected ests
    theta.hat[i] <- 1/temp
    delta.hat[i] <- max(diff(means),0)*temp
    if(diff(means)<0){
      theta.hat.check[i]<-delta.hat.check[i]<-0
      theta.hat.deps[i]<-delta.hat.deps[i]<-0
      theta.hat.noeps[i]<-delta.hat.noeps[i]<-0
    }else{
      temp <- 1 + max(diff(vars),0)/(diff(means)^2 + eps1)
      theta.hat.check[i] <- 1/temp
      delta.hat.check[i] <- diff(means)*temp
      temp <- 1 + max(diff(vars),0)/(diff(means)^2)
      theta.hat.deps[i] <- 1/temp
      delta.hat.deps[i] <- diff(means)*(temp+eps1)
      theta.hat.noeps[i] <- 1/temp
      delta.hat.noeps[i] <- diff(means)*temp
    }
    # mom ests using adjusted iqr instead of sample variance
    if(diff(meds)<0){theta.hat.adj[i]<-delta.hat.adj[i]<-0}else{
      theta.hat.adj[i] <- 1/temp.adj
      delta.hat.adj[i] <- diff(meds)*temp.adj
    }
    
    min.ind <- which.min(c(obs.alt[1:(term*m),]))
    max.ind <- (term*m)+which.max(obs.alt[1:(term*m),2])
    cantLin<-matrix(c(min.ind,max.ind),nrow=1)
    
    term.obs <- obs.alt[1:(term*m),]
    term.x <- term.obs[,1]; term.y <- term.obs[,2]
    
    sind <- sample(1:(2*term*m));ckm.ind <- sind
    ckm <- conclust::ckmeans(data=c(term.obs)[sind], k=2, mustLink=t(combn(which(sind %in% 1:(term*m)),2)), cantLink = matrix(c(which.min(c(term.obs)[sind]),which(max(term.y)==c(term.obs)[sind])),nrow=1))
    ckm.ss <- sum(aggregate(c(term.obs),by=list(ckm[order(sind)]),FUN=function(x){sum((x-mean(x))^2)/length(x)})$x)
    for(i1 in 1:49){
      sind <- sample(1:(2*term*m))
      temp.ckm <- conclust::ckmeans(data=c(term.obs)[sind], k=2, mustLink=t(combn(which(sind %in% 1:(term*m)),2)), cantLink = matrix(c(which.min(c(term.obs)[sind]),which(max(term.y)==c(term.obs)[sind])),nrow=1))
      ss <- sum(aggregate(c(term.obs),by=list(temp.ckm[order(sind)]),FUN=function(x){sum((x-mean(x))^2)/length(x)})$x)
      if(ss<ckm.ss){
        ckm <- temp.ckm; ckm.ss <- ss; ckm.ind <- sind
      }
      rm(i1);rm(temp.ckm);rm(ss);rm(sind)
    }
    tckm.ind <- ifelse(ckm[order(ckm.ind)][1]==1,2,1) # which cluster indicator is the treatment group
    theta.hat.ckm[i] <- sum(ckm[order(ckm.ind)][(term*m+1):(2*term*m)]==tckm.ind)/(term*m)
    delta.hat.ckm[i] <- mean(c(term.obs)[ckm[order(ckm.ind)]==tckm.ind])-mean(c(term.obs)[ckm[order(ckm.ind)]!=tckm.ind])
    
    theta.hat.ckmx[i] <- delta.hat.ckmx[i] <- 0
    
    # mle
    comb.log.lik <- function(p,x,y){
      # add small value 1e-10 to avoid optimization problems
      -(sum(des.dpsi(x,mu=p[1],sig=p[2],lg=T)) + sum(log((1-p[4])*des.dpsi(y,mu=p[1],sig=p[2]) + p[4]*des.dpsi(y,mu=p[1]+p[3],sig=p[2]) + 1e-10)))
    }
    
    if(diff(means)<0){
      bc.theta.hat.np[i] <- 0
      bc.delta.hat.np[i] <- 0
    } else {
      if(bootstrap>0){ # bootstrap ests
        t.hat <- 1/temp
        d.hat <- max(diff(means),0)*temp
        t.hat.boot.p <- d.hat.boot.p <- vector("numeric",length=bootstrap)
        t.hat.boot.np <- d.hat.boot.np <- vector("numeric",length=bootstrap)
        t.hat.boot.mle <- d.hat.boot.mle <- vector("numeric",length=bootstrap)
        for(j in 1:bootstrap){
          
          # parametric bootstrap
          ind.boot <- sample(1:2,prob=c(1-t.hat,t.hat),size=k*m,replace = TRUE)
          mus.boot <- c(0,d.hat)
          y.boot <- f(n=k*m,mu=mus.boot[ind.boot],sig=sqrt(vars)[1])
          obs.alt.boot.p <- matrix(c(f(n=k*m,sig=sqrt(vars)[1]),y.boot),ncol=2,nrow=(k*m),byrow=F)
          # nonparametric bootstrap with ogive and exp tails to avoid ties
          x.vals <- obs.alt[1:(term*m),1]
          fhat <- sapply(x.vals,function(x){sum(x.vals<=x)/m})
          probs <- runif(m)
          boot.samp <- ifelse(probs >= min(fhat) & probs <= sort(fhat)[m-1],
                              approx(x=sort(unique(fhat)), y=sort(unique(x.vals)), xout=probs)$y,
                              ifelse(r < min(fhat), -log(probs)/(-log(min(fhat))/min(x.vals)),
                                     -log(1-probs)/(-log(1-sort(fhat)[m-1])/sort(x.vals)[m-1])))
          obs.alt.boot.np <- matrix(boot.samp+c(rep(0,k*m),mus.boot[ind]),ncol=2,nrow=(k*m),byrow=F)
          
          ts.alt.boot.p <- ts.alt.boot.np <- vector("numeric",length=k)
          for(i1 in 1:k){
            if(method=="rerank"){
              if(ties=="average"||ties=="min"||ties=="max"||ties=="random"){
                # standardize with variance for ties?
                r <- rank(obs.alt.boot[1:(i1*m),],ties.method = ties)[(i1*m+1):(i1*m+i1*m)]
                ts.alt.boot[i1] <- sum(r)
                # getting variance with ties
                t.r <- table(r)
                sig[i1] <- sqrt(((1:i1)*m)^2*(2*(1:i1)*m + 1)/12 - (((1:i1)*m)^2/(12*2*(1:i1)*m*(2*(1:i1)*m-1))*sum((t.r-1)*t.r*(t.r+1))))
              } else if(ties=="jitter"){ # eliminated ties, this is what's used
                # for jittered data
                ts.alt.boot[i1] <- sum(rank(obs.alt.boot[1:(i1*m),])[(i1*m+1):(i1*m+i1*m)])
              }
            } else if(method=="sar"){
              
              # not dealing with ties anymore
              ts.alt.boot.p[i1] <- sum(rank(obs.alt.boot.p[((i1-1)*m+1):(i1*m),])[(n+1):(n+m)])
              ts.alt.boot.np[i1] <- sum(rank(obs.alt.boot.np[((i1-1)*m+1):(i1*m),])[(n+1):(n+m)])
              
            } else if(method=="average"){
              ts.alt.boot[i1] <- (1/sqrt(2*m*i1*sigma^2))*(diff(colSums(obs.alt[1:(i1*m),])))
            }
          }
          if(method=="rerank"){
            mu <- (1:k)*m*((1:k)*(m+m)+1)/2
            sig <- sqrt((1:k)^2*m*m*((1:k)*(m+m)+1)/12)
            ts.alt.boot <- (ts.alt.boot-mu)/sig
          } else if(method=="sar"){
            mu <- m*((m+m)+1)/2
            sig <- sqrt(m*m*((m+m)+1)/12)/sqrt(1:k)
            ts.alt.boot.p <- cumsum(ts.alt.boot.p)/seq_along(ts.alt.boot.p)
            ts.alt.boot.np <- cumsum(ts.alt.boot.np)/seq_along(ts.alt.boot.np)
            ts.alt.boot.p <- (ts.alt.boot.p-mu)/sig
            ts.alt.boot.np <- (ts.alt.boot.np-mu)/sig
          }
          # for parametric
          if(is.na(which(ts.alt.boot.p[1:(k-1)]>=r[1:(k-1)])[1]) & is.na(which(ts.alt.boot.p[1:(k-1)]<=a[1:(k-1)])[1])){
            term.boot.p <- k
          } else if(is.na(which(ts.alt.boot.p[1:(k-1)]>=r[1:(k-1)])[1])){
            term.boot.p <- which(ts.alt.boot.p[1:(k-1)]<=a[1:(k-1)])[1]
          } else {
            term.boot.p <- which(ts.alt.boot.p[1:(k-1)]>=r[1:(k-1)])[1]
          }
          # for nonparametric
          if(is.na(which(ts.alt.boot.np[1:(k-1)]>=r[1:(k-1)])[1]) & is.na(which(ts.alt.boot.np[1:(k-1)]<=a[1:(k-1)])[1])){
            term.boot.np <- k
          } else if(is.na(which(ts.alt.boot.np[1:(k-1)]>=r[1:(k-1)])[1])){
            term.boot.np <- which(ts.alt.boot.np[1:(k-1)]<=a[1:(k-1)])[1]
          } else {
            term.boot.np <- which(ts.alt.boot.np[1:(k-1)]>=r[1:(k-1)])[1]
          }
          # parametric bootstrap ests
          if(eps==0){eps.p <- var(obs.alt.boot.p[1:(term.boot.p*m),1])*log((2*term.boot.p*m)^2)/(2*term.boot.p*m)}
          means.boot.p <- colMeans(obs.alt.boot.p[1:(term.boot.p*m),])
          vars.boot.p <- apply(obs.alt.boot.p[1:(term.boot.p*m),],2,var)
          temp.boot.p <- 1 + max(diff(vars.boot.p),0)/(max(diff(means.boot.p),0)^2 + eps.p)
          t.hat.boot.p[j] <- 1/temp.boot.p
          d.hat.boot.p[j] <- max(diff(means.boot.p),0)*temp.boot.p
          # mle
          mle <- optim(c(0,1,0.5,0.5), comb.log.lik, x=obs.alt.boot.p[1:(term.boot.p*m),1], y=obs.alt.boot.p[1:(term.boot.p*m),2],
                       method = "L-BFGS-B",lower=c(-Inf, 1e-10, -Inf, 1e-10), upper = c(Inf,Inf,Inf,1-1e-10))$par
          t.hat.boot.mle[j] <- mle[4]
          d.hat.boot.mle[j] <- mle[3]
          
          # nonparametric bootstrap ests
          if(eps==0){eps.np <- var(obs.alt.boot.np[1:(term.boot.np*m),1])*log((2*term.boot.np*m)^2)/(2*term.boot.np*m)}
          means.boot.np <- colMeans(obs.alt.boot.np[1:(term.boot.np*m),])
          vars.boot.np <- apply(obs.alt.boot.np[1:(term.boot.np*m),],2,var)
          temp.boot.np <- 1 + max(diff(vars.boot.np),0)/(max(diff(means.boot.np),0)^2 + eps.np)
          t.hat.boot.np[j] <- 1/temp.boot.np
          d.hat.boot.np[j] <- max(diff(means.boot.np),0)*temp.boot.np
        }
        bc.theta.hat.p[i] <- 2*t.hat - mean(t.hat.boot.p)
        bc.delta.hat.p[i] <- 2*d.hat - mean(d.hat.boot.p)
        bc.theta.hat.mle[i] <- 2*t.hat - mean(t.hat.boot.mle)
        bc.delta.hat.mle[i] <- 2*d.hat - mean(d.hat.boot.mle)
        bc.theta.hat.np[i] <- 2*t.hat - mean(t.hat.boot.np)
        bc.delta.hat.np[i] <- 2*d.hat - mean(d.hat.boot.np)
      }
    }
  } # not the end
  
  theta.hat.ckmch <- ifelse(delta.hat.ckm<0,0,theta.hat.ckm)
  delta.hat.ckmch <- ifelse(delta.hat.ckm<0,0,delta.hat.ckm)
  theta.hat.ckmxch <- ifelse(delta.hat.ckmx<0,0,theta.hat.ckmx)
  delta.hat.ckmxch <- ifelse(delta.hat.ckmx<0,0,delta.hat.ckmx)
  
  theta.hat <- theta.hat.check
  delta.hat <- delta.hat.check
  t.bias <- mean(theta.hat) - theta
  d.bias <- mean(delta.hat) - delta
  t.mse <- var(theta.hat) + t.bias^2
  d.mse <- var(delta.hat) + d.bias^2
  t.out <- c(t.bias,t.mse)
  d.out <- c(d.bias,d.mse)
  names(t.out) <- names(d.out) <- c("bias","mse")
  cat("theta =",theta,"\n")
  cat("mean =",mean(theta.hat),"\n")
  cat("bias =",t.bias,"\n")
  
  # mom with adjusted iqrs
  t.bias.adj <- mean(theta.hat.adj) - theta
  d.bias.adj <- mean(delta.hat.adj) - delta
  t.mse.adj <- var(theta.hat.adj) + t.bias.adj^2
  d.mse.adj <- var(delta.hat.adj) + d.bias.adj^2
  
  # kmeans ests
  t.bias.km <- mean(theta.hat.km) - theta
  d.bias.km <- mean(delta.hat.km) - delta
  t.mse.km <- var(theta.hat.km) + t.bias.km^2
  d.mse.km <- var(delta.hat.km) + d.bias.km^2
  
  t.bias.ckmch <- mean(theta.hat.ckmch) - theta
  d.bias.ckmch <- mean(delta.hat.ckmch) - delta
  t.mse.ckmch <- var(theta.hat.ckmch) + t.bias.ckmch^2
  d.mse.ckmch <- var(delta.hat.ckmch) + d.bias.ckmch^2
  
  t.bias.ckmxch <- mean(theta.hat.ckmxch) - theta
  d.bias.ckmxch <- mean(delta.hat.ckmxch) - delta
  t.mse.ckmxch <- var(theta.hat.ckmxch) + t.bias.ckmxch^2
  d.mse.ckmxch <- var(delta.hat.ckmxch) + d.bias.ckmxch^2
  
  agg.fun.t <- function(x) c(bias = mean(x) - theta, root.mse = sqrt(var(x) + (mean(x)-theta)^2))
  agg.fun.d <- function(x) c(bias = mean(x) - delta, root.mse = sqrt(var(x) + (mean(x)-delta)^2))
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
  # comments are if going to do breakdown of stages for bootstraps
  # agg.fun <- function(x) c(n = length(x), bias = mean(x) - theta, root.mse = sqrt(var(x) + (mean(x)-theta)^2))
  t.df <- specify_decimal(rbind(aggregate(theta.hat,by=list(stage),FUN = agg.fun.t)[,-1],c(t.bias,sqrt(t.mse))),4)
  d.df <- specify_decimal(rbind(aggregate(delta.hat,by=list(stage),FUN = agg.fun.d)[,-1],c(d.bias,sqrt(d.mse))),4)
  # t.df <- round(rbind(aggregate(theta.hat,by=list(stage),FUN = agg.fun)[,-1],c(1000,t.bias,sqrt(t.mse))),4)
  # d.df <- round(rbind(aggregate(delta.hat,by=list(stage),FUN = agg.fun)[,-1],c(1000,d.bias,sqrt(d.mse))),4)
  # compile mom ests with adjusted iqr for output
  t.df.adj <- specify_decimal(rbind(aggregate(theta.hat.adj,by=list(stage),FUN = agg.fun.t)[,-1],c(t.bias.adj,sqrt(t.mse.adj))),4)
  d.df.adj <- specify_decimal(rbind(aggregate(delta.hat.adj,by=list(stage),FUN = agg.fun.d)[,-1],c(d.bias.adj,sqrt(d.mse.adj))),4)
  
  # compile km ests for output
  t.df.km <- specify_decimal(rbind(aggregate(theta.hat.km,by=list(stage),FUN = agg.fun.t)[,-1],c(t.bias.km,sqrt(t.mse.km))),4)
  d.df.km <- specify_decimal(rbind(aggregate(delta.hat.km,by=list(stage),FUN = agg.fun.d)[,-1],c(d.bias.km,sqrt(d.mse.km))),4)
  # compile kmeans ests for output
  t.df.ckmch <- specify_decimal(rbind(aggregate(theta.hat.ckmch,by=list(stage),FUN = agg.fun.t)[,-1],c(t.bias.ckmch,sqrt(t.mse.ckmch))),4)
  d.df.ckmch <- specify_decimal(rbind(aggregate(delta.hat.ckmch,by=list(stage),FUN = agg.fun.d)[,-1],c(d.bias.ckmch,sqrt(d.mse.ckmch))),4)
  t.df.ckmxch <- specify_decimal(rbind(aggregate(theta.hat.ckmxch,by=list(stage),FUN = agg.fun.t)[,-1],c(t.bias.ckmxch,sqrt(t.mse.ckmxch))),4)
  d.df.ckmxch <- specify_decimal(rbind(aggregate(delta.hat.ckmxch,by=list(stage),FUN = agg.fun.d)[,-1],c(d.bias.ckmxch,sqrt(d.mse.ckmxch))),4)
  
  # num.ests <- 6 # mle, mom, mom.adj, mom.mad, mad.adj, km
  # num.ests <- 4 # mle, mom, mom.adj, km
  num.ests <- 3
  temp.out <- data.frame(theta=rep(theta,1+2*num.ests),delta=rep(delta,1+2*num.ests),
                         est=c("n",rep(c("MoM","kmeans","ckmnch"),each=2)),
                         type=c("",rep(c("bias","rmse"),num.ests)))
  t.out <- cbind(temp.out,rbind(c(tabulate(stage),B),t(t.df),t(t.df.km),t(t.df.ckmch)))
  d.out <- cbind(temp.out,rbind(c(tabulate(stage),B),t(d.df),t(d.df.km),t(d.df.ckmch)))
  colnames(t.out)[ncol(t.out)-k:1] <- colnames(d.out)[ncol(d.out)-k:1] <- paste("Stage",1:k)
  colnames(t.out)[ncol(t.out)] <- colnames(d.out)[ncol(d.out)] <- "overall"
  
  if(bootstrap>0){
    bc.t.bias.p <- mean(bc.theta.hat.p) - theta
    bc.d.bias.p <- mean(bc.delta.hat.p) - delta
    bc.t.mse.p <- var(bc.theta.hat.p) + bc.t.bias.p^2
    bc.d.mse.p <- var(bc.delta.hat.p) + bc.d.bias.p^2
    bc.t.bias.np <- mean(bc.theta.hat.np) - theta
    bc.d.bias.np <- mean(bc.delta.hat.np) - delta
    bc.t.mse.np <- var(bc.theta.hat.np) + bc.t.bias.np^2
    bc.d.mse.np <- var(bc.delta.hat.np) + bc.d.bias.np^2
    
    bc.t.bias.mle <- mean(bc.theta.hat.mle) - theta
    bc.d.bias.mle <- mean(bc.delta.hat.mle) - delta
    bc.t.mse.mle <- var(bc.theta.hat.mle) + bc.t.bias.mle^2
    bc.d.mse.mle <- var(bc.delta.hat.mle) + bc.d.bias.mle^2
    
    t.df.bc.p <- specify_decimal(rbind(aggregate(bc.theta.hat.p,by=list(stage),FUN = agg.fun.t)[,-1],c(bc.t.bias.p,sqrt(bc.t.mse.p))),4)
    d.df.bc.p <- specify_decimal(rbind(aggregate(bc.delta.hat.p,by=list(stage),FUN = agg.fun.d)[,-1],c(bc.d.bias.p,sqrt(bc.d.mse.p))),4)
    
    t.df.bc.np <- specify_decimal(rbind(aggregate(bc.theta.hat.np,by=list(stage),FUN = agg.fun.t)[,-1],c(bc.t.bias.np,sqrt(bc.t.mse.np))),4)
    d.df.bc.np <- specify_decimal(rbind(aggregate(bc.delta.hat.np,by=list(stage),FUN = agg.fun.d)[,-1],c(bc.d.bias.np,sqrt(bc.d.mse.np))),4)
    
    t.df.bc.mle <- specify_decimal(rbind(aggregate(bc.theta.hat.mle,by=list(stage),FUN = agg.fun.t)[,-1],c(bc.t.bias.mle,sqrt(bc.t.mse.mle))),4)
    d.df.bc.mle <- specify_decimal(rbind(aggregate(bc.delta.hat.mle,by=list(stage),FUN = agg.fun.d)[,-1],c(bc.d.bias.mle,sqrt(bc.d.mse.mle))),4)
    
    t.mat.h <- rbind(t.mat.h, matrix(c(paste(" & & bc.np &",paste(t(t.df.bc.np),collapse = " & "),"\\"),
                                       paste(" & & bc.p &",paste(t(t.df.bc.p),collapse = " & "),"\\"),
                                       paste(" & & bc.mle &",paste(t(t.df.bc.mle),collapse = " & "),"\\")
    ),
    ncol=1))
    t.mat.v <- rbind(t.mat.v, matrix(c(paste(" & & bc.np &",paste(t.df.bc.np[,1],collapse = " & "),"\\"),
                                       paste0(" & & & (",paste(t.df.bc.np[,2],collapse = ") & ("),") \\"),
                                       paste(" & & bc.p &",paste(t.df.bc.p[,1],collapse = " & "),"\\"),
                                       paste0(" & & & (",paste(t.df.bc.p[,2],collapse = ") & ("),") \\"),
                                       paste(" & & bc.mle &",paste(t.df.bc.mle[,1],collapse = " & "),"\\"),
                                       paste0(" & & & (",paste(t.df.bc.mle[,2],collapse = ") & ("),") \\")
    ),
    ncol=1))
    d.mat.h <- rbind(d.mat.h, matrix(c(paste(" & & bc.np &",paste(t(d.df.bc.np),collapse = " & "),"\\"),
                                       paste(" & & bc.p &",paste(t(d.df.bc.p),collapse = " & "),"\\"),
                                       paste(" & & bc.mle &",paste(t(d.df.bc.mle),collapse = " & "),"\\")
    ),
    ncol=1))
    d.mat.v <- rbind(d.mat.v, matrix(c(paste(" & & bc.np &",paste(d.df.bc.np[,1],collapse = " & "),"\\"),
                                       paste0(" & & & (",paste(d.df.bc.np[,2],collapse = ") & ("),") \\"),
                                       paste(" & & bc.p &",paste(d.df.bc.p[,1],collapse = " & "),"\\"),
                                       paste0(" & & & (",paste(d.df.bc.p[,2],collapse = ") & ("),") \\"),
                                       paste(" & & bc.mle &",paste(d.df.bc.mle[,1],collapse = " & "),"\\"),
                                       paste0(" & & & (",paste(d.df.bc.mle[,2],collapse = ") & ("),") \\")
    ),
    ncol=1))
    
    if(out.tab=="t"){return(list(wide=format(data.frame(theta.tab=t.mat.h),justify="left"),
                                 long=format(data.frame(theta.tab=t.mat.v),justify="left")))
    } else if(out.tab=="d"){return(list(wide=format(data.frame(delta.tab=d.mat.h),justify="left"),
                                        long=format(data.frame(delta.tab=d.mat.v),justify="left")))
    } else{
      return(format(data.frame(tab=rbind(t.mat.h,paste("&  &  &  &  \\"),d.mat.h)),justify="left"))
    }
    
  } else {
    
    return(list(theta.tab=t.out,delta.tab=d.out))
  }
}



# robust estimation
# run a number of different estimation simulations
for(o1 in 1){
  if(o1==1) est.des <- list(k=2,theta=0.7,K=1.5,f="Normal") # m=17
  if(o1==2) est.des <- list(k=3,theta=0.6,K=0.5,f="Laplace") # m=37
  print(as.data.frame(est.des))
  x<-grp.seq.sim(k=est.des$k,alpha=0.01,beta=0.1,rho=2,theta=est.des$theta,K=est.des$K,sigma=1,
                 f=tolower(est.des$f),B=100000,method="sar",std=TRUE,norm.approx=TRUE)
  for(i1 in c("Normal","Logistic","Laplace","t")){
    for(i2 in est.des$theta + c(0)){
      for(i3 in est.des$K + c(1.5)){
        res <- est.sim(k=est.des$k,sigma=1,theta=i2,K=i3,f=tolower(i1),method="sar",des.f=tolower(est.des$f),m=x$m,r=x$r,a=x$a,
                       eps=0.00,bootstrap=0e3,nonpar=FALSE,B=1e3,ties="jitter",out.tab="b")
        
        save(res,file=paste0("C:/Users/Dylan/Documents/Prof Jeske/code/Cluster Estimation/Des-k",est.des$k,"-",est.des$f,"-Th",est.des$theta,"-K",est.des$K,"/",i1,"-theta-",i2,"-delta-",i3,".Rda"))
        cat("F:",i1,"\n")
        cat("theta:",i2,"\n")
        cat("delta:",i3,"\n")
        rm(i3);rm(res)
      }
    }
  }
  rm(x);rm(i1);rm(i2);rm(o1)
}


# function for help formatting tables for latex
empty_cells <- function(x){
  value <- x[1]
  for(i in 2:length(x)){
    if(x[i]==value){ x[i] <- ""
    } else { value <- x[i]}
  }
  return(x)
}


# compile tables
# read simulation results and create table for latex
for(o1 in 1){
  metric <- "Bias"
  # metric <- "RMSE"
  # if(o1==1) est.des <- list(k=2,theta=0.8,K=0.75,f="Normal",m=21) # m=21
  if(o1==1) est.des <- list(k=2,theta=0.7,K=1.5,f="Normal") # m=17
  # if(o1==2) est.des <- list(k=3,theta=0.6,K=0.5,f="Laplace",m=37) # m=37
  main.dir <- paste0("C:/Users/Dylan/Documents/Prof Jeske/code/Cluster Estimation/Des-k",est.des$k,"-",est.des$f,"-Th",est.des$theta,"-K",est.des$K,"/")
  for(o2 in c("theta","delta")){
    if(o2=="theta"){ind <- 9:13}
    if(o2=="delta"){ind <- 18:22}
    comp.table <- NULL
    for(i1 in c("Normal","Logistic","Laplace","t")){
      full.table <- NULL
      for(i2 in est.des$theta + c(-0.7)){
        for(i3 in est.des$K + c(-1.5)){
          # txt file
          # path.file <- paste0(i1,"-theta-",i2,"-delta-",i3,".txt")
          # full.table <- rbind(full.table,data.frame(x=readLines(paste0(main.dir,i1,"-theta-",i2,"-delta-",i3,".txt"))[ind]))
          # rda file
          path.file <- paste0(i1,"-theta-",i2,"-delta-",i3,".Rda")
          load(paste0(main.dir,path.file))
          full.table <- rbind(full.table,res[[paste0(o2,".tab")]])
        }
      }
      full.table <- full.table[with(full.table, order(theta,delta,factor(est,levels=c("n","MoM","kmeans","ckmnch")))),]
      
      if(is.null(comp.table)){ # start comparison table if doesn't exist
        comp.table <- full.table[full.table$type==tolower(metric),][c("theta","delta","est","type","overall")]
      } else { comp.table <- cbind(comp.table,full.table[full.table$est!="n" & full.table$type==tolower(metric),]["overall"]) }
      colnames(comp.table)[ncol(comp.table)] <- i1
      
      full.table <- full.table[full.table$type==tolower(metric),-which(colnames(full.table)=="type")]
      full.table$theta <- empty_cells(full.table$theta)
      full.table$delta <- empty_cells(full.table$delta)
      names(full.table)[1:3] <- c("reptheta","repdelta","Estimator")
      names(full.table)[ncol(full.table)] <- "Overall"
      full.lab <- paste0("tab:des",o1,"_",o2,"_",tolower(metric),"_",tolower(i1))
      if(i1=="t") i1 <- "$t_3$"
      full.cap <- paste0(metric," values for estimating $\\",o2,"$ where the data come from the ",i1,
                         " distribution with $\\theta$ and $\\delta$ values indicated in the table for a design scenario with ",
                         paste(c("$S","\\theta","\\delta","F","\\alpha","\\beta"),"=$",c(est.des[-5],0.05,0.2),collapse = ", $"),
                         ", and $\\rho = 2$ resulting in an arm size of ",est.des$m,".")
    }
    comp.table$theta <- empty_cells(comp.table$theta)
    comp.table$delta <- empty_cells(comp.table$delta)
    comp.cap <- paste0(metric," values for estimating $\\",o2,"$ where the data come from the distribution, $\\theta$, and $\\delta$ values indicated in the table for a design scenario with ",
                       paste(c("$S","\\theta","\\delta","F","\\alpha","\\beta"),"=$",c(est.des[-5],0.05,0.2),collapse = ", $"),
                       ", and $\\rho = 2$ resulting in an arm size of ",est.des$m,".")
    names(comp.table)[1:3] <- c("reptheta","repdelta","Estimator")
    names(comp.table)[which(names(comp.table)=="t")] <- "rept3"
    comp.lab <- paste0("tab:des",o1,"_",o2,"_",tolower(metric),"_comp")
    comp.table <- comp.table[,-which(colnames(comp.table)=="type")] # remove "type" column
    print(xtable::xtable(comp.table,caption = comp.cap,label=comp.lab,align=c("l",rep("l",3),rep("r",ncol(comp.table)-3))),
          include.rownames = FALSE)
    agg.table <- comp.table
    agg.table$Estimator <- factor(agg.table$Estimator, levels = c("MLE","PAM","MoM","BC MoM","Robust MoM"))
    agg.table[,4:7] <- apply(agg.table[,4:7],2,as.numeric)
    agg.cap <- paste0("Average ", metric," values for estimating $\\",o2,"$ where the data come from the distribution indicated in the table and averaged over each combination of $\\theta = (",
                      paste0(est.des$theta + c(-0.1,0,0.1),collapse=", "),")$ and $\\delta = (",paste0(est.des$K + c(-0.25,0,0.25),collapse=", "),")$ for the design scenario with ",
                      paste(c("$S","\\theta","\\delta","F","\\alpha","\\beta"),"=$",c(est.des[-5],0.05,0.2),collapse = ", $"),
                      ", and $\\rho = 2$ resulting in an arm size of ",est.des$m,".")
  }
  rm(o1);rm(o2);rm(est.des);rm(main.dir);rm(path.file);rm(i1);rm(i2);rm(i3);rm(res);rm(metric);rm(full.cap);rm(comp.cap)
  rm(agg.table)
}
# then replace all tabular:smalltabular, reptheta:$\theta$, repdelta:$\delta$, rept3:$t_3$


