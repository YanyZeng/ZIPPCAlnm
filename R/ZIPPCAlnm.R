#' @title ZIPPCAlnm
#' @description Zero-Inflated Probabilistic PCA framework with logistical normal multinomial distribution (ZIPPCA-LNM),
#' that extends probabilistic PCA from the Gaussian setting to multivariate abundance data, and an empirical Bayes approach was
#' proposed for inferring microbial compositions. An efficient VA algorithm, classification VA, has been developed for fitting this model.
#' @param X count matrix of observations.
#' @param V vector of sample covariate.
#' @param d_choice FALSE, “BIC or "CV". Indicating whether the rank or number of factors, or dimension after dimensional reduction,
#'                 will be chosen from 1 to 5.  Options are "BIC" (Bayesian information criterion), and "CV" (Cross-validation).Defaults to FALSE.
#' @param parallel logical, if TRUE, use parallel toolbox to accelerate.

#' @return
#'
#'
#'  \item{VLB }{ variational lower bound of log likelihood}
#'  \item{lvs}{list of latent variables
#'  \itemize{
#'    \item{pi }{ the probabilities of excess zeros}
#'    \item{factor_scores }{ coordinates or factor scores in low-dimensional subspace}
#'    }}
#'  \item{params}{list of model parameters
#'  \itemize{
#'    \item{factor_coefs_j }{ coefficients of latent variables fator scores or factor loadings}
#'    \item{factor_coefs_0 }{ taxon-specific intercepts}
#'    \item{gamma }{ coeffcients of sample covariate}
#'    }}
#'  \item{Q }{ the underlying composition of microbiome data}
#'  \item{bic}{ if d_choice is "BIC", the number of the rank or factors, or dimension, will be chosen by BIC type information criterion}
#'  \item{cv}{ if d_choice is "CV", the number of the rank or factors, or dimension, will be chosen by Cross-validation}

#' @examples
#' n.n = 50
#' n.w = 100
#' n.factors = 2
#' set.seed(1)
#' si <- diag(n.factors)
#' me <- c(0,0)
#' f <- matrix(0,nrow = n.n, ncol = n.factors)
#' for(i in 1:n.n){
#'  f[i,] <- mvrnorm(1,me,si)
#' }

#' betaj <- matrix(0,nrow = n.w, ncol = n.factors)
#' for(j in 1:n.w){
#'   betaj[j,] <- runif(n.factors,-3.5,3.5)
#' }
#' beta0 <- rep(2,n.w)
#' l <- matrix(beta0,n.n,n.w,byrow=TRUE)+f %*%  t(betaj)
#' eta_j <- rep(0.25,n.w)
#' z <- matrix(0,n.n,n.w,byrow = TRUE)
#' for(i in 1:n.n){
#'   z[i,] <- rbinom(n.w, size=1, prob=eta_j)
#' }
#' sum <- rowSums(exp(l))
#' Qn <- exp(l)/sum
#' sum <- matrix(rowSums((1-z)*exp(l)),n.n,n.w)
#' Qn_z <- matrix(I(z==0)*(exp(l)/sum)+I(z==1)*0,n.n,n.w)
#' X <- matrix(0,n.n,n.w,byrow = TRUE)
#' for(i in 1:n.n){
#'   X[i,] <- rmultinom(1, size = runif(1,800,1000), prob = Qn_z[i,])
#' }
#' zerorow <- which(rowSums(X)==0)
#' if(length(zerorow) >0 ){
#'   X <- X[-zerorow,];Qn <- Qn[-zerorow,];Qn_z <- Qn_z[-zerorow,];
#' }
#' zerocol <- which(colSums(X)==0)
#' if(length(zerocol) >0 ){
#'   X <- X[,-zerocol];Qn <- Qn[,-zerocol];Qn_z <- Qn_z[,-zerocol];
#' }
#' result <- ZIPPCAlnm::ZIPPCAlnm(X,V=NULL,d_choice=FALSE,parallel=TRUE)

#' @export

ZIPPCAlnm <- function(X,V=NULL,d_choice=FALSE,parallel=TRUE){

  ZILNMVA <- function(X,V,trace = FALSE,n.factors=2,maxit = 100,cv_group=NULL) {

    n.s<-nrow(X); n.f<-ncol(X);
    if(is.null(V)){Y <- 0
    }else if(is.numeric(V)){Y <- V
    }else{
      Y <- as.numeric(as.factor(V))-1}
    M <- rowSums(X)
    if (is.null(cv_group)) {
      cvsample <- matrix(0,n.s,n.f)
    }else{
      cvsample <- cv_group
    }
    out.list <- list()

    ### Initialization 1
    pzero.col <- apply(X, 2, function(x) {sum(x==0)/n.s})
    z.hat <- pi <- new.pi <- t(ifelse(t(X)==0, pzero.col, 0))
    eta <- new.eta <- round(apply(pi, 2, mean), 6)
    sigma <- new.sigma <- matrix(1,n.s,n.factors)
    factor_coefs_0 <-  new.factor_coefs_0 <-  rep(1,n.f)
    gamma <-  new.gamma <-  rep(1,n.f)

    X.rc <- scale(log(X+0.05),scale = T,center = T)
    re <- svd(X.rc,n.factors,n.factors)
    factor_coefs_j <- new.factor_coefs_j <- re$v
    if(n.factors==1){
      factor_scores <- new.factor_scores <- re$u * (re$d[1])
    }else{factor_scores <- new.factor_scores <- re$u %*% diag(re$d[1:n.factors])}

    ### Initialization 2
    cur.VLB <- -1e6; iter <- 1; ratio <- 10; diff=1e5;eps = 1e-4;max.iter = 100;
    b.cur.logfunc <- -1e6;b0.cur.logfunc <- -1e6;f.cur.logfunc <- -1e6;ga.cur.logfunc <- -1e6
    while((diff> eps*(abs(cur.VLB)+eps)) && iter <= max.iter) {
      if(trace) cat("Iteration:", iter, "\n")
      ## VLB
      VLB_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*ll*I(cvsample==0)
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))}
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))
        return(y)
      }

      ###optim f
      f_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lf*I(cvsample==0)
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        fun2 <- function(i) { 0.5 * (- sum(new.factor_scores[i,]^2))}
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))
        return(y)
      }
      f_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0)/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0)))
        fun2 <- function(i) { -new.factor_scores[i,]+((X*I(cvsample==0))[i,])%*%new.factor_coefs_j-(M[i]*sum[i,])%*%new.factor_coefs_j}
        f_grad <- t(sapply(1:n.s,fun2))
        return(c(f_grad))
      }
      q <- try(optim(c(factor_scores), x=new.factor_coefs_0,b=new.factor_coefs_j,s=new.sigma,g=new.gamma, method = "BFGS", fn = f_f_ini, gr = f_grad_f_ini, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_scores <- factor_scores;
      }else{
        if(iter > 1 && f.cur.logfunc > q$value){if(trace)
          cat("Optimization of m did not improve on iteration step ",iter,"\n");
          new.factor_scores <- factor_scores;
        }else{
          if(trace) cat("Variational parameters m updated","\n")
          new.factor_scores <- matrix(q$par,n.s,n.factors);
          if(q$convergence != 0) { if(trace) cat("Optimization of m did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update sigma
      beta2 <- new.factor_coefs_j^2
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
      sum <- M*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0)))
      for(i in 1:n.s){
        if(n.factors==1){
          new.sigma[i] <- 1/(1+(apply((sum)[i,]*beta2,2,sum)))
        }else{
          new.sigma[i,] <- 1/(1+(diag(diag(apply((sum)[i,]*beta2,2,sum)))))
        }}

      ###update beta
      beta_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lf*I(cvsample==0)
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        y <- sum(y1)+sum(y2)

        return(y)
      }
      beta_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        b3 <- NULL
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- M*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0)))
        for(p in 1:n.factors){
          b3 <- c(b3,(sweep((X*I(cvsample==0)),1,new.factor_scores[,p],"*") -sweep((new.sigma[,p]%*%t(new.factor_coefs_j[,p])),1,new.factor_scores[,p],"+") * sum))
        }

        b3 <- matrix(b3,n.s,n.f*n.factors)
        return(c(colSums(b3)))
      }

      q <- try(optim(c(factor_coefs_j), x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,g=new.gamma, method = "BFGS", fn = beta_f_ini, gr = beta_grad_f_ini, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_j <- factor_coefs_j;
      }else{
        if(iter > 1 && b.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_j <- factor_coefs_j;
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_j <- matrix(q$par,n.f,n.factors);
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update beta0
      b0_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lab*I(cvsample==0)
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        y <- sum(y1)+sum(y2)
        return(y)
      }
      b0_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X*I(cvsample==0)-M*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0)))
        return(c(colSums(grad)))
      }

      q <- try(optim(c(factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,g=new.gamma,method = "BFGS", fn = b0_f_ini, gr = b0_grad_f_ini, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_0 <- factor_coefs_0
      }else{
        if(iter > 1 && b0.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta0 did not improve on iteration step ",iter,"\n");
          new.factor_coefs_0 <- factor_coefs_0
        }else{
          if(trace) cat("Model parameters beta0 updated","\n")
          new.factor_coefs_0 <- q$par;
          if(q$convergence != 0) { if(trace) cat("Optimization of beta0 did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update alpha,beta0
      gam_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lab*I(cvsample==0)
        y2 <- -X*log(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        y <- sum(y1)+sum(y2)

        return(y)
      }
      gam_grad_f_ini <- function(x,b=NULL,f=NULL,s=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X*I(cvsample==0)*Y-M*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0)))*Y

        return(c(colSums(grad)))
      }

      q <- try(optim(c(gamma),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,x=new.factor_coefs_0,method = "BFGS", fn = gam_f_ini, gr = gam_grad_f_ini, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.gamma <- gamma
      }else{
        if(iter > 1 && ga.cur.logfunc > q$value){if(trace)
          cat("Optimization of gamma did not improve on iteration step ",iter,"\n");
          new.gamma <- gamma
        }else{
          if(trace) cat("Model parameters gamma updated","\n")
          new.gamma <- q$par;
          if(q$convergence != 0) { if(trace) cat("Optimization of gamma did not converge on iteration step ", iter,"\n") }
        }
      }

      q1 <- list(value = beta_f_ini(c(new.factor_coefs_j),x=new.factor_coefs_0,f=new.factor_scores,g=new.gamma,s=new.sigma))
      b.new.logfunc <- q1$value
      b.cur.logfunc <- b.new.logfunc

      q2 <- list(value = f_f_ini(c(new.factor_scores),x=new.factor_coefs_0, b=new.factor_coefs_j,g=new.gamma,s=new.sigma))
      new.f.cur.logfunc <- q2$value
      f.cur.logfunc <- new.f.cur.logfunc

      q3 <- list(value = b0_f_ini(c(new.factor_coefs_0), f=new.factor_scores,b=new.factor_coefs_j,g=new.gamma,s=new.sigma))
      new.b0.cur.logfunc <- q3$value
      b0.cur.logfunc <- new.b0.cur.logfunc

      q4 <- list(value = gam_f_ini(c(new.gamma), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,x=new.factor_coefs_0))
      new.ga.cur.logfunc <- q4$value
      ga.cur.logfunc <- new.ga.cur.logfunc

      ## Take values of VLB to define stopping rule
      q <- list(value = VLB_ini(c(new.factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,g=new.gamma,s=new.sigma))
      new.VLB <- q$value
      diff=abs(new.VLB-cur.VLB)
      ratio <- abs(new.VLB/cur.VLB);
      if(trace) cat("New VLB:", new.VLB,"cur VLB:", cur.VLB, "Ratio of VLB", ratio, ". Difference in VLB:",diff,"\n")
      cur.VLB <- new.VLB

      factor_coefs_0 <-new.factor_coefs_0
      factor_coefs_j <- new.factor_coefs_j
      factor_scores <- new.factor_scores
      sigma <- new.sigma
      gamma <- new.gamma
      iter <- iter + 1
    }

    ###VA iteration
    cur.VLB <- -1e6; iter <- 1; ratio <- 10; diff=1e5;eps = 1e-4;max.iter = 100;
    b.cur.logfunc <- -1e6;b0.cur.logfunc <- -1e6;f.cur.logfunc <- -1e6;ga.cur.logfunc <- -1e6

    while((diff> eps*(abs(cur.VLB)+eps)) && iter <= max.iter) {
      if(trace) cat("Iteration:", iter, "\n")
      ## VLB
      VLB <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,e=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.eta <- e
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*ll*I(cvsample==0)
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))}
        fun3 <- function(i) {
          new.eta <- new.eta+1e-8
          new.pi <- new.pi+1e-8
          pi.mat <-  (1-(new.pi*I(cvsample==0))[i,])*log((1-new.eta)/(1-(new.pi*I(cvsample==0))[i,]))+ (new.pi*I(cvsample==0))[i,]*log(new.eta/(new.pi*I(cvsample==0))[i,])
          sum(na.omit(pi.mat))
        }
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))+sum(sapply(1:n.s,fun3))
        return(y)
      }

      ## VLB_rept
      VLB_rept <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,e=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.eta <- e
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*ll*I(cvsample==1)
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==1)
        fun2 <- function(i) { 0.5 * (sum(log(t(new.sigma)[,i])) - sum(t(new.sigma)[,i]) - sum(new.factor_scores[i,]^2))}
        fun3 <- function(i) {
          new.eta <- new.eta+1e-8
          new.pi <- new.pi+1e-8
          pi.mat <-  (1-(new.pi*I(cvsample==1))[i,])*log((1-new.eta)/(1-(new.pi*I(cvsample==1))[i,]))+ (new.pi*I(cvsample==1))[i,]*log(new.eta/(new.pi*I(cvsample==1))[i,])
          sum(na.omit(pi.mat))
        }
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))+sum(sapply(1:n.s,fun3))
        return(y)
      }


      ###optim pi
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
      alp <- log(M/(rowSums((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))))))
      e.mat <- ll + 0.5 * (new.sigma) %*% t(new.factor_coefs_j^2)+matrix(alp,n.s,n.f)
      # e.mat <- ll
      sum1 <- exp( - exp(e.mat))*I(cvsample==0)
      for(i in 1:n.s){
        new.pi[i,] <- new.eta/(new.eta+(1-new.eta)*sum1[i,]+1e-8)
      }
      new.pi[X!=0]=0

      ###optim f
      f_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lf*I(cvsample==0)
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        fun2 <- function(i) { 0.5 * (- sum(new.factor_scores[i,]^2))}
        y <- sum(y1)+sum(y2)+sum(sapply(1:n.s,fun2))

        return(y)
      }
      f_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0)))
        fun2 <- function(i) { -new.factor_scores[i,]+((X*I(cvsample==0))[i,])%*%new.factor_coefs_j-(M[i]*sum[i,])%*%new.factor_coefs_j}
        f_grad <- t(sapply(1:n.s,fun2))
        return(c(f_grad))
      }

      q <- try(optim(c(factor_scores), x=new.factor_coefs_0,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,g=new.gamma, method = "BFGS", fn = f_f, gr = f_grad_f, control = list(trace = 0,  fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_scores <- factor_scores;
      }else{
        if(iter > 1 && f.cur.logfunc > q$value){if(trace)
          cat("Optimization of m did not improve on iteration step ",iter,"\n");
          new.factor_scores <- factor_scores;
        }else{
          if(trace) cat("Variational parameters m updated","\n")
          new.factor_scores <- matrix(q$par,n.s,n.factors);
          if(q$convergence != 0) { if(trace) cat("Optimization of m did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update sigma
      beta2 <- new.factor_coefs_j^2
      ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
      sum <- M*I((1-new.pi)>0.5)*((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0))))
      for(i in 1:n.s){
        if(n.factors==1){
          new.sigma[i] <- 1/(1+(apply((sum)[i,]*beta2,2,sum)))
        }else{
          new.sigma[i,] <- 1/(1+(diag(diag(apply((sum)[i,]*beta2,2,sum)))))
        }}

      ###update beta
      beta_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lf <- new.factor_scores %*% t(new.factor_coefs_j)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lf*I(cvsample==0)
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        y <- sum(y1)+sum(y2)
        return(y)
      }
      beta_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        b3 <- NULL
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        sum <- M*I((1-new.pi)>0.5)*((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0))))
        for(p in 1:n.factors){
          b3 <- c(b3,(sweep((X*I(cvsample==0)),1,new.factor_scores[,p],"*") -sweep((new.sigma[,p]%*%t(new.factor_coefs_j[,p])),1,new.factor_scores[,p],"+") * sum))
        }
        b3 <- matrix(b3,n.s,n.f*n.factors)
        return(c(colSums(b3)))
      }

      q <- try(optim(c(factor_coefs_j), x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,pi=new.pi,  g=new.gamma,method = "BFGS", fn = beta_f, gr = beta_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_j <- factor_coefs_j;
      }else{
        if(iter > 1 && b.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta did not improve on iteration step ",iter,"\n");
          new.factor_coefs_j <- factor_coefs_j;
        }else{
          if(trace) cat("Model parameters beta updated","\n")
          new.factor_coefs_j <- matrix(q$par,n.f,n.factors);
          if(q$convergence != 0) { if(trace) cat("Optimization of beta did not converge on iteration step ", iter,"\n") }
        }
      }

      ###update alpha,beta0
      b0_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        y1 <- X*lab*I(cvsample==0)
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)
        y <- sum(y1)+sum(y2)
        return(y)
      }
      b0_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X*I(cvsample==0)-M*I((1-new.pi)>0.5)*((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*((exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))*I(cvsample==0))))
        return(c(colSums(grad)))
      }

      q <- try(optim(c(factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma, pi=new.pi,g=new.gamma,method = "BFGS", fn = b0_f, gr = b0_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.factor_coefs_0 <- factor_coefs_0
      }else{
        if(iter > 1 && b0.cur.logfunc > q$value){if(trace)
          cat("Optimization of beta0 did not improve on iteration step ",iter,"\n");
          new.factor_coefs_0 <- factor_coefs_0
        }else{
          if(trace) cat("Model parameters beta0 updated","\n")
          new.factor_coefs_0 <- q$par
          if(q$convergence != 0) { if(trace) cat("Optimization of beta0 did not converge on iteration step ", iter,"\n") }
        }
      }
      ###update gamma
      ga_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        lab <- matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y
        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)

        y1 <- X*lab*I(cvsample==0)
        y2 <- -X*log(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2)))))*I(cvsample==0)

        y <- sum(y1)+sum(y2)

        return(y)
      }
      ga_grad_f <- function(x,b=NULL,f=NULL,s=NULL,pi=NULL,g=NULL) {
        new.factor_coefs_0 <- x[1:n.f];
        new.factor_coefs_j <- matrix(b,n.f,n.factors)
        new.factor_scores <- matrix(f,n.s,n.factors)
        new.sigma <- matrix(s,n.s,n.factors)
        new.pi <- matrix(pi,n.s,n.f)
        new.gamma <- g[1:n.f];

        ll <- matrix(new.factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(new.gamma,n.s,n.f,byrow=TRUE)*Y+new.factor_scores %*% t(new.factor_coefs_j)
        grad <- X*I(cvsample==0)*Y-M*I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))/(rowSums(I((1-new.pi)>0.5)*(exp(ll+0.5*(new.sigma) %*% t(new.factor_coefs_j^2))*I(cvsample==0))))*Y

        return(c(colSums(grad)))
      }


      q <- try(optim(c(gamma),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma, pi=new.pi,x=new.factor_coefs_0, method = "BFGS", fn = ga_f, gr = ga_grad_f, control = list(trace = 0, fnscale = -1, maxit = maxit)), silent = TRUE)
      if("try-error" %in% class(q)){ new.gamma <- gamma
      }else{
        if(iter > 1 && ga.cur.logfunc > q$value){if(trace)
          cat("Optimization of gamma did not improve on iteration step ",iter,"\n");
          new.gamma <- gamma
        }else{
          if(trace) cat("Model parameters gamma updated","\n")
          new.gamma <- q$par
          if(q$convergence != 0) { if(trace) cat("Optimization of gamma did not converge on iteration step ", iter,"\n") }
        }
      }

      new.eta <-  apply(I(new.pi>0.5)*new.pi*I(cvsample==0), 2, function(x) {sum(x)/n.s})

      q1 <- list(value = beta_f(c(new.factor_coefs_j),x=new.factor_coefs_0,f=new.factor_scores,s=new.sigma,pi=new.pi,g=new.gamma))
      b.new.logfunc <- q1$value
      b.cur.logfunc <- b.new.logfunc

      q2 <- list(value = f_f(c(new.factor_scores),x=new.factor_coefs_0, b=new.factor_coefs_j,s=new.sigma,pi=new.pi,g=new.gamma))
      new.f.cur.logfunc <- q2$value
      f.cur.logfunc <- new.f.cur.logfunc

      q3 <- list(value = b0_f(c(new.factor_coefs_0), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,g=new.gamma))
      new.b0.cur.logfunc <- q3$value
      b0.cur.logfunc <- new.b0.cur.logfunc

      q4 <- list(value = ga_f(c(new.gamma), f=new.factor_scores,b=new.factor_coefs_j,s=new.sigma,pi=new.pi,x=new.factor_coefs_0))
      new.ga.cur.logfunc <- q4$value
      ga.cur.logfunc <- new.ga.cur.logfunc

      ## Take values of VLB to define stopping rule
      q <- list(value = VLB(c(new.factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,pi=new.pi,e=new.eta,g=new.gamma))
      new.VLB <- q$value
      diff=abs(new.VLB-cur.VLB)
      ratio <- abs(new.VLB/cur.VLB);
      if(trace) cat("New VLB:", new.VLB,"cur VLB:", cur.VLB, "Ratio of VLB", ratio, ". Difference in VLB:",diff,"\n")
      cur.VLB <- new.VLB

      if(!is.null(cvsample)){
        q <- list(value = VLB_rept(c(new.factor_coefs_0),b=new.factor_coefs_j,f=new.factor_scores,s=new.sigma,pi=new.pi,e=new.eta,g=new.gamma))
        cur.VLB_rept <- q$value
      }
      factor_coefs_0 <-new.factor_coefs_0
      factor_coefs_j <- new.factor_coefs_j
      pi <- new.pi
      eta <- new.eta
      factor_scores <- new.factor_scores
      sigma <- new.sigma
      gamma <- new.gamma
      iter <- iter + 1
    }


    #ll <- matrix(factor_coefs_0,n.s,n.f,byrow=TRUE) +factor_scores %*% t(factor_coefs_j)
    ll <- matrix(factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(gamma,n.s,n.f,byrow=TRUE)*Y+factor_scores %*% t(factor_coefs_j)
    exp.mat <- exp(ll+0.5*(sigma) %*% t(factor_coefs_j^2))*I(cvsample==0)
    sum <- exp.mat/(rowSums(exp.mat))
    sum2 <- exp(ll)*I(cvsample==0)/rowSums(exp(ll)*I(cvsample==0))

    factor_scores <- scale(new.factor_scores)
    ll <- matrix(factor_coefs_0,n.s,n.f,byrow=TRUE)+matrix(gamma,n.s,n.f,byrow=TRUE)*Y+factor_scores %*% t(factor_coefs_j)
    sum3 <- exp(ll)*I(cvsample==0)/rowSums(exp(ll)*I(cvsample==0))


    if(iter > 99){
      print("ZILNMVA Not converging!")
    }

    ## print the output
    out.list$VLB <- cur.VLB
    if(!is.null(cvsample)){out.list$VLB_rept <- cur.VLB_rept}
    out.list$iter=iter-1
    out.list$lvs$pi <- pi
    out.list$lvs$factor_scores <- factor_scores
    out.list$lvs$sigma <- sigma
    out.list$params$eta <- eta
    out.list$params$gamma <- gamma
    out.list$params$factor_coefs_j <- factor_coefs_j
    out.list$params$factor_coefs_0 <- factor_coefs_0
    out.list$Q <- sum
    out.list$Q2 <- sum2
    out.list$Q3 <- sum3

    return(out.list)
  }

  if(d_choice==FALSE){
    re <- ZILNMVA(X,V)
  }else{
    if (parallel){
      #cl <- parallel::makeCluster(detectCores(logical = FALSE))
      cl <- parallel::makeCluster(getOption("cl.cores", 4))
      doParallel::registerDoParallel(cl)
    }
    out.list <- list()
    p <- ncol(X)
    n <- nrow(X)
    rank <- 5
    fold <- 5
    if(d_choice=="BIC"){

      beta <- list()
      beta0 <- list()
      f <- list()
      Q <- list()
      Q2 <- list()
      Q3 <- list()
      pi <- list()
      eta <- list()
      sigma <- list()
      gamma <- list()
      L <- rep(0,rank)
      G_w <- rep(0,rank);bic <- rep(0,rank);

      if (parallel){
        Mres <- foreach(w=1:rank) %dopar% {
          re <- ZILNMVA(X,V,n.factors=w)
          re
        }
        for(w in 1:rank){
          L[w] <- Mres[[w]]$VLB
          beta[[w]] <- Mres[[w]]$params$factor_coefs_j
          beta0[[w]] <- Mres[[w]]$params$factor_coefs_0
          eta[[w]] <- Mres[[w]]$params$eta
          gamma[[w]] <- Mres[[w]]$params$gamma
          f[[w]] <- Mres[[w]]$lvs$factor_scores
          sigma[[w]] <- Mres[[w]]$lvs$sigma
          Q[[w]] <- Mres[[w]]$Q
          Q2[[w]] <- Mres[[w]]$Q2
          Q3[[w]] <- Mres[[w]]$Q3
          pi[[w]] <- Mres[[w]]$lvs$pi
          G_w[w] <- w*p-w^2+2*w*n
          bic[w] <- -2*L[w]+(log(n)+log(p))*G_w[w]
        }
      }else{
        for(w in 1:rank){
          re <- ZILNMVA(X,V,n.factors=w)
          L[w] <- re$VLB
          beta[[w]] <- re$params$factor_coefs_j
          beta0[[w]] <- re$params$factor_coefs_0
          eta[[w]] <- re$params$eta
          gamma[[w]] <- re$params$gamma
          sigma[[w]] <- re$lvs$sigma
          f[[w]] <- re$lvs$factor_scores
          Q[[w]] <- re$Q
          Q2[[w]] <- re$Q2
          Q3[[w]] <- re$Q3
          pi[[w]] <- re$lvs$pi
          G_w[w] <- w*p-w^2+2*w*n
          bic[w] <- -2*L[w]+(log(n)+log(p))*G_w[w]
        }
      }
      out.list$VLB <- L[[(which.min(bic))]]
      out.list$lvs$pi <- pi[[(which.min(bic))]]
      out.list$lvs$factor_scores <- f[[(which.min(bic))]]
      out.list$lvs$sigma <- sigma[[(which.min(bic))]]
      out.list$params$factor_coefs_j <- beta[[(which.min(bic))]]
      out.list$params$factor_coefs_0 <- beta0[[(which.min(bic))]]
      out.list$params$eta <- eta[[(which.min(bic))]]
      out.list$params$gamma <- gamma[[(which.min(bic))]]
      out.list$Q <- Q[[(which.min(bic))]]
      out.list$Q2 <- Q2[[(which.min(bic))]]
      out.list$Q3 <- Q3[[(which.min(bic))]]
      out.list$bic <- which.min(bic)

    }
    if(d_choice=="CV"){

      cvs<- NULL
      for (i in 1:fold){
        cvs <- c(cvs, i * rep(1, floor(n*p/fold)))
      }
      cvs <- sample(cvs, length(cvs), replace = FALSE)
      cvs <- matrix(cvs, nrow = n, ncol = p)

      All_rept <- NULL
      for (rept in 1:fold){
        cvsample <- cvs
        cvsample[cvsample==rept] <- 1
        cvsample[cvsample!=1] <- 0
        # X_cv <- X[cvsample==0]
        # X_rept <- X[cvsample==1]
        L_rept <- matrix(0, nrow = 1, ncol = rank)
        if (parallel){
          Mres <- foreach(w=1:rank)%dopar% {
            re <- ZILNMVA(X,V, n.factors=w, cv_group=cvsample)
            re$VLB_rept
          }
          L_rept[1,(1:rank)] <- Mres
        }else{
          for(w in 1:rank){
            re <- ZILNMVA(X,V, n.factors=w, cv_group=cvsample)
            L_rept[1,w] <- re$VLB_rept
          }
        }
        All_rept <- rbind(All_rept, L_rept)
      }
      cv <- which.max(colSums(matrix(unlist(All_rept),fold,rank)))

      re <- tryCatch({ZILNMVA(X,V,n.factors=cv)},error=function(e){NaN})

      out.list$cv <- cv
      out.list$VLB <- re$VLB[[cv]]
      out.list$params$factor_coefs_j <- re$params$factor_coefs_j[[cv]]
      out.list$params$factor_coefs_0 <-  re$params$factor_coefs_0[[cv]]
      out.list$params$eta <- re$params$eta[[cv]]
      out.list$params$gamma <- re$params$gamma[[cv]]
      out.list$lvs$sigma <- re$lvs$sigma[[cv]]
      out.list$lvs$factor_scores <- re$lvs$factor_scores[[cv]]
      out.list$lvs$pi <-  re$lvs$pi[[cv]]
      out.list$Q <-re$Q[[cv]]
      out.list$Q2 <-re$Q2[[cv]]
      out.list$Q3 <-re$Q3[[cv]]

    }
    if (parallel){
      parallel::stopCluster(cl = cl)
    }
    return(out.list)

  }
}
