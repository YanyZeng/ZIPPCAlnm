# ZIPPCAlnm
Zero-Inflated Probabilistic PCA framework with logistical normal multinomial distribution (ZIPPCA-LNM)

# Installation
```r
install.packages("devtools")  
devtools::install_github("YanyZeng/ZIPPCAlnm")  
library(ZIPPCAlnm) 
```

# Description
Zero-Inflated Probabilistic PCA framework with logistical normal multinomial distribution (ZIPPCA-LNM), that extends probabilistic PCA from the Gaussian setting to multivariate abundance data, and an empirical Bayes approach was proposed for inferring microbial compositions. An efficient VA algorithm, classification VA, has been developed for fitting this model.

# Usage
```r
ZIPPCAlnm(X, V = NULL, d_choice = FALSE, parallel = TRUE)
```
* X: count matrix of observations.
* V: vector of sample covariate.
* d_choice: FALSE, "BIC" or "CV". Indicating whether the rank or number of factors, is chosen from 1 to 5. Options are "BIC" (Bayesian information criterion), and "CV" (Cross-validation). BIC is recommended. Defaults to FALSE.
* parallel: logical, if TRUE, use parallel toolbox to accelerate.


# Examples
```r
n.n = 50
n.w = 100
n.factors = 2
set.seed(1)
si <- diag(n.factors)
me <- c(0,0)
f <- matrix(0,nrow = n.n, ncol = n.factors)
for(i in 1:n.n){
 f[i,] <- mvrnorm(1,me,si)
}
betaj <- matrix(0,nrow = n.w, ncol = n.factors)
for(j in 1:n.w){
  betaj[j,] <- runif(n.factors,-3.5,3.5)
}
beta0 <- rep(2,n.w)
l <- matrix(beta0,n.n,n.w,byrow=TRUE)+f \%*\%  t(betaj)
eta_j <- rep(0.25,n.w)
z <- matrix(0,n.n,n.w,byrow = TRUE)
for(i in 1:n.n){
  z[i,] <- rbinom(n.w, size=1, prob=eta_j)
}
sum <- rowSums(exp(l))
Qn <- exp(l)/sum
sum <- matrix(rowSums((1-z)*exp(l)),n.n,n.w)
Qn_z <- matrix(I(z==0)*(exp(l)/sum)+I(z==1)*0,n.n,n.w)
X <- matrix(0,n.n,n.w,byrow = TRUE)
for(i in 1:n.n){
  X[i,] <- rmultinom(1, size = runif(1,800,1000), prob = Qn_z[i,])
}
zerorow <- which(rowSums(X)==0)
if(length(zerorow) >0 ){
  X <- X[-zerorow,];Qn <- Qn[-zerorow,];Qn_z <- Qn_z[-zerorow,];
}
zerocol <- which(colSums(X)==0)
if(length(zerocol) >0 ){
  X <- X[,-zerocol];Qn <- Qn[,-zerocol];Qn_z <- Qn_z[,-zerocol];
}
result <- ZIPPCAlnm::ZIPPCAlnm(X,V=NULL,d_choice=FALSE,parallel=TRUE)
 ```
