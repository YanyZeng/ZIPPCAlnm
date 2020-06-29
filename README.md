# ZIPPCAlnm
Zero-Inflated Probabilistic PCA framework with logistical normal multinomial distribution (ZIPPCA-LNM)

# Installation
```r
install.packages("devtools")  
devtools::install_github("YanyZeng/ZIPPCAlnm")  
library(ZIPPCAlnm) 
```r
# Description
Zero-Inflated Probabilistic PCA framework with logistical normal multinomial distribution (ZIPPCA-LNM), that extends probabilistic PCA from the Gaussian setting to multivariate abundance data, and an empirical Bayes approach was proposed for inferring microbial compositions. An efficient VA algorithm, classification VA, has been developed for fitting this model.

# Usage
```r
ZIPPCAlnm(X, d_choice = FALSE)
```r
* X: count matrix of observations.
* d_choice: logical, if TRUE the rank or number of factors, or dimension after dimensional reduction, will be chosen from 1 to 5. Defaults to FALSE.

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
 #f[i,] <- rnorm(n.factors, mean = 0, sd = 1)
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
X <- matrix(0,n.n,n.w,byrow = TRUE)
for(i in 1:n.n){
  X[i,] <- rmultinom(1, size = runif(1,800,1000), prob = Qn[i,])
}
X[z==1] <-0
zerorow <- which(rowSums(X)==0)
if(length(zerorow) >0 ){
  X <- X[-zerorow,];Qn <- Qn[-zerorow,];
}
zerocol <- which(colSums(X)==0)
if(length(zerocol) >0 ){
  X <- X[,-zerocol];Qn <- Qn[,-zerocol];
}
result <- ZIPPCAlnm::ZIPPCAlnm(X,d_choice=FALSE)
 ```
