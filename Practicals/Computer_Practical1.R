#Computer Practical 1
#Introduction to high dimensional statistics
################################################################################

#Divergence of the pairwise Euclidean distances in high dimensions

set.seed(1)

n <- 100
p <- 2
X <- matrix(0,n,p)
for(i in 1:n)
 {
 for(j in 1:p)
 {
  X[i,j] <- runif(1,0,1)  
 }
}
pairwise_dist <- dist(X)
hist(pairwise_dist, xlim=c(0,max(pairwise_dist)))

modified_pairwise_dist <- dist(X)/sqrt(p)
hist(modified_pairwise_dist, xlim=c(0,max(modified_pairwise_dist)))

################################################################################

#Increasing mean squared error of the OLS method as dimension p gets larger

set.seed(1)
n <- 100
p <- 2
X <- matrix(0,n,p)
Y <- numeric(n)
beta <- numeric(p)
for(i in 1:n)
{
 for(j in 1:p)
 {
   X[i,j] <- 3*sin(-pi*j*i/n) 
 }
}

for(j in 1:p)
{
  beta[j] <- rnorm(1,0,j^(-4)) 
}

for(i in 1:n)
{
  Y[i] <- t(X[i,])%*%beta+rnorm(1,0,1)
}

#an alternative way of getting simulated responses Y in matrix form
Y <- X%*%beta+rnorm(n,0,1)

LS <- lm(Y~-1+X)
summary(LS)
beta_LS <- LS$coefficients

Yhat <- numeric(n)
Ytrue <- numeric(n)
for(i in 1:n)
{
  Yhat[i] <- t(X[i,])%*%beta_LS
  Ytrue[i] <- t(X[i,])%*%beta
}

#a simpler way of calculating the fitted and true response values
Yhat <- X%*%beta_LS
Ytrue <- X%*%beta

t <- (1:n)/n
plot(x=t,y=Y, ylim=c(-4,5))
points (x=t, y=Ytrue, type="l", lty=2, col="blue")
points (x=t, y=Yhat, type="l", lty=1, col="red")
legend(x="bottomleft",          # Position
       legend=c("OLS","signal"),  # Legend texts
       lty=c(1, 2),           # Line types
       col=c(2, 4),           # Line colors
       lwd=2)                 # Line width
  
################################################################################  



