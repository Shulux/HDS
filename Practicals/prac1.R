#Q1

# X1, ..., X_n
# X_i = (X_i(1), ... X_i(p))
# i = 1, ..., n

# pairwise distances with L2 norm

n <- 100
p <- 1000
X <- matrix(0, n, p)
for (i in 1:n){
  for (j in 1:p){
    X[i, j] <- runif(1,0,1)
  }
}

pairwise_distances <- dist(X, method = "euclidean")

hist(pairwise_distances, xlim=c(0,max(pairwise_distances)))

modified_pairwise_dist <- dist(X)/sqrt(p)
hist(modified_pairwise_dist, xlim=c(0,1))

# we can see with the modified distance, the pairwise distances converge rather than diverge.


#Q2
set.seed(1)
n <- 100
p <- 2

X <- matrix(0, n, p)
for (i in 1:n){
  for (j in 1:p){
    X[i, j] <- 3*sin(-pi*j*i/n)
  }
}

Y <- rep(0, n)
beta <- c(rep(0, p))
epsilon <- rep(0, n)
for (i in 1:n){
  epsilon[i] <- rnorm(1, 0, 1)
}
for (j in 1:p){
  beta[j] <- rnorm(1, 0, j^(-4))
}


for (i in 1:n){
  sum_val <- 0
  for (j in 1:p){
    sum_val <- sum_val + 3 * beta[j] * sin(-pi*j*i/n) + epsilon[i]
  }
  Y[i] <- sum_val
}

LS <- lm(Y~-1+X)
summary(LS)
beta_LS <- LS$coefficients
