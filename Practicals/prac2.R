# Prac 2

library(glmnet)

# Q2.1
# eyedata.rda
load("Practicals/eyedata.rda")

dim(eyedata)
# 120x201
# 200 genes for 120 samples plus a response

# first column Y = TRIM32 gene response variable
# find out which covariates are associated with Y, and predict Y from covariates

Y = eyedata[,1]
X = eyedata[,-1]

# scale so that they have zero mean

centered.Y <- as.vector(scale(Y, center = TRUE, scale = FALSE))
centered.X <- as.matrix(scale(X, center = TRUE, scale = FALSE))

QRdecomp <- qr(centered.X)
rankX <- QRdecomp$rank
# solve(t(centered.X)%*%centered.X)
# X is singular, so we can not calculate (t(X) X)^{-1}

summary(lm(centered.Y~-1+centered.X))
# we have 81 undefined coefficients due to singularity
# 81 because n < p by 81. n = 120, p = 201, p-n = 81

# alpha = 0 -> ridge only
# lambda = 2 defined in question

#apply ridge regression with Cross Validation for selection of lambda
M1 <- glmnet(centered.X, centered.Y, alpha=0, lambda=2)
M1$beta
CV1 <- cv.glmnet(centered.X, centered.Y, alpha=0)
lambda_ridge <- CV1$lambda.min
plot(CV1)
M1 <- glmnet(centered.X, centered.Y, alpha=0)
betahat_ridge <- coef(M1, s=lambda_ridge)
betahat_ridge_nointercept <- betahat_ridge[-1]

#apply lasso with Cross Validation for selection of lambda
CV2 <- cv.glmnet(centered.X, centered.Y, alpha=1)
lambda_lasso <- CV2$lambda.min
plot(CV2)
plot(CV2$glmnet.fit, xvar="lambda")
M2 <- glmnet(centered.X, centered.Y, alpha=1)
betahat_lasso <- coef(M2, s=lambda_lasso)
betahat_lasso_nointercept <- betahat_lasso[-1]


set.seed(2)
n <- nrow(eyedata)
train_rows <- sample(1:n, 0.7*n)
X.train <- centered.X[train_rows, ]
X.test <- centered.X[-train_rows, ]
Y.train <- centered.Y[train_rows]
Y.test <- centered.Y[-train_rows]

CV1 <- cv.glmnet(X.train, Y.train, alpha=0)
lambda_ridge_train <- CV1$lambda.min
M1 <- glmnet(X.train, Y.train, alpha=0)
Yhat_ridge <- predict(M1,X.test,s=lambda_ridge_train)
MSPE_ridge <- mean((Y.test - Yhat_ridge)^2)

CV2 <- cv.glmnet(X.train, Y.train, alpha=1)
lambda_lasso_train <- CV2$lambda.min
M2 <- glmnet(X.train, Y.train, alpha=1)
Yhat_lasso <- predict(M2,X.test,s=lambda_lasso_train)
MSPE_lasso <- mean((Y.test - Yhat_lasso)^2)


c(MSPE_ridge, MSPE_lasso)
# ridge has a lower means square prediction error so its better


# Q2.2
# riboflavin.rda
load("Practicals/riboflavin.rda")

dim(riboflavin)
# 71x2

# 71 samples with 4088 covariates
# set up so col 2 has all 4088 covariates within it.

Y <- riboflavin$y

X <- riboflavin$x

length(Y) # response 71 samples
dim(X) # design matrix X covariates 71x4088


X <- scale(X, center=TRUE, scale=FALSE)
round(colMeans(X)[1:5])

n <- nrow(X)
p <- ncol(X)

mean(Y)


# pairwise correlation between X_i, X_j

pairwise_cor <- cor(X)

pairwise_cor <- round(pairwise_cor, 4)
hist(pairwise_cor)

mean((pairwise_cor))

corrX <- pairwise_cor

library(reshape2)
corrX <- corrX[1:50,1:50]
melted_corrX <- melt(corrX)
head(melted_corrX)
library(ggplot2)
ggplot(data = melted_corrX, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 8, hjust = 0))+
  coord_fixed()



#apply ridge regression with Cross Validation for selection of lambda
library(glmnet)
CV1 <- cv.glmnet(X, Y, alpha=0)
lambda_ridge <- CV1$lambda.min
plot(CV1)
M1 <- glmnet(X, Y, alpha=0)
betahat_ridge <- coef(M1, s=lambda_ridge)

#apply lasso with Cross Validation for selection of lambda
CV2 <- cv.glmnet(X, Y, alpha=1)
lambda_lasso <- CV2$lambda.min
plot(CV2)
M2 <- glmnet(X, Y, alpha=1)
betahat_lasso <- coef(M2, s=lambda_lasso)
betahat_lasso_final <- betahat_lasso[betahat_lasso[,1]!=0,]
length(betahat_lasso_final)


# partitioning data into train, test splits.

set.seed(1)
train_rows <- sample(1:n, 0.7*n)
X.train <- X[train_rows, ]
X.test <- X[-train_rows, ]
Y.train <- Y[train_rows]
Y.test <- Y[-train_rows]

CV1 <- cv.glmnet(X.train, Y.train, alpha=0)
lambda_ridge_train <- CV1$lambda.min
M1 <- glmnet(X.train, Y.train, alpha=0)
Yhat_ridge <- predict(M1,X.test,s=lambda_ridge_train)
MSPE_ridge <- mean((Y.test - Yhat_ridge)^2)

CV2 <- cv.glmnet(X.train, Y.train, alpha=1)
lambda_lasso_train <- CV2$lambda.min
M2 <- glmnet(X.train, Y.train, alpha=1)
Yhat_lasso <- predict(M2,X.test,s=lambda_lasso_train)
MSPE_lasso <- mean((Y.test - Yhat_lasso)^2)


print(c(MSPE_ridge, MSPE_lasso))


#apply elastic net with a specific value for lambda (here lambda=0.5)
CV3 <- cv.glmnet(X, Y, alpha=0.5)
lambda_elastic <- CV3$lambda.min
plot(CV3)
M3 <- glmnet(X, Y, alpha=0.5)
betahat_elastic <- coef(M3, s=lambda_elastic)
betahat_elastic_final <- betahat_elastic[betahat_elastic[,1]!=0,]
length(betahat_elastic_final)

#apply elastic net with Cross Validation for selection of lambda
models <- list()
for (i in 0:20)
{
  name <- paste0("alpha", i/20)
  models[[name]] <- cv.glmnet(X,Y,alpha=i/20)
}

results <- data.frame()
for (i in 0:20) {
  name <- paste0("alpha", i/20)
  predicted <- predict(models[[name]],
                       s=models[[name]]$lambda.1se, newx=X.test)
  mse <- mean((Y.test - predicted)^2)
  temp <- data.frame(alpha=i/20, mse=mse, name=name)
  results <- rbind(results, temp)
}
plot(results$alpha, results$mse)

elastic_net_est <- predict(models[["alpha0.8"]], type = "coef")
elastic_net_est <- elastic_net_est[elastic_net_est[,1]!=0,]
length(elastic_net_est)

