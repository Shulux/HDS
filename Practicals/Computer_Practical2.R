#Computer Practical 2
#Regularised regression
################################################################################

#Applying the OLS, ridge and lasso regression to eyedata

set.seed(1)

load(file="eyedata.rda")
View(eyedata)
dim(eyedata)
Y <- eyedata[,1]
X <- eyedata[,-1]
X <- as.matrix(X)
X <- scale(X, center=TRUE, scale=FALSE) 
Y <- scale(Y, center=TRUE, scale=FALSE)

QRdecomp <- qr(X)
rankX <- QRdecomp$rank
solve(t(X)%*%X)

summary(lm(Y~-1+X))

#install.packages("glmnet")
library(glmnet)

#apply ridge regression with Cross Validation for selection of lambda
M1 <- glmnet(X, Y, alpha=0, lambda=2)
M1$beta
CV1 <- cv.glmnet(X, Y, alpha=0)
lambda_ridge <- CV1$lambda.min
plot(CV1)
M1 <- glmnet(X, Y, alpha=0)
betahat_ridge <- coef(M1, s=lambda_ridge)
betahat_ridge_nointercept <- betahat_ridge[-1]

#apply lasso with Cross Validation for selection of lambda
CV2 <- cv.glmnet(X, Y, alpha=1)
lambda_lasso <- CV2$lambda.min
plot(CV2)
plot(CV2$glmnet.fit, xvar="lambda")
M2 <- glmnet(X, Y, alpha=1)
betahat_lasso <- coef(M2, s=lambda_lasso)
betahat_lasso_nointercept <- betahat_lasso[-1]

set.seed(2)
n <- nrow(eyedata)
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


################################################################################

#Regularised regression for analysing the riboflavin data

set.seed(1)

load(file="riboflavin.rda")
Y <- riboflavin$y
X <- riboflavin$x
n <- nrow(X)
p <- ncol(X)

X <- scale(X,center=TRUE,scale=FALSE)

corrX <- cor(X)
corrX <- round(corrX,4)
hist(corrX)

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
  theme(axis.text.x = element_text(angle = 90, vjust = 0, 
                                   size = 8, hjust = 0))+
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

################################################################################  



