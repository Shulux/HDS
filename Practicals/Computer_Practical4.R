#Computer Practical 4
#High dimensional cluster analysis
#####################################################################
#High dimensional clustering for Human tumor microarray data

library(sparsepca)
library(fossil)
#library(foreign)
#library(pdfCluster)

set.seed(1)

#Load human tumor microarray data
data <- read.csv("Practicals/Human_Tumor_Microarray.csv",header=TRUE)
Tumor <- data[,-1]
Tumor <- t(Tumor)#to make it as an nxp data matrix with p variables in columns

Tumor <- as.data.frame(Tumor)

#from the file of labels available on Ultra
row_names <- c("CNS", "CNS", "CNS", "RENAL", "BREAST", "CNS", "CNS", "BREAST", "NSCLC",
                      "NSCLC", "RENAL", "RENAL", "RENAL", "RENAL", "RENAL", "RENAL", "RENAL",
                      "BREAST", "NSCLC", "RENAL", "UNKNOWN", "OVARIAN", "MELANOMA",
                      "PROSTATE", "OVARIAN", "OVARIAN", "OVARIAN", "OVARIAN",
                      "OVARIAN", "PROSTATE", "NSCLC", "NSCLC", "NSCLC", "LEUKEMIA",
                      "K562A-repro", "K562A-repro", "LEUKEMIA", "LEUKEMIA",
                      "LEUKEMIA", "LEUKEMIA", "LEUKEMIA", "COLON", "COLON", "COLON",
                      "COLON", "COLON", "COLON", "COLON", "K562A-repro", "BREAST",
                      "K562A-repro", "BREAST", "NSCLC", "NSCLC", "NSCLC", "MELANOMA",
                      "BREAST", "BREAST", "MELANOMA", "MELANOMA", "MELANOMA",
                      "MELANOMA", "MELANOMA", "MELANOMA")

rownames(Tumor) <- make.names(row_names, unique=TRUE)

Tumor <- scale(Tumor, center=TRUE, scale=FALSE)

#K-means clustering
km.out <- kmeans(x=Tumor, centers=3, iter.max=100, nstart=1)
km.out
print(km.out$cluster)
summary(km.out)

plot(Tumor, col=km.out$cluster, main="K-means with 3 clusters", xlab="", ylab="")

#Handling random algorithms
# Set up a 2x3 plotting grid
par(mfrow = c(2, 3))

for(i in 1:6) {
  # Run kmeans() on x with three clusters and one start
  km.out <- kmeans(x=Tumor, centers=3, iter.max=100, nstart=1)
  # Plot clusters
  plot(x=Tumor, col=km.out$cluster, main=km.out$tot.withinss)
}

#Selecting an optimal number of clusters
# Initialize total within sum of squares error (wss)
wss <- 0

# For 1 to 15 cluster centers
for (i in 1:10) {
  km.out <- kmeans(x=Tumor, centers=i, iter.max=100, nstart=20)
  # Save total within sum of squares to wss variable
  wss[i] <- km.out$tot.withinss

}

# Plot total within sum of squares vs. number of clusters
plot(1:10, wss, type = "b",
     xlab = "Number of clusters",
     ylab = "Sum of squares")

k <- 3
plot(Tumor[, c(5,24)],
     col = km.out$cluster,
     main = paste("k-means clustering of Tumor with", k, "clusters"),
     xlab = "Defense", ylab = "Speed")

#Hierarchical clustering
# Clustering using complete linkage
hc.out1 <- hclust(dist(Tumor), method="complete")

# Clustering using average linkage
hc.out2 <- hclust(dist(Tumor), method="average")

# Clustering using single linkage
hc.out3 <- hclust(dist(Tumor), method="single")

# Plot dendrogram of Clustering using average linkage
plot(hc.out2, main="Average")

#hierarchical clustering using dimensionality reduction via PCA
pr.out <- prcomp(x=Tumor, center=FALSE, scale=FALSE)
summary(pr.out)
proj.data <- pr.out$x[,1:42] #The first 42 leading PCs
proj.data <- Tumor%*%pr.out$rotation[,1:42]
hc.out4 <- hclust(dist(proj.data), method="average")
plot(hc.out4, main="Average")

#K-means clustering using dimensionality reduction via PCA
pr.out <- prcomp(x=Tumor, center=FALSE, scale=FALSE)
summary(pr.out)
proj.data <- Tumor%*%pr.out$rotation[,1:42]
km.out.pca <- kmeans(x=proj.data, centers=3, iter.max=100, nstart=20)

#Sparse PCA
spca.out <- spca(Tumor,k=20,alpha=1e-3,beta=1e-3,center=TRUE,scale=FALSE,verbose=0)
print(spca.out$loadings)
View(spca.out$loadings)
summary(spca.out)

#library(elasticnet)
spca.out2 <- spca(x=Tumor, K=10, para=rep(0.5,10), type=c("predictor","Gram"),
                  sparse=c("penalty","varnum"), use.corr=FALSE, lambda=1e-6,
                  max.iter=200, trace=FALSE, eps.conv=1e-3)
print(spca.out2$loadings)

###################################################################

#High dimensional clustering for S&P500 stock data

set.seed(1)

load(file="SP500data.rda")
dimnames(SP500data) <- list(rownames(SP500data, do.NULL=FALSE, prefix="T"), colnames(SP500data, do.NULL=FALSE, prefix="S"))
SP500data <- as.data.frame(SP500data)
SP500data <- scale(SP500data,center=TRUE,scale=TRUE)

#K-means clustering
km.out <- kmeans(x=SP500data, centers=3, iter.max=100, nstart=20)
km.out
print(km.out$cluster)
summary(km.out)

plot(SP500data, col=km.out$cluster, main="K-means with 3 clusters", xlab="", ylab="")

#Handling random algorithms
# Set up a 2x3 plotting grid
par(mfrow = c(2, 3))

for(i in 1:6) {
  # Run kmeans() on x with three clusters and one start
  km.out <- kmeans(x=SP500data, centers=3, iter.max=100, nstart=1)
  # Plot clusters
  plot(x=SP500data, col=km.out$cluster, main=km.out$tot.withinss)
}

#Selecting an optimal number of clusters
# Initialize total within sum of squares error (wss)
wss <- 0

# For 1 to 15 cluster centers
for (i in 1:10) {
  km.out <- kmeans(x=SP500data, centers=i, iter.max=100, nstart=20)
  # Save total within sum of squares to wss variable
  wss[i] <- km.out$tot.withinss

}

# Plot total within sum of squares vs. number of clusters
plot(1:10, wss, type = "b",
     xlab = "Number of clusters",
     ylab = "Sum of squares")

k <- 3
plot(SP500data[, c(5,24)],
     col = km.out$cluster,
     main = paste("k-means clustering of SP500data with", k, "clusters"),
     xlab = "Defense", ylab = "Speed")

#hierarchical clustering
# Clustering using complete linkage
hc.out1 <- hclust(dist(SP500data), method="complete")

# Clustering using average linkage
hc.out2 <- hclust(dist(SP500data), method="average")

# Clustering using single linkage
hc.out3 <- hclust(dist(SP500data), method="single")

# Plot dendrogram of Clustering using average linkage
plot(hc.out2, main="Average")

#clustering using the dimensionality reduction by PCA
pr.out <- prcomp(x=SP500data, center=FALSE, scale=FALSE)
summary(pr.out)
proj.data <- SP500data%*%pr.out$rotation[,1:7]#leading PCs
proj.data <-pr.out$x[,1:7]#leading PCs, the same as above line
#K-means
km.out.pca <- kmeans(x=proj.data, centers=3, iter.max=100, nstart=20)
plot(proj.data, col=km.out.pca$cluster, main="K-means with 3 clusters", xlab="", ylab="")
#hierarchical
hc.out.pca <- hclust(dist(proj.data), method="average")
plot(hc.out.pca, main="Average")

####################################################################
#Dimensionality reduction with PCA
# Perform scaled PCA: pr.out
pr.out <- prcomp(x=SP500data, center=FALSE, scale=FALSE)
summary(pr.out)
biplot(pr.out)
biplot(pr.out,choices=1:2,expand=20,xlim=c(-1.5,1.0), ylim=c(-1.0,1.5))

#pca results with a smaller dimension
pr.out2 <- prcomp(x=SP500data[,sample(1:496,10)], center=FALSE, scale=FALSE)
summary(pr.out2)
biplot(pr.out2)
#Note that R function "princomp" does not work for high dimensional data, to see this try princomp(x=SP500data)

#PCA with dimensionality reduction, here only keeping first 5 components
pr.out2 <- prcomp(x=SP500data, center=FALSE, scale=FALSE, rank.=5)
summary(pr.out2)
biplot(pr.out2)

# Variability of each principal component: pr.var
pr.var <- (pr.out$sdev)^2

# Variance explained by each principal component: pve
pve <- pr.var / sum(pr.var)

# Plot variance explained for each principal component
plot(pve, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     ylim = c(0, 1), type = "b")

# Plot cumulative proportion of variance explained
plot(cumsum(pve), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     ylim = c(0, 1), type = "b")


#Rand index
rand.index(as.vector(km.out$cluster), as.vector(km.out.pca$cluster))
adj.rand.index(as.vector(km.out$cluster), as.vector(km.out.pca$cluster))

