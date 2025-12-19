##### HDS Mini Project Report
##### author: qvns53 - Nic Salmon
##### github: https://github.com/Shulux/HDS/tree/main


library(sparsepca)
library(fossil)

# save all plots to a pdf file automatically
pdf(file = "coursework/figures/hds_plots.pdf", height=7, width=8)
#png(file="coursework/figures/hds_plots_%d.png",type="cairo",width=1800,height=1600,res=200)

set.seed(1)

#Load human tumor microarray data
data <- read.csv("Practicals/Human_Tumor_Microarray.csv",header=TRUE)
Tumor <- data[,-1]
Tumor <- t(Tumor) #to make it as an nxp data matrix with p variables in columns

Tumor <- as.data.frame(Tumor)

for (i in (1:dim(Tumor)[2])) {colnames(Tumor)[i] <- paste0("g", i)}

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

labels <- factor(row_names, levels = names(sort(table(row_names), TRUE)))
str(labels)
table(labels)
# the classes are not balanced, eg rarer tumor type example is "PROSTATE", a common type is "RENAL".

library(ggplot2)

#barplot(prop.table(table(labels)), main="Tumor Label Distribution", xlab="Tumor Label", ylab="Proportion", las=2)

labels_df <- data.frame(labels)

ggplot(labels_df, aes(x = labels)) +
  geom_bar(fill = "steelblue") + # Adds bars, colors them
  scale_y_continuous(labels = function(x) paste0(round(x/64*100,2), "%"), breaks = seq(0, 9, 1)) +
  labs(title = "Tumor Label Distribution",
       x = "Tumor Label", y = "Proportion")


rownames(Tumor) <- make.names(row_names, unique=TRUE)

###### Data Exploration
######
dim(Tumor) # n=64, p=6830
###### since the data is very high dimensional,
###### i am only looking at 10 randomly sampled dimensions to reduce output bloat

sum(is.na(Tumor)) # no missing values


sampleRows <- sample(1:dim(Tumor)[1], size=5, replace=FALSE)
sampleCols <- sample(1:dim(Tumor)[2], size=10, replace=FALSE)
sampleColsBig <- sample(1:dim(Tumor)[2], size=20, replace=FALSE)

sampledTumor <- Tumor[,sampleCols]

head(sampledTumor[sampleRows,])  # show only the sampled 5 rows, sampled 10 columns
# values seem to be small, continous, positive and negative

str(sampledTumor) # check structure (sampled 10 columns)
# all numeric variables g1->g6380

summary(sampledTumor) # summary statistics of sampled 10 columns
# values are small, seem to be centered roughly around zero, and small in magnitude


### always center the data by subtracting the mean

boxplot(Tumor[, sampleColsBig], main="Subset Gene Expression marginalised over patients", xlab="Gene", ylab="Gene Expression Value", las=2)



Tumor <- scale(Tumor, center=TRUE, scale=FALSE)
copyTumor <- as.data.frame(Tumor)


centered_sampledTumor <- copyTumor[,sampleCols]

head(centered_sampledTumor[sampleRows,])
str(centered_sampledTumor)
summary(centered_sampledTumor)
# centered correctly



## the variable columns g1->g6380 represent the expression values for the genes (for the corresponding sample)
## which is an experimental observation done in a lab setting to quantify how active the gene is
## this reflects the amount of mRNA produced, which indicates functionality.
## positive values indicated more expression and negative indicated less expression,
## which is still true after the centering.
##
## the observation/sample rows represent a human patient with a tumor
## each sample has a corresponding label which identifies the type of tumor, eg. "MELANOMA".


###### Basic Data Analysis

# variable intervals and summary statistics of the random subset
summary(centered_sampledTumor)
sum_stats <- as.data.frame(apply(centered_sampledTumor, 2, summary))

# visualisation of summary statistics through boxplots
boxplot(Tumor[, sampleColsBig], main="Subset Gene Expression marginalised over patients (centered)", xlab="Gene", ylab="Gene Expression Value", las=2)

boxplot(scale(Tumor, center=TRUE, scale=TRUE)[, sampleColsBig], main="Subset Gene Expression marginalised over patients (centered and scaled)", xlab="Gene", ylab="Gene Expression Value", las=2)



# distribution plots (histograms) of a few of the randomly selected variables (genes)
par(mfrow = c(2, 3))
for(i in 1:6) {
  randCol <- centered_sampledTumor[, i]
  randColName <- colnames(centered_sampledTumor)[i]
  hist(randCol, xlab= "", ylab="Density", main=paste(randColName , "Distribution"), freq=FALSE)
  dens <- density(randCol)
  lines(dens)
}
par(mfrow = c(1, 1))

# we can see that this subset of genes, seem to have unimodal distributions,
# which are centered on zero as expected.


# pairs plot of the random subset (reduced to 6 columns for clarity)
pairs_Tumor <- centered_sampledTumor[, 1:6]
pairs(pairs_Tumor)
# does not seem to be any clear pattern or distribution within the subset data
# although


# using the library pheatmap, we can see that for some specific patient samples,
# eg for BREAST.2, we see a very high gene expression for g597, potentially indicating
# a correlation between gene 597 expression and breast tumors. or g5307 with OVARIAN.5
# similarly, some genes seem to have (almost) invariate expressions across all patient samples,
# eg. g2347 has only a few samples with a moderate gene expression, which may suggest
# that it is an 'inactive' variable for tumor presence.
library(pheatmap)
tumor_scaled <- as.matrix(scale(centered_sampledTumor))
pheatmap::pheatmap(tumor_scaled, treeheight_row = 0, treeheight_col = 0)



gene_var <- apply(Tumor, 2, var)
summary(gene_var)
hist(gene_var, breaks=50, main="Gene-wise variance distribution")
abline(v = quantile(gene_var, 0.9), col="red", lwd=2)


# many genes have very little variance and so are likely uninformative for clustering
# only a few of the genes have higher and more substantial variances and so are potentially the most informative genes

# this gives evidence to the idea that the signal is concentrated
# in a low-dimensional subspace, with the genes that have high variance/information
# most genes contribute mostly noise
# a few genes contribute signal
# supports the use of PCA and sparse PCA


###### Basic clustering

X <- as.matrix(Tumor)          # centered, unscaled, full data
labels <- factor(row_names) # previously defined, but now without sorting
K_true <- length(unique(labels)) # we have 11 unique label classes (clusters)

# fossil-compatible label vector (integer cluster IDs)
labels_id <- as.integer(labels)


n <- nrow(X) #64
p <- ncol(X) #6830


# including all the genes in clustering may dilute the distances as the dimension is increased.
# as we know as p -> infinity, Euclidean distances diverge, so we want to lower p.

dists <- dist(X)
summary(dists)
hist(dists, breaks=50, main="Pairwise distances between samples", xlab="Euclidean distance")

### K-Means

set.seed(1)
k_max <- 15
wss <- numeric(k_max)
ari_vals <- numeric(k_max)

for (k in 1:k_max) {
  km <- kmeans(X, centers = k, iter.max=100, nstart = 20)
  wss[k] <- km$tot.withinss
  ari_vals[k] <- adj.rand.index(km$cluster, labels_id)
}

plot(1:k_max, wss, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow plot for k-means clustering",
     xaxt = 'n')
axis(1, at=c(1:k_max))

plot(1:k_max, ari_vals[1:k_max], type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Adjusted Rand Index",
     main = "ARI vs k clusters")
axis(1, at=c(1:k_max))

# the elbow/scree plot does not show a clear value to choose for k.
# roughly between 4-7
# went further by using the adjust rand index values plot,
# the peak is at k=5, (potential other peaks at k=3, k=6)
# this helped the decision for optimal k.
# chosen to be k=5 from looking at both plots.
k_chosen <- 5

set.seed(1)
km_raw <- kmeans(X, centers = k_chosen, nstart = 50)

clusters_km <- km_raw$cluster

# evaluating k means clustering with adjusted rand index
ari_km_raw <- adj.rand.index(clusters_km, labels_id)
ari_km_raw



# 2D g1533 g4567 plot of k means clustering
plot(X[, sampleCols[1:2]],
     col = clusters_km,
     pch = 19,
     main = paste("k-means clusters with", k_chosen, "clusters (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)

#using PCA (first 2 PCs) to get a visualisation, not for dim reduction purposes yet.
plot(prcomp(X)$x[,1:2],
     col = clusters_km,
     pch = 19,
     main = "K-means clusters (first 2 PCs, raw clustering)")
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)


### Hierarchical clustering
set.seed(1)
## Average linkage
hc.outA <- hclust(dist(X), method="average")
## simple linkage
hc.outS <- hclust(dist(X), method="single")
## complete linkage
hc.outC <- hclust(dist(X), method="complete")

# Plot dendrograms of Clustering
plot(hc.outA, main="Average")
plot(hc.outS, main="Single")
plot(hc.outC, main="Complete")

clusters_hcA <- cutree(hc.outA, k = k_chosen)
clusters_hcS <- cutree(hc.outS, k = k_chosen)
clusters_hcC <- cutree(hc.outC, k = k_chosen)

ari_hcA_raw <- adj.rand.index(clusters_hcA, labels_id)
ari_hcA_raw

ari_hcS_raw <- adj.rand.index(clusters_hcS, labels_id)
ari_hcS_raw

ari_hcC_raw <- adj.rand.index(clusters_hcC, labels_id)
ari_hcC_raw

round(c(ari_hcA_raw, ari_hcS_raw, ari_hcC_raw),5)


plot(X[, sampleCols[1:2]],
     col = clusters_hcA,
     pch = 19,
     main = paste("Average linkage clusters with", k_chosen, "cluster cutoff (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hcS,
     pch = 19,
     main = paste("Single linkage clusters with", k_chosen, "cluster cutoff (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hcC,
     pch = 19,
     main = paste("Complete linkage clusters with", k_chosen, "cluster cutoff (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)


plot(prcomp(X)$x[,1:2],
     col = clusters_hcA,
     pch = 19,
     main = "Average linkage clusters (first 2 PCs, 5 clusters cut-off)")
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)


plot(prcomp(X)$x[,1:2],
     col = clusters_hcS,
     pch = 19,
     main = "Single linkage clusters (first 2 PCs, 5 clusters cut-off)")
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)

plot(prcomp(X)$x[,1:2],
     col = clusters_hcC,
     pch = 19,
     main = "Complete linkage clusters (first 2 PCs, 5 clusters cut-off)")
legend("topright",
       legend = paste("Cluster", 1:k_chosen),
       col = 1:k_chosen,
       pch = 19,
       cex = 0.8)


###### Dimensionality reduction

# PCA on centered data
pca_raw <- prcomp(X, center = FALSE, scale. = FALSE)

# basic outputs
pca_raw$sdev        # singular values
pca_raw$rotation    # loadings
pca_raw$x           # scores


# variance explained
var_explained <- pca_raw$sdev^2
pve <- var_explained / sum(var_explained)
cum_pve <- cumsum(pve)

# scree plot to choose the number of PCs to use in further analysis
plot(cum_pve, type = "b",
     xlab = "Number of principal components",
     ylab = "Cumulative proportion of variance explained",
     main = "Cumulative variance explained by PCA",
     yaxt = 'n',
     xaxt = 'n')
axis(1, at=seq(0, 60, 5))
axis(1, at=c(25, 42, 50, 64))
axis(2, at=seq(0, 1, 0.1))
abline(h = c(0.75,0.9, 0.95), col = c("red", "blue", "green"), lty = 2)
abline(v = c(25, 42, 50), col = c("red", "blue", "green"), lty = 2)

summary(pca_raw)
print(cum_pve)

# choose ks such that 90/95% of the variance is explained
k75 <- which(cum_pve >= 0.75)[1]
k80 <- which(cum_pve >= 0.80)[1]
k85 <- which(cum_pve >= 0.85)[1]
k90 <- which(cum_pve >= 0.90)[1]
k95 <- which(cum_pve >= 0.95)[1]
k975 <- which(cum_pve >= 0.975)[1]
k75; k80; k85; k90; k95; k975

# picking k90, which is k=42 clusters, as we want to reduce from available 64 PCs

# PCA-reduced data to use for clustering
X_pca_75 <- pca_raw$x[, 1:k75]
X_pca_90 <- pca_raw$x[, 1:k90]
X_pca_95 <- pca_raw$x[, 1:k95]



#Sparse PCA
set.seed(1)
K_sp <- 20
spca.out <- spca(X, k=K_sp, alpha=1e-3, beta=1e-3, center=TRUE, scale=FALSE, verbose=0)

summary(spca.out)

Z_sp <- X %*% spca.out$loadings


#################################### Clustering after applying dim reduction on reduced data Z

# picked proportion of variance
k75 # first 25 PCs hold 75% of the total variance, which is a lot lower than even 42
k90 # first 42 PCs hold 90% of the total variance
Z <- X_pca_90
Z_75 <- X_pca_75

nrow(Z) #64
ncol(Z) #42 PCs heavily reduced dimensionality from 6830 >> 42 > 25
ncol(Z_75)


head(Z[sampleRows,1:5]) # each patient observation is represented on the new basis (PCs)
head(Z_75[sampleRows,1:5])


dists_pca_sp <- dist(Z_sp)
summary(dists_pca_sp)
hist(dists_pca_sp, breaks=50, main="Pairwise distances between samples on reduced data Z_sp", xlab="Euclidean distance")

dists_pca75 <- dist(X_pca_75)
summary(dists_pca75)
hist(dists_pca75, breaks=50, main="Pairwise distances between samples on reduced data with Z_75", xlab="Euclidean distance")


dists_pca90 <- dist(X_pca_90)
summary(dists_pca90)
hist(dists_pca90, breaks=50, main="Pairwise distances between samples on reduced data with Z_90", xlab="Euclidean distance")

dists_pca95 <- dist(X_pca_95)
summary(dists_pca95)
hist(dists_pca95, breaks=50, main="Pairwise distances between samples on reduced data with Z_95", xlab="Euclidean distance")



####### kmeans after PCA

set.seed(1)
wss_pca <- numeric(k_max)
ari_vals_pca <- numeric(k_max)

for (k in 1:k_max) {
  km <- kmeans(Z, centers = k, iter.max=100, nstart = 20)
  wss_pca[k] <- km$tot.withinss
  ari_vals_pca[k] <- adj.rand.index(km$cluster, labels_id)
}

plot(1:k_max, wss_pca, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow plot for k-means clustering on reduced data Z",
     xaxt = 'n')
axis(1, at=c(1:k_max))

plot(1:k_max, ari_vals_pca[1:k_max], type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Adjusted Rand Index",
     main = "ARI vs k clusters on reduced data Z")
axis(1, at=c(1:k_max))

which(wss != wss_pca)
# always a lower WSS_PCA than WSS
which(ari_vals != ari_vals_pca) # the same clustering results for k=1,2,3,4,5,6,7,10
# not the same for k=8,9,11,12,14,15


# repeated but for Z_75 instead

set.seed(1)
wss_pca_75 <- numeric(k_max)
ari_vals_pca_75 <- numeric(k_max)

for (k in 1:k_max) {
  km <- kmeans(Z_75, centers = k, iter.max=100, nstart = 20)
  wss_pca_75[k] <- km$tot.withinss
  ari_vals_pca_75[k] <- adj.rand.index(km$cluster, labels_id)
}

plot(1:k_max, wss_pca_75, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow plot for k-means clustering on reduced data Z_75",
     xaxt = 'n')
axis(1, at=c(1:k_max))

plot(1:k_max, ari_vals_pca_75[1:k_max], type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Adjusted Rand Index",
     main = "ARI vs k clusters on reduced data Z_75")
axis(1, at=c(1:k_max))

which(wss != wss_pca_75)
which(wss_pca != wss_pca_75)
# always a lower WSS_PCA than WSS
which(ari_vals != ari_vals_pca_75) # the same clustering results for some k
which(ari_vals_pca != ari_vals_pca_75)
# not the same for k=5  6  7  8  9 11 12 13 14 15



# repeated but for Z_75 instead

set.seed(1)
wss_pca_sp <- numeric(k_max)
ari_vals_pca_sp <- numeric(k_max)

for (k in 1:k_max) {
  km <- kmeans(Z_sp, centers = k, iter.max=100, nstart = 20)
  wss_pca_sp[k] <- km$tot.withinss
  ari_vals_pca_sp[k] <- adj.rand.index(km$cluster, labels_id)
}

plot(1:k_max, wss_pca_sp, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow plot for k-means clustering on reduced data Z_sp",
     xaxt = 'n')
axis(1, at=c(1:k_max))

plot(1:k_max, ari_vals_pca_sp[1:k_max], type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Adjusted Rand Index",
     main = "ARI vs k clusters on reduced data Z_sp")
axis(1, at=c(1:k_max))

which(wss != wss_pca_sp)
which(wss_pca != wss_pca_sp)
# always a lower WSS_PCA than WSS
which(ari_vals != ari_vals_pca_sp) # the same clustering results for k <= 6
which(ari_vals_pca != ari_vals_pca_sp)

# we are getting the almost the same scree and ARI plots on the reduced data with k90,
# so PCA made only a small change here.

k_chosen_PCA <- 5
k_chosen_PCA_75 <- 5
k_chosen_PCA_sp <- 5

set.seed(1)
km_pca_90 <- kmeans(Z, centers = k_chosen_PCA, nstart = 50)
km_pca_75 <- kmeans(Z_75, centers = k_chosen_PCA_75, nstart = 50)
km_pca_sp <- kmeans(Z_sp, centers = k_chosen_PCA_sp, nstart = 50)
km_pca_95 <- kmeans(X_pca_95, centers = k_chosen_PCA, nstart = 50)

clusters_km_90 <- km_pca_90$cluster
clusters_km_75 <- km_pca_75$cluster
clusters_km_sp <- km_pca_sp$cluster
clusters_km_95 <- km_pca_95$cluster

# checking if there is a difference in actual clustering after PCA
sum(clusters_km_90 == clusters_km)  ##  == 64, which means no change for k=5 clusters
sum(clusters_km_75 == clusters_km)  ##  == 0, which means all of them have changed
sum(clusters_km_sp == clusters_km)  ##  == 8, which means only a few have not changed
# there is change for chosen k values though

km_pca_75$size
km_pca_90$size
km_pca_sp$size
km_raw$size

ari_km_pca_90 <- adj.rand.index(
  as.integer(clusters_km_90),
  labels_id
)

ari_km_pca_75 <- adj.rand.index(
  as.integer(clusters_km_75),
  labels_id
)

ari_km_pca_sp <- adj.rand.index(
  as.integer(clusters_km_sp),
  labels_id
)

ari_km_pca_95 <- adj.rand.index(
  as.integer(clusters_km_95),
  labels_id
)

ari_km_raw
ari_km_pca_90
ari_km_pca_95
ari_km_pca_75
ari_km_pca_sp

# all the same ARI as the clustering results for k=5 is the same before and after PCA

# hierarchical after PCA
set.seed(1)
## Average linkage
hc.pcaA <- hclust(dist(Z), method="average")
## simple linkage
hc.pcaS <- hclust(dist(Z), method="single")
## complete linkage
hc.pcaC <- hclust(dist(Z), method="complete")

# Plot dendrograms of Clustering
plot(hc.pcaA, main="Average linkage on reduced data Z")
plot(hc.pcaS, main="Single linkage on reduced data Z")
plot(hc.pcaC, main="Complete linkage on reduced data Z")

clusters_hc_pca_A <- cutree(hc.pcaA, k = k_chosen_PCA)
clusters_hc_pca_S <- cutree(hc.pcaS, k = k_chosen_PCA)
clusters_hc_pca_C <- cutree(hc.pcaC, k = k_chosen_PCA)


ari_hcA_pca <- adj.rand.index(clusters_hc_pca_A, labels_id)
ari_hcA_pca

ari_hcS_pca <- adj.rand.index(clusters_hc_pca_S, labels_id)
ari_hcS_pca

ari_hcC_pca <- adj.rand.index(clusters_hc_pca_C, labels_id)
ari_hcC_pca

round(c(ari_hcA_pca, ari_hcS_pca, ari_hcC_pca),5) # complete is again the best performing ARI

## REPEAT FOR Z_75

set.seed(1)
## Average linkage
hc.pcaA_75 <- hclust(dist(Z_75), method="average")
## simple linkage
hc.pcaS_75 <- hclust(dist(Z_75), method="single")
## complete linkage
hc.pcaC_75 <- hclust(dist(Z_75), method="complete")

# Plot dendrograms of Clustering
plot(hc.pcaA_75, main="Average linkage on reduced data Z_75")
plot(hc.pcaS_75, main="Single linkage on reduced data Z_75")
plot(hc.pcaC_75, main="Complete linkage on reduced data Z_75")

clusters_hc_pca_A_75 <- cutree(hc.pcaA_75, k = k_chosen_PCA_75)
clusters_hc_pca_S_75 <- cutree(hc.pcaS_75, k = k_chosen_PCA_75)
clusters_hc_pca_C_75 <- cutree(hc.pcaC_75, k = k_chosen_PCA_75)


ari_hcA_pca_75 <- adj.rand.index(clusters_hc_pca_A_75, labels_id)
ari_hcA_pca_75

ari_hcS_pca_75 <- adj.rand.index(clusters_hc_pca_S_75, labels_id)
ari_hcS_pca_75

ari_hcC_pca_75 <- adj.rand.index(clusters_hc_pca_C_75, labels_id)
ari_hcC_pca_75

round(c(ari_hcA_pca_75, ari_hcS_pca_75, ari_hcC_pca_75),5) # complete is again the best performing ARI

## REPEAT FOR Z_sp

set.seed(1)
## Average linkage
hc.pcaA_sp <- hclust(dist(Z_sp), method="average")
## simple linkage
hc.pcaS_sp <- hclust(dist(Z_75), method="single")
## complete linkage
hc.pcaC_sp <- hclust(dist(Z_sp), method="complete")

# Plot dendrograms of Clustering
plot(hc.pcaA_sp, main="Average linkage on reduced data Z_sp")
plot(hc.pcaS_sp, main="Single linkage on reduced data Z_sp")
plot(hc.pcaC_sp, main="Complete linkage on reduced data Z_sp")

clusters_hc_pca_A_sp <- cutree(hc.pcaA_sp, k = k_chosen_PCA_sp)
clusters_hc_pca_S_sp <- cutree(hc.pcaS_sp, k = k_chosen_PCA_sp)
clusters_hc_pca_C_sp <- cutree(hc.pcaC_sp, k = k_chosen_PCA_sp)


ari_hcA_pca_sp <- adj.rand.index(clusters_hc_pca_A_sp, labels_id)
ari_hcA_pca_sp

ari_hcS_pca_sp <- adj.rand.index(clusters_hc_pca_S_sp, labels_id)
ari_hcS_pca_sp

ari_hcC_pca_sp <- adj.rand.index(clusters_hc_pca_C_sp, labels_id)
ari_hcC_pca_sp

round(c(ari_hcA_pca_sp, ari_hcS_pca_sp, ari_hcC_pca_sp),5) # average is the best performing ARI for Z_sp




##################### visualisations of reduced data Z clustering algorithm results

# the randomly sampled genes plot for the reduced kmeans clustering using Z
plot(X[, sampleCols[1:2]],
     col = clusters_km_90,
     pch = 19,
     main = paste("Reduced data Z k-means clusters with", k_chosen_PCA, "clusters (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)

# the randomly sampled genes plot for the reduced kmeans clustering using Z_75
plot(X[, sampleCols[1:2]],
     col = clusters_km_75,
     pch = 19,
     main = paste("Reduced data Z_75 k-means clusters with", k_chosen_PCA_75, "clusters (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)

# the randomly sampled genes plot for the reduced kmeans clustering using Z_sp
plot(X[, sampleCols[1:2]],
     col = clusters_km_sp,
     pch = 19,
     main = paste("Reduced data Z_sp k-means clusters with", k_chosen_PCA_sp, "clusters (using 2 randomly sampled genes)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)

# along the first 2 PCs again for Z
plot(Z[,1:2],
     col = clusters_km_90,
     pch = 19,
     main = "K-means clustering (k=5) on reduced data Z",
     xlab = "PC1", ylab = "PC2")

legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA, pch = 19, cex = 0.8)

# along the first 2 PCs again for Z_75
plot(Z[,1:2],
     col = clusters_km_75,
     pch = 19,
     main = "K-means clustering (k=5) on reduced data Z_75",
     xlab = "PC1", ylab = "PC2")

legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA, pch = 19, cex = 0.8)


# along the first 2 PCs again for Z_sp
plot(Z[,1:2],
     col = clusters_km_sp,
     pch = 19,
     main = "K-means clustering (k=5) on reduced data Z_sp",
     xlab = "PC1", ylab = "PC2")

legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp, pch = 19, cex = 0.8)
#================================================================================ hierarchical visualisations


plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_A,
     pch = 19,
     main = paste("Average linkage with", k_chosen_PCA, "cluster cutoff on Z (sample)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_S,
     pch = 19,
     main = paste("Single linkage with", k_chosen_PCA, "cluster cutoff on Z (sample)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_C,
     pch = 19,
     main = paste("Complete linkage with", k_chosen_PCA, "cluster cutoff on Z (sample)"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)


plot(Z[,1:2],
     col = clusters_hc_pca_A,
     pch = 19,
     main = "Average linkage clusters (first 2 PCs, 5 cluster cut-off) on Z")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)


plot(Z[,1:2],
     col = clusters_hc_pca_S,
     pch = 19,
     main = "Single linkage clusters (first 2 PCs, 5 cluster cut-off) on Z")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)

plot(Z[,1:2],
     col = clusters_hc_pca_C,
     pch = 19,
     main = "Complete linkage clusters (first 2 PCs, 5 cluster cut-off) on Z")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA),
       col = 1:k_chosen_PCA,
       pch = 19,
       cex = 0.8)



############ REPEATED FOR Z_75

plot(X[, sampleCols[1:2]],
       col = clusters_hc_pca_A_75,
       pch = 19,
       main = paste("Average linkage with", k_chosen_PCA_75, "cluster cutoff on Z_75"),
       xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_S_75,
     pch = 19,
     main = paste("Single linkage with", k_chosen_PCA_75, "cluster cutoff on Z_75"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_C_75,
     pch = 19,
     main = paste("Complete linkage with", k_chosen_PCA_75, "cluster cutoff on Z_75"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)


plot(Z[,1:2],
     col = clusters_hc_pca_A_75,
     pch = 19,
     main = "Average linkage clusters (first 2 PCs, 5 cluster cut-off) on Z_75")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)


plot(Z[,1:2],
     col = clusters_hc_pca_S_75,
     pch = 19,
     main = "Single linkage clusters (first 2 PCs, 5 cluster cut-off) on Z_75")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)

plot(Z[,1:2],
     col = clusters_hc_pca_C_75,
     pch = 19,
     main = "Complete linkage clusters (first 2 PCs, 5 cluster cut-off) on Z_75")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_75),
       col = 1:k_chosen_PCA_75,
       pch = 19,
       cex = 0.8)



############ REPEATED FOR Z_sp

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_A_sp,
     pch = 19,
     main = paste("Average linkage with", k_chosen_PCA_sp, "cluster cutoff on Z_sp"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_S_sp,
     pch = 19,
     main = paste("Single linkage with", k_chosen_PCA_sp, "cluster cutoff on Z_sp"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)

plot(X[, sampleCols[1:2]],
     col = clusters_hc_pca_C_sp,
     pch = 19,
     main = paste("Complete linkage with", k_chosen_PCA_sp, "cluster cutoff on Z_sp"),
     xlab = paste0("g", sampleCols[1]), ylab = paste0("g", sampleCols[2]))
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)


plot(Z[,1:2],
     col = clusters_hc_pca_A_sp,
     pch = 19,
     main = "Average linkage clusters (first 2 PCs, 5 cluster cut-off) on Z_sp")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)


plot(Z[,1:2],
     col = clusters_hc_pca_S_sp,
     pch = 19,
     main = "Single linkage clusters (first 2 PCs, 5 cluster cut-off) on Z_sp")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)

plot(Z[,1:2],
     col = clusters_hc_pca_C_sp,
     pch = 19,
     main = "Complete linkage clusters (first 2 PCs, 5 cluster cut-off) on Z_sp")
legend("topright",
       legend = paste("Cluster", 1:k_chosen_PCA_sp),
       col = 1:k_chosen_PCA_sp,
       pch = 19,
       cex = 0.8)


#================================================================================


## Sparse PCA Interpretation

loadings <- spca.out$loadings
rownames(loadings) <- colnames(X)

nonzero_genes <- function(loadings, j) {
  v <- loadings[, j]
  names(v)[v != 0]
}

genes_pc1 <- nonzero_genes(loadings, 1)
genes_pc2 <- nonzero_genes(loadings, 2)
genes_pc3 <- nonzero_genes(loadings, 3)

gene_list <- c(genes_pc1, genes_pc2, genes_pc3)

gene_freq <- sort(table(gene_list), decreasing = TRUE)
head(gene_freq, 20)

shared_genes <- gene_freq[gene_freq >= 2]
shared_genes

unique_pc1 <- setdiff(genes_pc1, union(genes_pc2, genes_pc3))
unique_pc2 <- setdiff(genes_pc2, union(genes_pc1, genes_pc3))
unique_pc3 <- setdiff(genes_pc3, union(genes_pc1, genes_pc2))

abs_loadings <- abs(loadings[, 1:3])
abs_sum <- rowSums(abs_loadings)

top_genes <- sort(abs_sum, decreasing = TRUE)
head(top_genes, 20)

gene_summary <- data.frame(
  gene = names(gene_freq),
  frequency = as.integer(gene_freq),
  abs_loading_sum = abs_sum[names(gene_freq)]
)

gene_summary <- gene_summary[order(-gene_summary$frequency,
                                   -gene_summary$abs_loading_sum), ]

head(gene_summary, 10)


######################## putting scree plots together
k_vals <- 1:length(wss)

wss_norm      <- wss / wss[1]
wss_pca_norm  <- wss_pca / wss_pca[1]
wss_75_norm   <- wss_pca_75 / wss_pca_75[1]
wss_sp_norm   <- wss_pca_sp / wss_pca_sp[1]

plot(k_vals, wss_norm,
     type="b", pch=19, lwd=2,
     xlab="Number of clusters (k)",
     ylab="Relative within-cluster sum of squares",
     main="Normalised elbow plots for k-means")

lines(k_vals, wss_pca_norm,  type="b", pch=17, lwd=2, lty=2)
lines(k_vals, wss_75_norm,   type="b", pch=15, lwd=2, lty=3)
lines(k_vals, wss_sp_norm,   type="b", pch=18, lwd=2, lty=4)

legend("topright",
       legend=c("Raw data",
                "PCA (90%)",
                "PCA (75%)",
                "Sparse PCA"),
       lty=c(1,2,3,4),
       pch=c(19,17,15,18),
       lwd=2)
axis(1, at=1:15)

k_vals2 <- 2:length(ari_vals)
col_raw   <- "black"
col_pca75 <- "red"


plot(k_vals2, ari_vals[k_vals2],
     type = "b",
     col = col_raw,
     lwd = 2,
     pch = 19,
     ylim = c(0, 1),
     xlab = "Number of clusters (k)",
     ylab = "Adjusted Rand Index",
     main = "ARI vs number of clusters for k-means")

lines(k_vals[k_vals2], ari_vals_pca[k_vals2],
      type = "b",
      col = col_raw,
      lwd = 1,
      lty = 2,
      pch = 19)

lines(k_vals[k_vals2], ari_vals_pca_75[k_vals2],
      type = "b",
      col = col_pca75,
      lwd = 2,
      lty = 3,
      pch = 19)

lines(k_vals[k_vals2], ari_vals_pca_sp[k_vals2],
      type = "b",
      col = col_raw,
      lwd = 2,
      lty = 4,
      pch = 19)

legend("topright",
       legend = c("Raw data",
                  "PCA (90% variance)",
                  "PCA (75% variance)",
                  "Sparse PCA"),
       col = c(col_raw, col_raw, col_pca75, col_raw),
       lty = c(1, 2, 3, 4),
       lwd = 2,
       pch = 19,
       cex = 0.9)

axis(1, at=2:15)
abline(v = 5, lty = 3, col = "grey40")
text(5, 0.95, "k = 5", pos = 4, cex = 0.8)


########### all important results collated to report

head(X[sampleRows,sampleCols])
summary(pca_raw)
summary(spca.out)
paste0("k", c(seq(75,95,5),975, 100), ": ", c(k75, k80, k85, k90, k95, k975, 64))
head(Z[sampleRows,1:5])
head(Z_75[sampleRows,1:5])
head(Z_sp[sampleRows,1:5])
wss
wss_pca
wss_pca_75
wss_pca_sp
which(wss != wss_pca)
which(wss != wss_pca_75)
which(wss != wss_pca_sp)
ari_vals
ari_vals_pca
ari_vals_pca_75
ari_vals_pca_sp
which(ari_vals != ari_vals_pca)
which(ari_vals != ari_vals_pca_75)
which(ari_vals != ari_vals_pca_sp)
km_pca_75$size
km_pca_sp$size
km_pca_90$size
km_raw$size
c(ari_km_raw, ari_km_pca_90, ari_km_pca_75, ari_km_pca_sp)
c(ari_hcA_raw, ari_hcS_raw, ari_hcC_raw)
c(ari_hcA_pca, ari_hcS_pca, ari_hcC_pca)
c(ari_hcA_pca_sp, ari_hcS_pca_sp, ari_hcC_pca_sp)
c(ari_hcA_pca_75, ari_hcS_pca_75, ari_hcC_pca_75)
head(gene_summary, 10)


top10 <- head(gene_summary, 10)

# make frequency a factor for discrete colouring
top10$frequency <- factor(top10$frequency,
                          levels = c(1, 2, 3),
                          labels = c("1 PC", "2 PCs", "3 PCs"))

ggplot(top10,
       aes(x = reorder(gene, abs_loading_sum),
           y = abs_loading_sum,
           fill = frequency)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "Gene",
    y = "Sum of absolute loadings (PC1–PC3)",
    fill = "Appears in",
    title = "Top 10 genes from sparse PCA",
    subtitle = "Ranked by total absolute loading across first three sparse PCs"
  ) +
  scale_fill_manual(
    values = c("1 PC" = "grey70",
               "2 PCs" = "steelblue",
               "3 PCs" = "darkred")
  ) +
  theme_minimal(base_size = 12)


# end saving plots to pdf
dev.off()
