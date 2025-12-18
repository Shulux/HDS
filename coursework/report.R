##### HDS Mini Project Report
##### author: qvns53 - Nic Salmon
##### github: https://github.com/Shulux/HDS/tree/main


library(sparsepca)
library(fossil)

# save all plots to a pdf file automatically
#pdf(file = "coursework/figures/hds_plots.pdf", height=7, width=8)
png(file="coursework/figures/hds_plots_%d.png",type="cairo",width=1800,height=1600,res=200)

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


n <- nrow(X) #64
p <- ncol(X) #6830


# including all the genes in clustering may dilute the distances as the dimension is increased.
# as we know as p -> infinity, Euclidean distances diverge, so we want to lower p.

dists <- dist(X)
summary(dists)
hist(dists, breaks=50, main="Pairwise distances between samples", xlab="Euclidean distance")

### K-Means

k_max <- 15
wss <- numeric(k_max)
ari_vals <- numeric(k_max)

for (k in 1:k_max) {
  km <- kmeans(X, centers = k, iter.max=100, nstart = 20)
  wss[k] <- km$tot.withinss
  ari_vals[k] <- adj.rand.index(km$cluster, labels)
}

plot(1:k_max, wss, type = "b",
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow plot for k-means clustering",
     xaxt = 'n')
axis(1, at=c(1:k_max))

#plot(2:k_max, ari_vals[2:k_max], type = "b",
#     xlab = "Number of clusters (k)",
#     ylab = "Adjusted Rand Index",
#     main = "ARI vs k clusters")
#axis(1, at=c(2:k_max))

# chosen to be k=6, but honestly the elbow/scree plot does not show a clear value to choose for k.
# went further by using the adjust rand index values plot, this helped the decision for optimal k.
k_chosen <- 6

set.seed(1)
km_raw <- kmeans(X, centers = k_chosen, nstart = 50)

clusters_km <- km_raw$cluster

# evaluating k means clustering with adjusted rand index
ari_km_raw <- adj.rand.index(clusters_km, labels)
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

ari_hcA_raw <- adj.rand.index(clusters_hcA, labels)
ari_hcA_raw

ari_hcS_raw <- adj.rand.index(clusters_hcS, labels)
ari_hcS_raw

ari_hcC_raw <- adj.rand.index(clusters_hcC, labels)
ari_hcC_raw


### K-Medians


###### Basic clustering evaluation


###### Basic dimensionality reduction


###### Basic dimensionality reduction evaluation


###### Clustering after dim reduction


###### Clustering after dim reduction evaluation



# end saving plots to pdf
dev.off()
