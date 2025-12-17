##### HDS Mini Project Report
##### author: qvns53 - Nic Salmon
##### github: https://github.com/Shulux/HDS/tree/main


library(sparsepca)
library(fossil)

# save all plots to a pdf file automatically
pdf(file = "coursework/hds_plots.pdf", height=7, width=8)

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

rownames(Tumor) <- make.names(row_names, unique=TRUE)

###### Data Exploration
######
dim(Tumor) # n=64, p=6830
###### since the data is very high dimensional,
###### i am only looking at 10 randomly sampled dimensions to reduce output bloat

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




#K-means clustering
#km.out <- kmeans(x=Tumor, centers=3, iter.max=100, nstart=1)
#km.out
#print(km.out$cluster)
#summary(km.out)

#plot(Tumor, col=km.out$cluster, main="K-means with 3 clusters", xlab="", ylab="")


# end saving plots to pdf
dev.off()
