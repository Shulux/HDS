##### HDS Mini Project Report
##### author: qvns53 - Nic Salmon
##### github:


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
