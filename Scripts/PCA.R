#   Script to perform Principal Components Analysis (PCA)
#   This data is 50K SNP data consisting of 3764 markers that have no missing
#   data, the first column is the PI name.
#   Written by Michael B. Kantar

Data <- read.table("file path to marker data", header=TRUE)
#Make sure the data is correct (names should be chromosome positions)
names(Data)

#Make sure data is correct
head(Data)
length(Data[1,])
Data[1:4,1:3]

#Remove labels from data inorder to be able to do the PCA
Data1 <- Data[, -(1:2)]

# Pricipal Components Analysis
# entering raw data and extracting PCs 
# from the correlation matrix 
fit <- prcomp(Data1, cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit, type="lines") # scree plot
fit$scores # the principal components
biplot(fit)

#Make country a factor so it can be colored in a plot
country <- as.factor(Data2[,2])

#extract first two principle components
comp <- data.frame(fit$x[,1:2])

#plot PC1 x PC2 colored by country
plot(comp, pch=19, col=country)

#Make Legend (Modified for cultivars addition)
legend(
    "bottomleft",
    title="Country",
    cex=0.75,
    pch=16,
    col=c("Black","Green","Red","Blue","Cyan","Purple"),
    legend=c("China","Japan","Cultivars","Korea","Russia", "Taiwan"),
    ncol=6)

#load data for major cluster
major_cluster <- read.table(
    "C:/Users/Michael/Desktop/Forever_Green_Grant/major_cluster.txt",
    header=T)

#Make major cluster a factor for coloring
major <- as.factor(major_cluster[,5])

#plot PC1 x PC2 colored by structure output
plot(comp, pch=19, col=country)

#Make legend
legend(
    "topright",
    title="Structure Cluster",
    cex=0.75,
    pch=16,
    col=c("Black","Green","Red"),
    legend=c("cluster 1","cluster 2","cluster 3"),
    ncol=3)

# make plot of PC1 x PC2 x PC3
comp <- data.frame(fit$x[,1:3])
plot(comp, pch=19, col=country)

# make plot of PC1 x PC2 x PC3 x PC4 x PC5
comp <- data.frame(fit3$x[,1:5])
plot(comp, pch=19, col=country)

#verify PCA with new R package for PCA
library(FactoMineR)
result <- PCA(Data1)
res.mca <- MCA(Data1)

plot.PCA(result, axes=c(1, 2), choix="ind")

res.hpc2 <- HCPC(result)

plotellipses(result)


#######################enviromnetal pca
#load bioclimatic & biphysical data
dat <- read.table("file name for data", header=TRUE)

## Pricipal Components Analysis
# entering raw data and extracting PCs 
# from the correlation matrix 

#Remove labels from data inorder to be able to do the PCA
Dat1 <- dat[, -(2:5)]

fit <- prcomp(Dat1, cor=TRUE)
summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit, type="lines") # scree plot
fit$scores # the principal components
biplot(fit)

#Make country a factor so it can be colored in a plot
country <- as.factor(dat[,3])

#extract first 2 principle components
comp <- data.frame(fit$x[,1:2])

#plot PC1 x PC2 colored by country output
plot(comp, pch=19, col=major)

library(FactoMineR)
result <- PCA(Dat1, ncp=3) # graphs generated automatically
plot(result)
get_variables <- HCPC(result)
pca_soja <- get_variables$desc.var

write.csv(
    pca_soja,
    "file path to output vaiable clusters",
    sep="\t",
    col.names=FALSE)
