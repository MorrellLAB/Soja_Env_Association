#   R script to create a folded SFS for all markers on the SoySNP50K for
#   all genotyped soja and per cluster
#   Written by Thomas JY Kono

#   A function to calculate the minor allele frequency
#   The files are already in 0/1 format, so this should be easy!
MAF <- function(x) {
    #   First, we strip the missing data
    no_missing <- x[!is.na(x)]
    #   Then count up how many of each allele there is
    counts <- table(no_missing)
    #   Which is the minimum?
    minor_allele <- min(counts)
    #   And then the frequency
    maf <- minor_allele / length(no_missing)
    #   Return it
    #   There are some that are monomorphic
    #   Return NA in this case
    if(maf == 1) {
        return(NA)
    }
    else {
        return(maf)
    }
}

#   Read in all the separate matrices
#   Markers are columns, individuals are rows
#   The first column has an integer that describes population of origin
soja_all <- read.table("32k_ATCG.txt", header=T, row.names=2)
soja_pop1 <- soja_all[soja_all$Pop == 1,]
soja_pop2 <- soja_all[soja_all$Pop == 2,]
soja_pop3 <- soja_all[soja_all$Pop == 3,]


#   Get the minor allele frequencies for all markers in all partitions
all_mafs <- apply(soja_all, 2, MAF)
pop1_mafs <- apply(soja_pop1, 2, MAF)
pop2_mafs <- apply(soja_pop2, 2, MAF)
pop3_mafs <- apply(soja_pop3, 2, MAF)

#   Then, we bin them up
bins <- seq(0, 0.5, by=0.05)
all_mafs_bins <- cut(all_mafs, breaks=bins, include.lowest=TRUE)
pop1_mafs_bins <- cut(pop1_mafs, breaks=bins, include.lowest=TRUE)
pop2_mafs_bins <- cut(pop2_mafs, breaks=bins, include.lowest=TRUE)
pop3_mafs_bins <- cut(pop3_mafs, breaks=bins, include.lowest=TRUE)

#   Count up how many are in each bin
all_mafs_SFS <- table(all_mafs_bins) / length(all_mafs)
pop1_mafs_SFS <- table(pop1_mafs_bins) / length(pop1_mafs)
pop2_mafs_SFS <- table(pop2_mafs_bins) / length(pop2_mafs)
pop3_mafs_SFS <- table(pop3_mafs_bins) / length(pop3_mafs)

#   Stick them together into a data frame so we can plot side-by-side
SFSData <- as.data.frame(
    cbind(
        all_mafs_SFS,
        pop1_mafs_SFS,
        pop2_mafs_SFS,
        pop3_mafs_SFS)
    )

#   Open a handle to a PDF file for the plot
pdf(
    file="Soja_FoldedSFS_ByCluster.pdf",
    width=8,
    height=6,
    family="Helvetica",
    pointsize=16)
#   Make the barplot, we have to do a separate x-axis later
barplot(t(SFSData),
    ylim=c(0,0.20),
    beside=TRUE,
    axisnames=F,
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    col=c("black", "green2", "red", "blue"))
labels <- c(
    "0-0.05",
    "0.05-0.10",
    "0.10-0.15",
    "0.15-0.20",
    "0.20-0.25",
    "0.25-0.30",
    "0.30-0.35",
    "0.35-0.40",
    "0.40-0.45",
    "0.45-0.50")

#   You'll have to play with these numbers to get the labels to line up with
#   the tick marks. As far as I know, there isn't a good way to figure this out
#   other than trial and error
#at <- c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48)
#at <- c(3, 8, 13, 18, 23)
#   Create an axis
#axis(side=1,
#    at=at,
#    las=1,
#    labels=labels,
#    font=1,
#    padj=1,
#    cex.axis=0.8)

#   And a legend
legend(inset=0,
    cex=1.1,
    "topright",
    c("All", "Mainland South", "Island", "Mainland North"),
    fill=c("black", "green2", "red", "blue"))
dev.off ()
