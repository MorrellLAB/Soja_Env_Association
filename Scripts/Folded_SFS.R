#   R script to create a folded SFS for all markers on the SoySNP50K for
#   all genotyped soja and per cluster

#   A function to calculate the minor allele frequency
#   The files are already in 0/1 format, so this should be easy!
MAF <- function(x)
{
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
    if(maf == 1)
    {
        return(NA)
    }
    else
    {
        return(maf)
    } 
}

#   Read in all the separate matrices
#   Markers are columns, individuals are rows
#   The first column has an integer that describes population of origin
soja.all <- read.table('32k_ATCG.txt', header=T, row.names=2)
soja.pop1 <- soja.all[soja.all$Pop == 1,]
soja.pop2 <- soja.all[soja.all$Pop == 2,]
soja.pop3 <- soja.all[soja.all$Pop == 3,]


#   Get the minor allele frequencies for all markers in all partitions
all.mafs <- apply(soja.all, 2, MAF)
pop1.mafs <- apply(soja.pop1, 2, MAF)
pop2.mafs <- apply(soja.pop2, 2, MAF)
pop3.mafs <- apply(soja.pop3, 2, MAF)

#   Then, we bin them up
bins <- seq(0, 0.5, by=0.05)
all.mafs.bins <- cut(all.mafs, breaks=bins, include.lowest=TRUE)
pop1.mafs.bins <- cut(pop1.mafs, breaks=bins, include.lowest=TRUE)
pop2.mafs.bins <- cut(pop2.mafs, breaks=bins, include.lowest=TRUE)
pop3.mafs.bins <- cut(pop3.mafs, breaks=bins, include.lowest=TRUE)

#   Count up how many are in each bin
all.mafs.SFS <- table(all.mafs.bins)/length(all.mafs)
pop1.mafs.SFS <- table(pop1.mafs.bins)/length(pop1.mafs)
pop2.mafs.SFS <- table(pop2.mafs.bins)/length(pop2.mafs)
pop3.mafs.SFS <- table(pop3.mafs.bins)/length(pop3.mafs)

#   Stick them together into a data frame so we can plot side-by-side
SFSData <- as.data.frame(cbind(all.mafs.SFS, pop1.mafs.SFS, pop2.mafs.SFS, pop3.mafs.SFS))

#   Open a handle to a PDF file for the plot
pdf(file="Soja_FoldedSFS_ByCluster.pdf", width=8, height=6, family="Helvetica", pointsize=16)
#   Make the barplot, we have to do a separate x-axis later
barplot(t(SFSData),
    ylim=c(0,0.20),
    beside=TRUE,
    axisnames=F,
    xlab="Minor Allele Frequency",
    ylab="Proportion",
    col=c("black", "green2", "red", "blue"))
labels <- c("0-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50")

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
