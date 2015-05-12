#   Script to build density plots of the per-SNP FST in pericentromeric
#   regions and euchromatic regions.
#   Written by Thomas JY Kono

#   Read in the data with the pericentromere boundaries
pericentromeres <<- read.table("Pericentromeres.txt", header=TRUE)
#   Read in the per-SNP FST values
FSTs <- read.table("FSTs_With_Pos.txt", header=TRUE)

#   A function to ask whether or not a marker is pericentromeric
centromeric <- function(x) {
    chromosome <- x[1]
    position <- x[2]
    peri_start <- pericentromeres[pericentromeres$Chromosome == chromosome, 2]
    peri_end <- pericentromeres[pericentromeres$Chromosome == chromosome, 3]
    if( position > peri_start || position < peri_end) {
        return(TRUE)
    }
    else {
        return(FALSE)
    }
    }

#   Then, we get all the markers that are pericentromeric
peri_markers_bool <- apply(FSTs, 1, centromeric)
#   These are pericentromeric
peri_markers <- FSTs[peri_markers_bool,]
#   And these are euchromatic
euchr_markers <- FSTs[!peri_markers_bool,]

#   Summarize each partition
summary(peri_markers$Total)
summary(peri_markers$MS_V_I)
summary(peri_markers$I_V_MN)
summary(peri_markers$MS_V_MN)
summary(euchr_markers$Total)
summary(euchr_markers$MS_V_I)
summary(euchr_markers$I_V_MN)
summary(euchr_markers$MS_V_MN)

#   Write the euchromatic densities into a plot
pdf(file="Euchromatic_FST_Densities.pdf", 8, 8)
#   These are the limits of the data for the plot
x <- c(0, 1)
y <- c(0, 14)
#   Put an empty plot down
plot(
    x,
    y,
    type="n",
    xlab="FST",
    ylab="Density",
    main="Density of FST in Euchromatic Regions, by Population")
#   And then add density curves for each partition
lines(density(euchr_markers$Total, na.rm=TRUE), lwd=2, lty=1, col="black")
lines(density(euchr_markers$MS_V_I, na.rm=TRUE), lwd=2, lty=2, col="red")
lines(density(euchr_markers$I_V_MN, na.rm=TRUE), lwd=2, lty=3, col="blue")
lines(density(euchr_markers$MS_V_MN, na.rm=TRUE), lwd=2, lty=4, col="green")
#   Place a legend
legend(
    0.5,
    14,
    c(
        "Total",
        "Mainland South V. Island",
        "Island V. Mainland North",
        "Mainland South V. Mainland North"),
    lty=c(1, 2, 3, 4),
    lwd=c(2, 2, 2, 2),
    col=c("black", "red", "blue", "green"))
dev.off()

#   Do the same for pericentromeric densities
pdf(file="Pericentromeric_FST_Densities.pdf", 8, 8)
plot(
    x,
    y,
    type="n",
    xlab="FST",
    ylab="Density",
    main="Density of FST in Pericentromeric Regions, by Population")
lines(density(peri_markers$Total, na.rm=TRUE), lwd=2, lty=1, col="black")
lines(density(peri_markers$MS_V_I, na.rm=TRUE), lwd=2, lty=2, col="red")
lines(density(peri_markers$I_V_MN, na.rm=TRUE), lwd=2, lty=3, col="blue")
lines(density(peri_markers$MS_V_MN, na.rm=TRUE), lwd=2, lty=4, col="green")
legend(
    0.5,
    14,
    c(
        "Total",
        "Mainland South V. Island",
        "Island V. Mainland North",
        "Mainland South V. Mainland North"),
    lty=c(1, 2, 3, 4),
    lwd=c(2, 2, 2, 2),
    col=c("black", "red", "blue", "green"))
dev.off()

#   For each population, plot the densities for pericentromeric and euchromatic
pdf(file="All_FST_Peri_v_Euch.pdf", 8, 8)
plot(
    density(peri_markers$Total, na.rm=TRUE),
    lwd=2,
    lty=1,
    xlab="FST",
    ylab="Density",
    main=paste(
        "Density of Pericetromeric FST Values and Euchromatic FST Values",
        "All Populations", sep="\n"))
lines(
    density(euchr_markers$Total, na.rm=TRUE),
    lwd=2,
    lty=2,
    col="red")
legend(
    0.5,
    7,
    c("Pericentromeric", "Euchromatic"),
    lty=c(1, 2),
    lwd=c(2, 2),
    col=c("black", "red"))
dev.off()

pdf(file="MSvI_FST_Peri_v_Euch.pdf", 8, 8)
plot(
    density(peri_markers$MS_V_I, na.rm=TRUE),
    lwd=2,
    lty=1,
    xlab="FST",
    ylab="Density",
    main=paste(
        "Density of Pericetromeric FST Values and Euchromatic FST Values",
        "Mainland South V. Island", sep="\n"))
lines(
    density(euchr_markers$MS_V_I, na.rm=TRUE),
    lwd=2,
    lty=2,
    col="red")
legend(
    0.5,
    12,
    c("Pericentromeric", "Euchromatic"),
    lty=c(1, 2),
    lwd=c(2, 2),
    col=c("black", "red"))
dev.off()

pdf(file="IvMN_FST_Peri_v_Euch.pdf", 8, 8)
plot(
    density(peri_markers$I_V_MN, na.rm=TRUE),
    lwd=2,
    lty=1,
    xlab="FST",
    ylab="Density",
    main=pate("Density of Pericetromeric FST Values and Euchromatic FST Values",
        "Island V. Mainland North", sep="\n"))
lines(
    density(euchr_markers$I_V_MN, na.rm=TRUE),
    lwd=2,
    lty=2,
    col="red")
legend(
    0.6,
    4,
    c("Pericentromeric", "Euchromatic"),
    lty=c(1, 2),
    lwd=c(2, 2),
    col=c("black", "red"))
dev.off()

pdf(file="MSvMN_FST_Peri_v_Euch.pdf", 8, 8)
plot(
    density(peri_markers$MS_V_MN, na.rm=TRUE),
    lwd=2,
    lty=1,
    xlab="FST",
    ylab="Density",
    main=pate("Density of Pericetromeric FST Values and Euchromatic FST Values",
        "Mainland South V. Mainland North", sep="\n"))
lines(
    density(euchr_markers$MS_V_MN, na.rm=TRUE),
    lwd=2,
    lty=2,
    col="red")
legend(
    0.5,
    7,
    c("Pericentromeric", "Euchromatic"),
    lty=c(1, 2),
    lwd=c(2, 2),
    col=c("black", "red"))
dev.off()
