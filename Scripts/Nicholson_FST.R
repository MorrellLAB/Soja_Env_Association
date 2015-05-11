#   Load the R package to do the estimate
library(popgen)
#   To do density plots
library(ggplot2)
#  Take the arguments
args <- commandArgs(TRUE)

#   Set the arguments to the correct positional args
numalleles_file <- args[1]
numgenotypes_file <- args[2]
popalleles_file <- args[3]

#   We read in the relevant data
#   We will set these to the variables that are defined in the popgen
#   pacakge. These are not nice names, but I am sacrificing that in favor
#   of consistency with the package docs.
#       NUMA: Number of alleles per locus
#       N: Number of genotypes per population (non-missing individuals)
#       X: Count of each allelic class in each population
NUMA <- read.table(numalleles_file, header=F)
#   Number of alleles per locus should be a vector
NUMA <- as.vector(NUMA$V1)
#   The other two should be matrices
N <- read.table(numgenotypes_file, header=F)
N <- as.matrix(N)
X <- read.table(popalleles_file, header=F)
X <- as.matrix(X)

#   Set the parameters of our run here
#       L: number of loci
#       P: number of pops
#       burnin: length of burn-in iterations for MCMC
#       iter: number of sampled iterations for MCMC
#       csd: sd of the proposal distribution
#       outc: whether or not to return samples of c (FST)
L <- length(NUMA)
P <- ncol(N)
#   Set these values high so we can get a "hot" chain
csd <- 0.1
m <- 600
burnin <- 5000
iter <- 5000
outc <- TRUE

#   Run the estimator
fst <- popdiv(
    L,
    P,
    NUMA,
    N,
    X,
    m=m,
    csd=csd,
    burnin=burnin,
    iter=iter,
    outc=outc
    )

#   Plot the MCMC chain
pdf(file="Nicholson_FST_MCMC.pdf", 8, 8)
#   This creates an empty plot with the dimensions we set above
plot(c(0, iter), c(0, 1), type="n", xlab="MCMC Step", ylab="FST")
#   Add lines
lines(1:iter, fst$C.sample[,1], col="green")
lines(1:iter, fst$C.sample[,2], col="red")
lines(1:iter, fst$C.sample[,3], col="blue")
legend(
    "topright",
    c(
        "Mainland South",
        "Island",
        "Mainland North"),
    fill=c("green", "red", "blue"))
dev.off()
#   Save the chain data, so we can plot it later
#   We throw in a 0 and 0.4 just to extend the densities out to the limits.
#   This affects the chain a bit, but hopefully 2/15000 isn't a big deal
nicholson_fst <- data.frame(
    MainlandSouth=c(0, fst$C.sample[,1], 0.4),
    Island=c(0, fst$C.sample[,2], 0.4),
    MainlandNorth=c(0, fst$C.sample[,3], 0.4)
    )
summaries <- data.frame(
  Pops=c("Mainland South", "Island", "Mainland North"),
  Means=fst$muc,
  StdDevs=fst$sdc
  )
#   Write the data later, so we don't have to re-run the whole simulation to
#   change the plot.
write.table(
    nicholson_fst,
    file="FST_MCMC.txt",
    sep="\t",
    col.names=TRUE,
    row.names=FALSE)
write.table(
    summaries,
    file="MCMC_Summaries.txt",
    sep="\t",
    col.names=TRUE,
    row.names=FALSE)

#   Plot it out!
#   Create a dataframe for plotting things with ggplot2
plot_data <- data.frame(
    FST=c(
            c(0, fst$C.sample[,1], 0.4),
            c(0, fst$C.sample[,2], 0.4),
            c(0, fst$C.sample[,3], 0.4)
        ),
    lab=c(
        rep("Mainland South", 2 + iter),
        rep("Island", 2 + iter),
        rep("Mainland North", 2 + iter))
  )

#   Just in case you wanted to see it, print out the means and standard
#   deviations of the FST estimates for each population.
print(fst$muc)
print(fst$sdc)

#   Uncomment this (and comment everything above) if you have already run the
#   estimator and saved the chain data and just want to re-make the plots.
# fsts <- read.table("FST_MCMC.txt", header=T)
# summaries <- read.table("MCMC_Summaries.txt", header=T)
# plot_data <- data.frame(
#   FST=c(fsts$MainlandSouth, fsts$Island, fsts$MainlandNorth),
#   lab=c(
#     rep("Mainland South", length(fsts$MainlandSouth)),
#     rep("Island", length(fsts$Island)),
#     rep("Mainland North", length(fsts$MainlandNorth)))
#   )

#   Create the plot object!
plt <- ggplot(
    plot_data,
    aes(x=FST, color=lab)) +
    geom_density(size=0.6) +
    #   Plot the mean of the Mainland South Cluster
    geom_vline(
        aes(xintercept=summaries$Means[1]),
        color="green2",
        linetype="solid",
        size=0.3,
        alpha=0.6) +
    #   And +/- 2 * the stddev
    geom_vline(
        aes(xintercept=summaries$Means[1] - 2 * summaries$StdDevs[1]),
        color="green2",
        linetype="dotted",
        size=0.2) +
    geom_vline(
        aes(xintercept=summaries$Means[1] + 2 * summaries$StdDevs[1]),
        color="green2",
        linetype="dotted",
        size=0.2) +
    #   The same for the Island Cluster
    geom_vline(
        aes(xintercept=summaries$Means[2]),
        color="red",
        linetype="solid",
        size=0.3,
        alpha=0.6) +
    geom_vline(
        aes(xintercept=summaries$Means[2] - 2 * summaries$StdDevs[2]),
        color="red",
        linetype="dotted",
        size=0.2) +
    geom_vline(
        aes(xintercept=summaries$Means[2] + 2 * summaries$StdDevs[2]),
        color="red",
        linetype="dotted",
        size=0.2) +
    #   And finally, the Mainland North Cluster
    geom_vline(
        aes(xintercept=summaries$Means[3]),
        color="blue",
        linetype="solid",
        size=0.3,
        alpha=0.6) +
    geom_vline(
        aes(xintercept=summaries$Means[3] - 2 * summaries$StdDevs[3]),
        color="blue",
        linetype="dotted",
        size=0.2) +
    geom_vline(
        aes(xintercept=summaries$Means[3] + 2 * summaries$StdDevs[3]),
        color="blue",
        linetype="dotted",
        size=0.2) +
    #   Set the theme variables to give us a white background with a blank
    #   x-axis label.
    theme(
        axis.text=element_text(size=12, colour="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        legend.title=element_blank(),
        axis.ticks=element_line(colour="black"),
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none") +
    #   Set the y axis label and title
    labs(title="", y="Density") +
    #   Give a proper label for FST: Italic F, subscript ST
    xlab(expression(paste(italic(F)[ST]))) +
    #   Set the colors. Create a named vector so that the population labels
    #   are correctly assigned to their colors.
    scale_color_manual(
        values=c(
            "Mainland South"="green2",
            "Island"="red",
            "Mainland North"="blue")
        )

#   And finally, write the plot into a file.
pdf("Nicholson_FST_Densities.pdf", 8, 8)
plt
dev.off()
