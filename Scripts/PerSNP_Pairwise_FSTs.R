#   load up the hierfstat library
library(hierfstat)
#   and ggplot2 for making pretty pictures
library(ggplot2)
library(grid)
library(gridExtra)
#   How big to make the windows?
windowsize <- 10
#   and how many to slide over?
step <- 5

#   Create a new function to calculate sliding averages and generate the midpoints
SlidingWindow <- function(fsts, coords, winsize, winstep)
{
    #   Create a vector of window start positions
    points <- seq(from=1, to=length(fsts)-winsize, by=winstep)
    #   Two new vectors to grow that contain the average FST and position
    fst.averages <- c()
    fst.midpoints <- c()
    for(i in 1:length(points))
    {
        if( (i+winsize) <= length(fsts))
        {
            fst.averages[i] <- mean(fsts[points[i]:(points[i]+winsize)], na.rm=T)
            #   A little risky, but I am going to assume that the FST values and the
            #   coords vector are the same size
            fst.midpoints[i] <- mean(coords[points[i]:(points[i]+winsize)])
        }
        else
        {
            #   If we are pushed outside the window, we just take until the end
            fst.averages[i] <- mean(fsts[points[i]:length(fsts)], na.rm=T)
            fst.midpoints[i] <- mean(coords[points[i]:length(fsts)])
        }
    }
    sliding_windows <- data.frame(FST=fst.averages, Coords=fst.midpoints)
    return(sliding_windows)
}
#   Function yanked from
#       https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(p){
tmp <- ggplotGrob(p)
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
return(legend)}

#   Read in all the separate matrices
#   Markers are columns, individuals are rows
#   The first column has an integer that describes population of origin
dat <- read.table("32k_Numeric.txt", header=T, row.names=2)
pop1 <- dat[dat$Pop == 1,]
pop2 <- dat[dat$Pop == 2,]
pop3 <- dat[dat$Pop == 3,]

#   Get the marker names, since we will use this to do per-chromosome FST
#   plots later
marker.names <- names(pop1)[2:ncol(pop1)]

#   Put the pairwise comparisons into separate dataframes
combined.12 <- rbind(pop1, pop2)
combined.13 <- rbind(pop1, pop3)
combined.23 <- rbind(pop2, pop3)

#   calculate those FSTs!
fsts.12 <- wc(combined.12, diploid=F)
fsts.13 <- wc(combined.13, diploid=F)
fsts.23 <- wc(combined.23, diploid=F)
#   And the total FST
fst.total <- wc(dat, diploid=F)

#   What is overall FST?
print(fsts.12$FST)
print(fsts.13$FST)
print(fsts.23$FST)

#   Save all the FST comparisons in a table
fst_table <- data.frame(
    SNPName=marker.names,
    MS_V_I=as.vector(unlist(fsts.12$per.loc)),
    I_V_MN=as.vector(unlist(fsts.23$per.loc)),
    MS_V_MN=as.vector(unlist(fsts.13$per.loc)),
    Total=as.vector(unlist(fst.total$per.loc))
    )

write.table(fst_table, file='Pairwise_WC_FST_PerSNP.txt', quote=FALSE, sep="\t", row.names=FALSE)

# #   Make some pretty plots
# pdf(file="FST_12.pdf", width=6, height=6)
# hist(unlist(fsts.12$per.loc), main="Distribution of FST per SNP Between Clusters 1 and 2", xlab="Per-SNP FST", ylab="Density", freq=F, col="lightgreen")
# dev.off()
# pdf(file="FST_13.pdf", width=6, height=6)
# hist(unlist(fsts.13$per.loc), main="Distribution of FST per SNP Between Clusters 1 and 3", xlab="Per-SNP FST", ylab="Density", freq=F, col="lightgreen")
# dev.off()
# pdf(file="FST_23.pdf", width=6, height=6)
# hist(unlist(fsts.23$per.loc), main="Distribution of FST per SNP Between Clusters 2 and 3", xlab="Per-SNP FST", ylab="Density", freq=F, col="lightgreen")
# dev.off()

# #   Next, we get the per-chromosome FST values and plot those
# #   get the SNPs on a certain chromosome
# #       The names
# chr1.names <- grep("Chr01", marker.names, value=T)
# #       And the indices
# chr1.index <- grep("Chr01", marker.names, value=F)
# #       Strip out the coordinates
# chr1.coords <- as.numeric(gsub(pattern="Chr01_", replacement="", chr1.names))
# chr2.names <- grep("Chr02", marker.names, value=T)
# chr2.index <- grep("Chr02", marker.names, value=F)
# chr2.coords <- as.numeric(gsub(pattern="Chr02_", replacement="", chr2.names))
# chr3.names <- grep("Chr03", marker.names, value=T)
# chr3.index <- grep("Chr03", marker.names, value=F)
# chr3.coords <- as.numeric(gsub(pattern="Chr03_", replacement="", chr3.names))
# chr4.names <- grep("Chr04", marker.names, value=T)
# chr4.index <- grep("Chr04", marker.names, value=F)
# chr4.coords <- as.numeric(gsub(pattern="Chr04_", replacement="", chr4.names))
# chr5.names <- grep("Chr05", marker.names, value=T)
# chr5.index <- grep("Chr05", marker.names, value=F)
# chr5.coords <- as.numeric(gsub(pattern="Chr05_", replacement="", chr5.names))
# chr6.names <- grep("Chr06", marker.names, value=T)
# chr6.index <- grep("Chr06", marker.names, value=F)
# chr6.coords <- as.numeric(gsub(pattern="Chr06_", replacement="", chr6.names))
# chr7.names <- grep("Chr07", marker.names, value=T)
# chr7.index <- grep("Chr07", marker.names, value=F)
# chr7.coords <- as.numeric(gsub(pattern="Chr07_", replacement="", chr7.names))
# chr8.names <- grep("Chr08", marker.names, value=T)
# chr8.index <- grep("Chr08", marker.names, value=F)
# chr8.coords <- as.numeric(gsub(pattern="Chr08_", replacement="", chr8.names))
# chr9.names <- grep("Chr09", marker.names, value=T)
# chr9.index <- grep("Chr09", marker.names, value=F)
# chr9.coords <- as.numeric(gsub(pattern="Chr09_", replacement="", chr9.names))
# chr10.names <- grep("Chr10", marker.names, value=T)
# chr10.index <- grep("Chr10", marker.names, value=F)
# chr10.coords <- as.numeric(gsub(pattern="Chr10_", replacement="", chr10.names))
# chr11.names <- grep("Chr11", marker.names, value=T)
# chr11.index <- grep("Chr11", marker.names, value=F)
# chr11.coords <- as.numeric(gsub(pattern="Chr11_", replacement="", chr11.names))
# chr12.names <- grep("Chr12", marker.names, value=T)
# chr12.index <- grep("Chr12", marker.names, value=F)
# chr12.coords <- as.numeric(gsub(pattern="Chr12_", replacement="", chr12.names))
# chr13.names <- grep("Chr13", marker.names, value=T)
# chr13.index <- grep("Chr13", marker.names, value=F)
# chr13.coords <- as.numeric(gsub(pattern="Chr13_", replacement="", chr13.names))
# chr14.names <- grep("Chr14", marker.names, value=T)
# chr14.index <- grep("Chr14", marker.names, value=F)
# chr14.coords <- as.numeric(gsub(pattern="Chr14_", replacement="", chr14.names))
# chr15.names <- grep("Chr15", marker.names, value=T)
# chr15.index <- grep("Chr15", marker.names, value=F)
# chr15.coords <- as.numeric(gsub(pattern="Chr15_", replacement="", chr15.names))
# chr16.names <- grep("Chr16", marker.names, value=T)
# chr16.index <- grep("Chr16", marker.names, value=F)
# chr16.coords <- as.numeric(gsub(pattern="Chr16_", replacement="", chr16.names))
# chr17.names <- grep("Chr17", marker.names, value=T)
# chr17.index <- grep("Chr17", marker.names, value=F)
# chr17.coords <- as.numeric(gsub(pattern="Chr17_", replacement="", chr17.names))
# chr18.names <- grep("Chr18", marker.names, value=T)
# chr18.index <- grep("Chr18", marker.names, value=F)
# chr18.coords <- as.numeric(gsub(pattern="Chr18_", replacement="", chr18.names))
# chr19.names <- grep("Chr19", marker.names, value=T)
# chr19.index <- grep("Chr19", marker.names, value=F)
# chr19.coords <- as.numeric(gsub(pattern="Chr19_", replacement="", chr19.names))
# chr20.names <- grep("Chr20", marker.names, value=T)
# chr20.index <- grep("Chr20", marker.names, value=F)
# chr20.coords <- as.numeric(gsub(pattern="Chr20_", replacement="", chr20.names))


# #   Next, we split up the fst values by the chromosome they are on
# #   For pops 1 and 2...
# fst.12.gm01 <- fsts.12$per.loc$FST[chr1.index]
# fst.12.gm02 <- fsts.12$per.loc$FST[chr2.index]
# fst.12.gm03 <- fsts.12$per.loc$FST[chr3.index]
# fst.12.gm04 <- fsts.12$per.loc$FST[chr4.index]
# fst.12.gm05 <- fsts.12$per.loc$FST[chr5.index]
# fst.12.gm06 <- fsts.12$per.loc$FST[chr6.index]
# fst.12.gm07 <- fsts.12$per.loc$FST[chr7.index]
# fst.12.gm08 <- fsts.12$per.loc$FST[chr8.index]
# fst.12.gm09 <- fsts.12$per.loc$FST[chr9.index]
# fst.12.gm10 <- fsts.12$per.loc$FST[chr10.index]
# fst.12.gm11 <- fsts.12$per.loc$FST[chr11.index]
# fst.12.gm12 <- fsts.12$per.loc$FST[chr12.index]
# fst.12.gm13 <- fsts.12$per.loc$FST[chr13.index]
# fst.12.gm14 <- fsts.12$per.loc$FST[chr14.index]
# fst.12.gm15 <- fsts.12$per.loc$FST[chr15.index]
# fst.12.gm16 <- fsts.12$per.loc$FST[chr16.index]
# fst.12.gm17 <- fsts.12$per.loc$FST[chr17.index]
# fst.12.gm18 <- fsts.12$per.loc$FST[chr18.index]
# fst.12.gm19 <- fsts.12$per.loc$FST[chr19.index]
# fst.12.gm20 <- fsts.12$per.loc$FST[chr20.index]
# #   Get the sliding windows!
# fst.12.gm01.windows <- SlidingWindow(fst.12.gm01, chr1.coords, windowsize, step)
# fst.12.gm02.windows <- SlidingWindow(fst.12.gm02, chr2.coords, windowsize, step)
# fst.12.gm03.windows <- SlidingWindow(fst.12.gm03, chr3.coords, windowsize, step)
# fst.12.gm04.windows <- SlidingWindow(fst.12.gm04, chr4.coords, windowsize, step)
# fst.12.gm05.windows <- SlidingWindow(fst.12.gm05, chr5.coords, windowsize, step)
# fst.12.gm06.windows <- SlidingWindow(fst.12.gm06, chr6.coords, windowsize, step)
# fst.12.gm07.windows <- SlidingWindow(fst.12.gm07, chr7.coords, windowsize, step)
# fst.12.gm08.windows <- SlidingWindow(fst.12.gm08, chr8.coords, windowsize, step)
# fst.12.gm09.windows <- SlidingWindow(fst.12.gm09, chr9.coords, windowsize, step)
# fst.12.gm10.windows <- SlidingWindow(fst.12.gm10, chr10.coords, windowsize, step)
# fst.12.gm11.windows <- SlidingWindow(fst.12.gm11, chr11.coords, windowsize, step)
# fst.12.gm12.windows <- SlidingWindow(fst.12.gm12, chr12.coords, windowsize, step)
# fst.12.gm13.windows <- SlidingWindow(fst.12.gm13, chr13.coords, windowsize, step)
# fst.12.gm14.windows <- SlidingWindow(fst.12.gm14, chr14.coords, windowsize, step)
# fst.12.gm15.windows <- SlidingWindow(fst.12.gm15, chr15.coords, windowsize, step)
# fst.12.gm16.windows <- SlidingWindow(fst.12.gm16, chr16.coords, windowsize, step)
# fst.12.gm17.windows <- SlidingWindow(fst.12.gm17, chr17.coords, windowsize, step)
# fst.12.gm18.windows <- SlidingWindow(fst.12.gm18, chr18.coords, windowsize, step)
# fst.12.gm19.windows <- SlidingWindow(fst.12.gm19, chr19.coords, windowsize, step)
# fst.12.gm20.windows <- SlidingWindow(fst.12.gm20, chr20.coords, windowsize, step)

# #   The same for the other two comparisons...
# fst.13.gm01 <- fsts.13$per.loc$FST[chr1.index]
# fst.13.gm02 <- fsts.13$per.loc$FST[chr2.index]
# fst.13.gm03 <- fsts.13$per.loc$FST[chr3.index]
# fst.13.gm04 <- fsts.13$per.loc$FST[chr4.index]
# fst.13.gm05 <- fsts.13$per.loc$FST[chr5.index]
# fst.13.gm06 <- fsts.13$per.loc$FST[chr6.index]
# fst.13.gm07 <- fsts.13$per.loc$FST[chr7.index]
# fst.13.gm08 <- fsts.13$per.loc$FST[chr8.index]
# fst.13.gm09 <- fsts.13$per.loc$FST[chr9.index]
# fst.13.gm10 <- fsts.13$per.loc$FST[chr10.index]
# fst.13.gm11 <- fsts.13$per.loc$FST[chr11.index]
# fst.13.gm12 <- fsts.13$per.loc$FST[chr12.index]
# fst.13.gm13 <- fsts.13$per.loc$FST[chr13.index]
# fst.13.gm14 <- fsts.13$per.loc$FST[chr14.index]
# fst.13.gm15 <- fsts.13$per.loc$FST[chr15.index]
# fst.13.gm16 <- fsts.13$per.loc$FST[chr16.index]
# fst.13.gm17 <- fsts.13$per.loc$FST[chr17.index]
# fst.13.gm18 <- fsts.13$per.loc$FST[chr18.index]
# fst.13.gm19 <- fsts.13$per.loc$FST[chr19.index]
# fst.13.gm20 <- fsts.13$per.loc$FST[chr20.index]
# #   Get the sliding windows!
# fst.13.gm01.windows <- SlidingWindow(fst.13.gm01, chr1.coords, windowsize, step)
# fst.13.gm02.windows <- SlidingWindow(fst.13.gm02, chr2.coords, windowsize, step)
# fst.13.gm03.windows <- SlidingWindow(fst.13.gm03, chr3.coords, windowsize, step)
# fst.13.gm04.windows <- SlidingWindow(fst.13.gm04, chr4.coords, windowsize, step)
# fst.13.gm05.windows <- SlidingWindow(fst.13.gm05, chr5.coords, windowsize, step)
# fst.13.gm06.windows <- SlidingWindow(fst.13.gm06, chr6.coords, windowsize, step)
# fst.13.gm07.windows <- SlidingWindow(fst.13.gm07, chr7.coords, windowsize, step)
# fst.13.gm08.windows <- SlidingWindow(fst.13.gm08, chr8.coords, windowsize, step)
# fst.13.gm09.windows <- SlidingWindow(fst.13.gm09, chr9.coords, windowsize, step)
# fst.13.gm10.windows <- SlidingWindow(fst.13.gm10, chr10.coords, windowsize, step)
# fst.13.gm11.windows <- SlidingWindow(fst.13.gm11, chr11.coords, windowsize, step)
# fst.13.gm12.windows <- SlidingWindow(fst.13.gm12, chr12.coords, windowsize, step)
# fst.13.gm13.windows <- SlidingWindow(fst.13.gm13, chr13.coords, windowsize, step)
# fst.13.gm14.windows <- SlidingWindow(fst.13.gm14, chr14.coords, windowsize, step)
# fst.13.gm15.windows <- SlidingWindow(fst.13.gm15, chr15.coords, windowsize, step)
# fst.13.gm16.windows <- SlidingWindow(fst.13.gm16, chr16.coords, windowsize, step)
# fst.13.gm17.windows <- SlidingWindow(fst.13.gm17, chr17.coords, windowsize, step)
# fst.13.gm18.windows <- SlidingWindow(fst.13.gm18, chr18.coords, windowsize, step)
# fst.13.gm19.windows <- SlidingWindow(fst.13.gm19, chr19.coords, windowsize, step)
# fst.13.gm20.windows <- SlidingWindow(fst.13.gm20, chr20.coords, windowsize, step)

# fst.23.gm01 <- fsts.23$per.loc$FST[chr1.index]
# fst.23.gm02 <- fsts.23$per.loc$FST[chr2.index]
# fst.23.gm03 <- fsts.23$per.loc$FST[chr3.index]
# fst.23.gm04 <- fsts.23$per.loc$FST[chr4.index]
# fst.23.gm05 <- fsts.23$per.loc$FST[chr5.index]
# fst.23.gm06 <- fsts.23$per.loc$FST[chr6.index]
# fst.23.gm07 <- fsts.23$per.loc$FST[chr7.index]
# fst.23.gm08 <- fsts.23$per.loc$FST[chr8.index]
# fst.23.gm09 <- fsts.23$per.loc$FST[chr9.index]
# fst.23.gm10 <- fsts.23$per.loc$FST[chr10.index]
# fst.23.gm11 <- fsts.23$per.loc$FST[chr11.index]
# fst.23.gm12 <- fsts.23$per.loc$FST[chr12.index]
# fst.23.gm13 <- fsts.23$per.loc$FST[chr13.index]
# fst.23.gm14 <- fsts.23$per.loc$FST[chr14.index]
# fst.23.gm15 <- fsts.23$per.loc$FST[chr15.index]
# fst.23.gm16 <- fsts.23$per.loc$FST[chr16.index]
# fst.23.gm17 <- fsts.23$per.loc$FST[chr17.index]
# fst.23.gm18 <- fsts.23$per.loc$FST[chr18.index]
# fst.23.gm19 <- fsts.23$per.loc$FST[chr19.index]
# fst.23.gm20 <- fsts.23$per.loc$FST[chr20.index]
# #   Get the sliding windows!
# fst.23.gm01.windows <- SlidingWindow(fst.23.gm01, chr1.coords, windowsize, step)
# fst.23.gm02.windows <- SlidingWindow(fst.23.gm02, chr2.coords, windowsize, step)
# fst.23.gm03.windows <- SlidingWindow(fst.23.gm03, chr3.coords, windowsize, step)
# fst.23.gm04.windows <- SlidingWindow(fst.23.gm04, chr4.coords, windowsize, step)
# fst.23.gm05.windows <- SlidingWindow(fst.23.gm05, chr5.coords, windowsize, step)
# fst.23.gm06.windows <- SlidingWindow(fst.23.gm06, chr6.coords, windowsize, step)
# fst.23.gm07.windows <- SlidingWindow(fst.23.gm07, chr7.coords, windowsize, step)
# fst.23.gm08.windows <- SlidingWindow(fst.23.gm08, chr8.coords, windowsize, step)
# fst.23.gm09.windows <- SlidingWindow(fst.23.gm09, chr9.coords, windowsize, step)
# fst.23.gm10.windows <- SlidingWindow(fst.23.gm10, chr10.coords, windowsize, step)
# fst.23.gm11.windows <- SlidingWindow(fst.23.gm11, chr11.coords, windowsize, step)
# fst.23.gm12.windows <- SlidingWindow(fst.23.gm12, chr12.coords, windowsize, step)
# fst.23.gm13.windows <- SlidingWindow(fst.23.gm13, chr13.coords, windowsize, step)
# fst.23.gm14.windows <- SlidingWindow(fst.23.gm14, chr14.coords, windowsize, step)
# fst.23.gm15.windows <- SlidingWindow(fst.23.gm15, chr15.coords, windowsize, step)
# fst.23.gm16.windows <- SlidingWindow(fst.23.gm16, chr16.coords, windowsize, step)
# fst.23.gm17.windows <- SlidingWindow(fst.23.gm17, chr17.coords, windowsize, step)
# fst.23.gm18.windows <- SlidingWindow(fst.23.gm18, chr18.coords, windowsize, step)
# fst.23.gm19.windows <- SlidingWindow(fst.23.gm19, chr19.coords, windowsize, step)
# fst.23.gm20.windows <- SlidingWindow(fst.23.gm20, chr20.coords, windowsize, step)

# #   Build the data frames to plot
# fst.gm01.total <- data.frame(
#     FST=c(fst.12.gm01.windows$FST, fst.13.gm01.windows$FST, fst.23.gm01.windows$FST),
#     pos=c(fst.12.gm01.windows$Coords, fst.13.gm01.windows$Coords, fst.23.gm01.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm01.windows$FST)), rep("1-3", length(fst.13.gm01.windows$FST)), rep("2-3", length(fst.23.gm01.windows$FST))),
#     p_start=5320000,
#     p_end=47200000
# )
# fst.gm02.total <- data.frame(
#     FST=c(fst.12.gm02.windows$FST, fst.13.gm02.windows$FST, fst.23.gm02.windows$FST),
#     pos=c(fst.12.gm02.windows$Coords, fst.13.gm02.windows$Coords, fst.23.gm02.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm02.windows$FST)), rep("1-3", length(fst.13.gm02.windows$FST)), rep("2-3", length(fst.23.gm02.windows$FST))),
#     p_start=16900000,
#     p_end=38100000
# )
# fst.gm03.total <- data.frame(
#     FST=c(fst.12.gm03.windows$FST, fst.13.gm03.windows$FST, fst.23.gm03.windows$FST),
#     pos=c(fst.12.gm03.windows$Coords, fst.13.gm03.windows$Coords, fst.23.gm03.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm03.windows$FST)), rep("1-3", length(fst.13.gm03.windows$FST)), rep("2-3", length(fst.23.gm03.windows$FST))),
#     p_start=5480000,
#     p_end=31600000
# )
# fst.gm04.total <- data.frame(
#     FST=c(fst.12.gm04.windows$FST, fst.13.gm04.windows$FST, fst.23.gm04.windows$FST),
#     pos=c(fst.12.gm04.windows$Coords, fst.13.gm04.windows$Coords, fst.23.gm04.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm04.windows$FST)), rep("1-3", length(fst.13.gm04.windows$FST)), rep("2-3", length(fst.23.gm04.windows$FST))),
#     p_start=9650000,
#     p_end=42800000
# )
# fst.gm05.total <- data.frame(
#     FST=c(fst.12.gm05.windows$FST, fst.13.gm05.windows$FST, fst.23.gm05.windows$FST),
#     pos=c(fst.12.gm05.windows$Coords, fst.13.gm05.windows$Coords, fst.23.gm05.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm05.windows$FST)), rep("1-3", length(fst.13.gm05.windows$FST)), rep("2-3", length(fst.23.gm05.windows$FST))),
#     p_start=7450000,
#     p_end=27000000
# )
# fst.gm06.total <- data.frame(
#     FST=c(fst.12.gm06.windows$FST, fst.13.gm06.windows$FST, fst.23.gm06.windows$FST),
#     pos=c(fst.12.gm06.windows$Coords, fst.13.gm06.windows$Coords, fst.23.gm06.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm06.windows$FST)), rep("1-3", length(fst.13.gm06.windows$FST)), rep("2-3", length(fst.23.gm06.windows$FST))),
#     p_start=17900000,
#     p_end=47300000
# )
# fst.gm07.total <- data.frame(
#     FST=c(fst.12.gm07.windows$FST, fst.13.gm07.windows$FST, fst.23.gm07.windows$FST),
#     pos=c(fst.12.gm07.windows$Coords, fst.13.gm07.windows$Coords, fst.23.gm07.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm07.windows$FST)), rep("1-3", length(fst.13.gm07.windows$FST)), rep("2-3", length(fst.23.gm07.windows$FST))),
#     p_start=18000000,
#     p_end=34700000
# )
# fst.gm08.total <- data.frame(
#     FST=c(fst.12.gm08.windows$FST, fst.13.gm08.windows$FST, fst.23.gm08.windows$FST),
#     pos=c(fst.12.gm08.windows$Coords, fst.13.gm08.windows$Coords, fst.23.gm08.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm08.windows$FST)), rep("1-3", length(fst.13.gm08.windows$FST)), rep("2-3", length(fst.23.gm08.windows$FST))),
#     p_start=23300000,
#     p_end=39900000
# )
# fst.gm09.total <- data.frame(
#     FST=c(fst.12.gm09.windows$FST, fst.13.gm09.windows$FST, fst.23.gm09.windows$FST),
#     pos=c(fst.12.gm09.windows$Coords, fst.13.gm09.windows$Coords, fst.23.gm09.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm09.windows$FST)), rep("1-3", length(fst.13.gm09.windows$FST)), rep("2-3", length(fst.23.gm09.windows$FST))),
#     p_start=6330000,
#     p_end=37300000
# )
# fst.gm10.total <- data.frame(
#     FST=c(fst.12.gm10.windows$FST, fst.13.gm10.windows$FST, fst.23.gm10.windows$FST),
#     pos=c(fst.12.gm10.windows$Coords, fst.13.gm10.windows$Coords, fst.23.gm10.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm10.windows$FST)), rep("1-3", length(fst.13.gm10.windows$FST)), rep("2-3", length(fst.23.gm10.windows$FST))),
#     p_start=6950000,
#     p_end=36800000
# )
# fst.gm11.total <- data.frame(
#     FST=c(fst.12.gm11.windows$FST, fst.13.gm11.windows$FST, fst.23.gm11.windows$FST),
#     pos=c(fst.12.gm11.windows$Coords, fst.13.gm11.windows$Coords, fst.23.gm11.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm11.windows$FST)), rep("1-3", length(fst.13.gm11.windows$FST)), rep("2-3", length(fst.23.gm11.windows$FST))),
#     p_start=13000000,
#     p_end=29000000
# )
# fst.gm12.total <- data.frame(
#     FST=c(fst.12.gm12.windows$FST, fst.13.gm12.windows$FST, fst.23.gm12.windows$FST),
#     pos=c(fst.12.gm12.windows$Coords, fst.13.gm12.windows$Coords, fst.23.gm12.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm12.windows$FST)), rep("1-3", length(fst.13.gm12.windows$FST)), rep("2-3", length(fst.23.gm12.windows$FST))),
#     p_start=9300000,
#     p_end=32700000
# )
# fst.gm13.total <- data.frame(
#     FST=c(fst.12.gm13.windows$FST, fst.13.gm13.windows$FST, fst.23.gm13.windows$FST),
#     pos=c(fst.12.gm13.windows$Coords, fst.13.gm13.windows$Coords, fst.23.gm13.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm13.windows$FST)), rep("1-3", length(fst.13.gm13.windows$FST)), rep("2-3", length(fst.23.gm13.windows$FST))),
#     p_start=0,
#     p_end=9390000
# )
# fst.gm14.total <- data.frame(
#     FST=c(fst.12.gm14.windows$FST, fst.13.gm14.windows$FST, fst.23.gm14.windows$FST),
#     pos=c(fst.12.gm14.windows$Coords, fst.13.gm14.windows$Coords, fst.23.gm14.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm14.windows$FST)), rep("1-3", length(fst.13.gm14.windows$FST)), rep("2-3", length(fst.23.gm14.windows$FST))),
#     p_start=9830000,
#     p_end=43600000
# )
# fst.gm15.total <- data.frame(
#     FST=c(fst.12.gm15.windows$FST, fst.13.gm15.windows$FST, fst.23.gm15.windows$FST),
#     pos=c(fst.12.gm15.windows$Coords, fst.13.gm15.windows$Coords, fst.23.gm15.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm15.windows$FST)), rep("1-3", length(fst.13.gm15.windows$FST)), rep("2-3", length(fst.23.gm15.windows$FST))),
#     p_start=15800000,
#     p_end=47100000
# )
# fst.gm16.total <- data.frame(
#     FST=c(fst.12.gm16.windows$FST, fst.13.gm16.windows$FST, fst.23.gm16.windows$FST),
#     pos=c(fst.12.gm16.windows$Coords, fst.13.gm16.windows$Coords, fst.23.gm16.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm16.windows$FST)), rep("1-3", length(fst.13.gm16.windows$FST)), rep("2-3", length(fst.23.gm16.windows$FST))),
#     p_start=7330000,
#     p_end=27800000
# )
# fst.gm17.total <- data.frame(
#     FST=c(fst.12.gm17.windows$FST, fst.13.gm17.windows$FST, fst.23.gm17.windows$FST),
#     pos=c(fst.12.gm17.windows$Coords, fst.13.gm17.windows$Coords, fst.23.gm17.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm17.windows$FST)), rep("1-3", length(fst.13.gm17.windows$FST)), rep("2-3", length(fst.23.gm17.windows$FST))),
#     p_start=14100000,
#     p_end=37300000
# )
# fst.gm18.total <- data.frame(
#     FST=c(fst.12.gm18.windows$FST, fst.13.gm18.windows$FST, fst.23.gm18.windows$FST),
#     pos=c(fst.12.gm18.windows$Coords, fst.13.gm18.windows$Coords, fst.23.gm18.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm18.windows$FST)), rep("1-3", length(fst.13.gm18.windows$FST)), rep("2-3", length(fst.23.gm18.windows$FST))),
#     p_start=17200000,
#     p_end=42700000
# )
# fst.gm19.total <- data.frame(
#     FST=c(fst.12.gm19.windows$FST, fst.13.gm19.windows$FST, fst.23.gm19.windows$FST),
#     pos=c(fst.12.gm19.windows$Coords, fst.13.gm19.windows$Coords, fst.23.gm19.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm19.windows$FST)), rep("1-3", length(fst.13.gm19.windows$FST)), rep("2-3", length(fst.23.gm19.windows$FST))),
#     p_start=4560000,
#     p_end=34300000
# )
# fst.gm20.total <- data.frame(
#     FST=c(fst.12.gm20.windows$FST, fst.13.gm20.windows$FST, fst.23.gm20.windows$FST),
#     pos=c(fst.12.gm20.windows$Coords, fst.13.gm20.windows$Coords, fst.23.gm20.windows$Coords),
#     lab=c(rep("1-2", length(fst.12.gm20.windows$FST)), rep("1-3", length(fst.13.gm20.windows$FST)), rep("2-3", length(fst.23.gm20.windows$FST))),
#     p_start=5390000,
#     p_end=32400000
# )

# #   Start building plots!
# p1 <- ggplot(
#     fst.gm01.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm01", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p2 <- ggplot(
#     fst.gm02.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm02", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p3 <- ggplot(
#     fst.gm03.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm03", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p4 <- ggplot(
#     fst.gm04.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm04", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p5 <- ggplot(
#     fst.gm05.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm05", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p6 <- ggplot(
#     fst.gm06.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm06", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p7 <- ggplot(
#     fst.gm07.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm07", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p8 <- ggplot(
#     fst.gm08.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm08", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p9 <- ggplot(
#     fst.gm09.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm09", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p10 <- ggplot(
#     fst.gm10.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm10", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p11 <- ggplot(
#     fst.gm11.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm11", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p12 <- ggplot(
#     fst.gm12.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm12", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p13 <- ggplot(
#     fst.gm13.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm13", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p14 <- ggplot(
#     fst.gm14.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm14", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p15 <- ggplot(
#     fst.gm15.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm15", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p16 <- ggplot(
#     fst.gm16.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm16", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p17 <- ggplot(
#     fst.gm17.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm17", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p18 <- ggplot(
#     fst.gm18.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm18", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p19 <- ggplot(
#     fst.gm19.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm19", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))
# p20 <- ggplot(
#     fst.gm20.total,
#     aes(x=pos, y=FST, color=lab, ymin=0, ymax=1, xmin=p_start, xmax=p_end)) +
#     geom_rect(aes(alpha="Pericentromere"), fill="grey90", color="grey90") +
#     geom_line(size=0.3) +
#     theme_bw() +
#     theme(axis.text=element_text(size=8), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
#     labs(title="Gm20", x="", y="") +
#     coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
#     scale_x_continuous(breaks=seq(0,60000000, 5000000)) +
#     scale_y_continuous(breaks=seq(0, 1, 0.2))

# #   Arrange and plot them!

# legend <- g_legend(p1)
# lwidth <- sum(legend$width)

# pdf(file="GM01-05.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p1 + theme(legend.position="none"),
#         p2 + theme(legend.position="none"),
#         p3 + theme(legend.position="none"),
#         p4 + theme(legend.position="none"),
#         p5 + theme(legend.position="none"),
#         main="FST Between Clusters",
#         left="FST",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
# pdf(file="GM06-10.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p6 + theme(legend.position="none"),
#         p7 + theme(legend.position="none"),
#         p8 + theme(legend.position="none"),
#         p9 + theme(legend.position="none"),
#         p10 + theme(legend.position="none"),
#         main="FST Between Clusters",
#         left="FST",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
# pdf(file="GM11-15.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p11 + theme(legend.position="none"),
#         p12 + theme(legend.position="none"),
#         p13 + theme(legend.position="none"),
#         p14 + theme(legend.position="none"),
#         p15 + theme(legend.position="none"),
#         main="FST Between Clusters",
#         left="FST",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
# pdf(file="GM16-20.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p16 + theme(legend.position="none"),
#         p17 + theme(legend.position="none"),
#         p18 + theme(legend.position="none"),
#         p19 + theme(legend.position="none"),
#         p20 + theme(legend.position="none"),
#         main="FST Between Clusters",
#         left="FST",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
