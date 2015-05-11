#   load up the hierfstat library
library(hierfstat)
#   and ggplot2 for making pretty pictures
library(ggplot2)
library(grid)
library(gridExtra)
#   How big to make the windows?
windowsize <- 5
#   and how many to slide over?
step <- 3

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
soja <- read.table('32k_Numeric.txt', header=T, row.names=2)
#   Read in the SPA scores
spa <- read.table('SPA_Scores.txt', header=T)
#   Read the cM/Mb data
recomb <- read.table('6K_cM_per_Mb.txt', header=T)
#   There are some colinearity problems, so we cannot just take the max of this
#   naively... Examining it, most of them seem to be below 15
#maxrecomb <- max(recomb$cM_per_Mb)
maxrecomb <- 15
#   Get the max SPA score
maxspa <- max(spa$Selction_score)

#   Calculate those FSTs
fsts <- wc(soja, diploid=FALSE)

#   Get the marker names, since we will use this to do per-chromosome FST
#   plots later
marker.names <- names(soja)[2:ncol(soja)]

#   Next, we get the per-chromosome FST values and plot those
#   get the SNPs on a certain chromosome
#       The names
chr1.names <- grep("Chr01", marker.names, value=T)
#       And the indices
chr1.index <- grep("Chr01", marker.names, value=F)
#       Strip out the coordinates
chr1.coords <- as.numeric(gsub(pattern="Chr01_", replacement="", chr1.names))
chr2.names <- grep("Chr02", marker.names, value=T)
chr2.index <- grep("Chr02", marker.names, value=F)
chr2.coords <- as.numeric(gsub(pattern="Chr02_", replacement="", chr2.names))
chr3.names <- grep("Chr03", marker.names, value=T)
chr3.index <- grep("Chr03", marker.names, value=F)
chr3.coords <- as.numeric(gsub(pattern="Chr03_", replacement="", chr3.names))
chr4.names <- grep("Chr04", marker.names, value=T)
chr4.index <- grep("Chr04", marker.names, value=F)
chr4.coords <- as.numeric(gsub(pattern="Chr04_", replacement="", chr4.names))
chr5.names <- grep("Chr05", marker.names, value=T)
chr5.index <- grep("Chr05", marker.names, value=F)
chr5.coords <- as.numeric(gsub(pattern="Chr05_", replacement="", chr5.names))
chr6.names <- grep("Chr06", marker.names, value=T)
chr6.index <- grep("Chr06", marker.names, value=F)
chr6.coords <- as.numeric(gsub(pattern="Chr06_", replacement="", chr6.names))
chr7.names <- grep("Chr07", marker.names, value=T)
chr7.index <- grep("Chr07", marker.names, value=F)
chr7.coords <- as.numeric(gsub(pattern="Chr07_", replacement="", chr7.names))
chr8.names <- grep("Chr08", marker.names, value=T)
chr8.index <- grep("Chr08", marker.names, value=F)
chr8.coords <- as.numeric(gsub(pattern="Chr08_", replacement="", chr8.names))
chr9.names <- grep("Chr09", marker.names, value=T)
chr9.index <- grep("Chr09", marker.names, value=F)
chr9.coords <- as.numeric(gsub(pattern="Chr09_", replacement="", chr9.names))
chr10.names <- grep("Chr10", marker.names, value=T)
chr10.index <- grep("Chr10", marker.names, value=F)
chr10.coords <- as.numeric(gsub(pattern="Chr10_", replacement="", chr10.names))
chr11.names <- grep("Chr11", marker.names, value=T)
chr11.index <- grep("Chr11", marker.names, value=F)
chr11.coords <- as.numeric(gsub(pattern="Chr11_", replacement="", chr11.names))
chr12.names <- grep("Chr12", marker.names, value=T)
chr12.index <- grep("Chr12", marker.names, value=F)
chr12.coords <- as.numeric(gsub(pattern="Chr12_", replacement="", chr12.names))
chr13.names <- grep("Chr13", marker.names, value=T)
chr13.index <- grep("Chr13", marker.names, value=F)
chr13.coords <- as.numeric(gsub(pattern="Chr13_", replacement="", chr13.names))
chr14.names <- grep("Chr14", marker.names, value=T)
chr14.index <- grep("Chr14", marker.names, value=F)
chr14.coords <- as.numeric(gsub(pattern="Chr14_", replacement="", chr14.names))
chr15.names <- grep("Chr15", marker.names, value=T)
chr15.index <- grep("Chr15", marker.names, value=F)
chr15.coords <- as.numeric(gsub(pattern="Chr15_", replacement="", chr15.names))
chr16.names <- grep("Chr16", marker.names, value=T)
chr16.index <- grep("Chr16", marker.names, value=F)
chr16.coords <- as.numeric(gsub(pattern="Chr16_", replacement="", chr16.names))
chr17.names <- grep("Chr17", marker.names, value=T)
chr17.index <- grep("Chr17", marker.names, value=F)
chr17.coords <- as.numeric(gsub(pattern="Chr17_", replacement="", chr17.names))
chr18.names <- grep("Chr18", marker.names, value=T)
chr18.index <- grep("Chr18", marker.names, value=F)
chr18.coords <- as.numeric(gsub(pattern="Chr18_", replacement="", chr18.names))
chr19.names <- grep("Chr19", marker.names, value=T)
chr19.index <- grep("Chr19", marker.names, value=F)
chr19.coords <- as.numeric(gsub(pattern="Chr19_", replacement="", chr19.names))
chr20.names <- grep("Chr20", marker.names, value=T)
chr20.index <- grep("Chr20", marker.names, value=F)
chr20.coords <- as.numeric(gsub(pattern="Chr20_", replacement="", chr20.names))


#   Next, we split up the fst values by the chromosome they are on
fsts.gm01 <- fsts$per.loc$FST[chr1.index]
fsts.gm02 <- fsts$per.loc$FST[chr2.index]
fsts.gm03 <- fsts$per.loc$FST[chr3.index]
fsts.gm04 <- fsts$per.loc$FST[chr4.index]
fsts.gm05 <- fsts$per.loc$FST[chr5.index]
fsts.gm06 <- fsts$per.loc$FST[chr6.index]
fsts.gm07 <- fsts$per.loc$FST[chr7.index]
fsts.gm08 <- fsts$per.loc$FST[chr8.index]
fsts.gm09 <- fsts$per.loc$FST[chr9.index]
fsts.gm10 <- fsts$per.loc$FST[chr10.index]
fsts.gm11 <- fsts$per.loc$FST[chr11.index]
fsts.gm12 <- fsts$per.loc$FST[chr12.index]
fsts.gm13 <- fsts$per.loc$FST[chr13.index]
fsts.gm14 <- fsts$per.loc$FST[chr14.index]
fsts.gm15 <- fsts$per.loc$FST[chr15.index]
fsts.gm16 <- fsts$per.loc$FST[chr16.index]
fsts.gm17 <- fsts$per.loc$FST[chr17.index]
fsts.gm18 <- fsts$per.loc$FST[chr18.index]
fsts.gm19 <- fsts$per.loc$FST[chr19.index]
fsts.gm20 <- fsts$per.loc$FST[chr20.index]

#   We get the SPA scores, too
spa.gm01 <- spa[spa$CHR=="Chr01",]
spa.gm02 <- spa[spa$CHR=="Chr02",]
spa.gm03 <- spa[spa$CHR=="Chr03",]
spa.gm04 <- spa[spa$CHR=="Chr04",]
spa.gm05 <- spa[spa$CHR=="Chr05",]
spa.gm06 <- spa[spa$CHR=="Chr06",]
spa.gm07 <- spa[spa$CHR=="Chr07",]
spa.gm08 <- spa[spa$CHR=="Chr08",]
spa.gm09 <- spa[spa$CHR=="Chr09",]
spa.gm10 <- spa[spa$CHR=="Chr10",]
spa.gm11 <- spa[spa$CHR=="Chr11",]
spa.gm12 <- spa[spa$CHR=="Chr12",]
spa.gm13 <- spa[spa$CHR=="Chr13",]
spa.gm14 <- spa[spa$CHR=="Chr14",]
spa.gm15 <- spa[spa$CHR=="Chr15",]
spa.gm16 <- spa[spa$CHR=="Chr16",]
spa.gm17 <- spa[spa$CHR=="Chr17",]
spa.gm18 <- spa[spa$CHR=="Chr18",]
spa.gm19 <- spa[spa$CHR=="Chr19",]
spa.gm20 <- spa[spa$CHR=="Chr20",]
#   And the recombination rates
recomb.gm01 <- recomb[recomb$Chrom=="Gm01",]
recomb.gm02 <- recomb[recomb$Chrom=="Gm02",]
recomb.gm03 <- recomb[recomb$Chrom=="Gm03",]
recomb.gm04 <- recomb[recomb$Chrom=="Gm04",]
recomb.gm05 <- recomb[recomb$Chrom=="Gm05",]
recomb.gm06 <- recomb[recomb$Chrom=="Gm06",]
recomb.gm07 <- recomb[recomb$Chrom=="Gm07",]
recomb.gm08 <- recomb[recomb$Chrom=="Gm08",]
recomb.gm09 <- recomb[recomb$Chrom=="Gm09",]
recomb.gm10 <- recomb[recomb$Chrom=="Gm10",]
recomb.gm11 <- recomb[recomb$Chrom=="Gm11",]
recomb.gm12 <- recomb[recomb$Chrom=="Gm12",]
recomb.gm13 <- recomb[recomb$Chrom=="Gm13",]
recomb.gm14 <- recomb[recomb$Chrom=="Gm14",]
recomb.gm15 <- recomb[recomb$Chrom=="Gm15",]
recomb.gm16 <- recomb[recomb$Chrom=="Gm16",]
recomb.gm17 <- recomb[recomb$Chrom=="Gm17",]
recomb.gm18 <- recomb[recomb$Chrom=="Gm18",]
recomb.gm19 <- recomb[recomb$Chrom=="Gm19",]
recomb.gm20 <- recomb[recomb$Chrom=="Gm20",]

#   calculate a correlation between the two, just out of curiosity

cor_coeff <- cor(fsts$per.loc$FST, spa$Selction_score, method="pearson")
print(cor_coeff)
print(cor_coeff**2)

#   Get the sliding windows!
fsts.gm01.windows <- SlidingWindow(fsts.gm01, chr1.coords, windowsize, step)
fsts.gm02.windows <- SlidingWindow(fsts.gm02, chr2.coords, windowsize, step)
fsts.gm03.windows <- SlidingWindow(fsts.gm03, chr3.coords, windowsize, step)
fsts.gm04.windows <- SlidingWindow(fsts.gm04, chr4.coords, windowsize, step)
fsts.gm05.windows <- SlidingWindow(fsts.gm05, chr5.coords, windowsize, step)
fsts.gm06.windows <- SlidingWindow(fsts.gm06, chr6.coords, windowsize, step)
fsts.gm07.windows <- SlidingWindow(fsts.gm07, chr7.coords, windowsize, step)
fsts.gm08.windows <- SlidingWindow(fsts.gm08, chr8.coords, windowsize, step)
fsts.gm09.windows <- SlidingWindow(fsts.gm09, chr9.coords, windowsize, step)
fsts.gm10.windows <- SlidingWindow(fsts.gm10, chr10.coords, windowsize, step)
fsts.gm11.windows <- SlidingWindow(fsts.gm11, chr11.coords, windowsize, step)
fsts.gm12.windows <- SlidingWindow(fsts.gm12, chr12.coords, windowsize, step)
fsts.gm13.windows <- SlidingWindow(fsts.gm13, chr13.coords, windowsize, step)
fsts.gm14.windows <- SlidingWindow(fsts.gm14, chr14.coords, windowsize, step)
fsts.gm15.windows <- SlidingWindow(fsts.gm15, chr15.coords, windowsize, step)
fsts.gm16.windows <- SlidingWindow(fsts.gm16, chr16.coords, windowsize, step)
fsts.gm17.windows <- SlidingWindow(fsts.gm17, chr17.coords, windowsize, step)
fsts.gm18.windows <- SlidingWindow(fsts.gm18, chr18.coords, windowsize, step)
fsts.gm19.windows <- SlidingWindow(fsts.gm19, chr19.coords, windowsize, step)
fsts.gm20.windows <- SlidingWindow(fsts.gm20, chr20.coords, windowsize, step)
spa.gm01.windows <- SlidingWindow(spa.gm01$Selction_score/maxspa, spa.gm01$BP, windowsize, step)
spa.gm02.windows <- SlidingWindow(spa.gm02$Selction_score/maxspa, spa.gm02$BP, windowsize, step)
spa.gm03.windows <- SlidingWindow(spa.gm03$Selction_score/maxspa, spa.gm03$BP, windowsize, step)
spa.gm04.windows <- SlidingWindow(spa.gm04$Selction_score/maxspa, spa.gm04$BP, windowsize, step)
spa.gm05.windows <- SlidingWindow(spa.gm05$Selction_score/maxspa, spa.gm05$BP, windowsize, step)
spa.gm06.windows <- SlidingWindow(spa.gm06$Selction_score/maxspa, spa.gm06$BP, windowsize, step)
spa.gm07.windows <- SlidingWindow(spa.gm07$Selction_score/maxspa, spa.gm07$BP, windowsize, step)
spa.gm08.windows <- SlidingWindow(spa.gm08$Selction_score/maxspa, spa.gm08$BP, windowsize, step)
spa.gm09.windows <- SlidingWindow(spa.gm09$Selction_score/maxspa, spa.gm09$BP, windowsize, step)
spa.gm10.windows <- SlidingWindow(spa.gm10$Selction_score/maxspa, spa.gm10$BP, windowsize, step)
spa.gm11.windows <- SlidingWindow(spa.gm11$Selction_score/maxspa, spa.gm11$BP, windowsize, step)
spa.gm12.windows <- SlidingWindow(spa.gm12$Selction_score/maxspa, spa.gm12$BP, windowsize, step)
spa.gm13.windows <- SlidingWindow(spa.gm13$Selction_score/maxspa, spa.gm13$BP, windowsize, step)
spa.gm14.windows <- SlidingWindow(spa.gm14$Selction_score/maxspa, spa.gm14$BP, windowsize, step)
spa.gm15.windows <- SlidingWindow(spa.gm15$Selction_score/maxspa, spa.gm15$BP, windowsize, step)
spa.gm16.windows <- SlidingWindow(spa.gm16$Selction_score/maxspa, spa.gm16$BP, windowsize, step)
spa.gm17.windows <- SlidingWindow(spa.gm17$Selction_score/maxspa, spa.gm17$BP, windowsize, step)
spa.gm18.windows <- SlidingWindow(spa.gm18$Selction_score/maxspa, spa.gm18$BP, windowsize, step)
spa.gm19.windows <- SlidingWindow(spa.gm19$Selction_score/maxspa, spa.gm19$BP, windowsize, step)
spa.gm20.windows <- SlidingWindow(spa.gm20$Selction_score/maxspa, spa.gm20$BP, windowsize, step)

#   Build the data frames to plot
fst.gm01.total <- data.frame(
    FST=c(fsts.gm01.windows$FST, spa.gm01.windows$FST, recomb.gm01$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm01.windows$Coords, spa.gm01.windows$Coords, recomb.gm01$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm01.windows$FST)), rep("SPA", length(spa.gm01.windows$FST)), rep("cM/Mb", length(recomb.gm01$cM_per_Mb))),
    p_start=5320000,
    p_end=47200000
)
fst.gm02.total <- data.frame(
    FST=c(fsts.gm02.windows$FST, spa.gm02.windows$FST, recomb.gm02$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm02.windows$Coords, spa.gm02.windows$Coords, recomb.gm02$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm02.windows$FST)), rep("SPA", length(spa.gm02.windows$FST)), rep("cM/Mb", length(recomb.gm02$cM_per_Mb))),
    p_start=16900000,
    p_end=38100000
)
fst.gm03.total <- data.frame(
    FST=c(fsts.gm03.windows$FST, spa.gm03.windows$FST, recomb.gm03$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm03.windows$Coords, spa.gm03.windows$Coords, recomb.gm03$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm03.windows$FST)), rep("SPA", length(spa.gm03.windows$FST)), rep("cM/Mb", length(recomb.gm03$cM_per_Mb))),
    p_start=5480000,
    p_end=31600000
)
fst.gm04.total <- data.frame(
    FST=c(fsts.gm04.windows$FST, spa.gm04.windows$FST, recomb.gm04$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm04.windows$Coords, spa.gm04.windows$Coords, recomb.gm04$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm04.windows$FST)), rep("SPA", length(spa.gm04.windows$FST)), rep("cM/Mb", length(recomb.gm04$cM_per_Mb))),
    p_start=9650000,
    p_end=42800000
)
fst.gm05.total <- data.frame(
    FST=c(fsts.gm05.windows$FST, spa.gm05.windows$FST, recomb.gm05$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm05.windows$Coords, spa.gm05.windows$Coords, recomb.gm05$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm05.windows$FST)), rep("SPA", length(spa.gm05.windows$FST)), rep("cM/Mb", length(recomb.gm05$cM_per_Mb))),
    p_start=7450000,
    p_end=27000000
)
fst.gm06.total <- data.frame(
    FST=c(fsts.gm06.windows$FST, spa.gm06.windows$FST, recomb.gm06$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm06.windows$Coords, spa.gm06.windows$Coords, recomb.gm06$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm06.windows$FST)), rep("SPA", length(spa.gm06.windows$FST)), rep("cM/Mb", length(recomb.gm06$cM_per_Mb))),
    p_start=17900000,
    p_end=47300000
)
fst.gm07.total <- data.frame(
    FST=c(fsts.gm07.windows$FST, spa.gm07.windows$FST, recomb.gm07$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm07.windows$Coords, spa.gm07.windows$Coords, recomb.gm07$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm07.windows$FST)), rep("SPA", length(spa.gm07.windows$FST)), rep("cM/Mb", length(recomb.gm07$cM_per_Mb))),
    p_start=18000000,
    p_end=34700000
)
fst.gm08.total <- data.frame(
    FST=c(fsts.gm08.windows$FST, spa.gm08.windows$FST, recomb.gm08$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm08.windows$Coords, spa.gm08.windows$Coords, recomb.gm08$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm08.windows$FST)), rep("SPA", length(spa.gm08.windows$FST)), rep("cM/Mb", length(recomb.gm08$cM_per_Mb))),
    p_start=23300000,
    p_end=39900000
)
fst.gm09.total <- data.frame(
    FST=c(fsts.gm09.windows$FST, spa.gm09.windows$FST, recomb.gm09$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm09.windows$Coords, spa.gm09.windows$Coords, recomb.gm09$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm09.windows$FST)), rep("SPA", length(spa.gm09.windows$FST)), rep("cM/Mb", length(recomb.gm09$cM_per_Mb))),
    p_start=6330000,
    p_end=37300000
)
fst.gm10.total <- data.frame(
    FST=c(fsts.gm10.windows$FST, spa.gm10.windows$FST, recomb.gm10$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm10.windows$Coords, spa.gm10.windows$Coords, recomb.gm10$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm10.windows$FST)), rep("SPA", length(spa.gm10.windows$FST)), rep("cM/Mb", length(recomb.gm10$cM_per_Mb))),
    p_start=6950000,
    p_end=36800000
)
fst.gm11.total <- data.frame(
    FST=c(fsts.gm11.windows$FST, spa.gm11.windows$FST, recomb.gm11$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm11.windows$Coords, spa.gm11.windows$Coords, recomb.gm11$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm11.windows$FST)), rep("SPA", length(spa.gm11.windows$FST)), rep("cM/Mb", length(recomb.gm11$cM_per_Mb))),
    p_start=13000000,
    p_end=29000000
)
fst.gm12.total <- data.frame(
    FST=c(fsts.gm12.windows$FST, spa.gm12.windows$FST, recomb.gm12$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm12.windows$Coords, spa.gm12.windows$Coords, recomb.gm12$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm12.windows$FST)), rep("SPA", length(spa.gm12.windows$FST)), rep("cM/Mb", length(recomb.gm12$cM_per_Mb))),
    p_start=9300000,
    p_end=32700000
)
fst.gm13.total <- data.frame(
    FST=c(fsts.gm13.windows$FST, spa.gm13.windows$FST, recomb.gm13$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm13.windows$Coords, spa.gm13.windows$Coords, recomb.gm13$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm13.windows$FST)), rep("SPA", length(spa.gm13.windows$FST)), rep("cM/Mb", length(recomb.gm13$cM_per_Mb))),
    p_start=0,
    p_end=9390000
)
fst.gm14.total <- data.frame(
    FST=c(fsts.gm14.windows$FST, spa.gm14.windows$FST, recomb.gm14$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm14.windows$Coords, spa.gm14.windows$Coords, recomb.gm14$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm14.windows$FST)), rep("SPA", length(spa.gm14.windows$FST)), rep("cM/Mb", length(recomb.gm14$cM_per_Mb))),
    p_start=9830000,
    p_end=43600000
)
fst.gm15.total <- data.frame(
    FST=c(fsts.gm15.windows$FST, spa.gm15.windows$FST, recomb.gm15$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm15.windows$Coords, spa.gm15.windows$Coords, recomb.gm15$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm15.windows$FST)), rep("SPA", length(spa.gm15.windows$FST)), rep("cM/Mb", length(recomb.gm15$cM_per_Mb))),
    p_start=15800000,
    p_end=47100000
)
fst.gm16.total <- data.frame(
    FST=c(fsts.gm16.windows$FST, spa.gm16.windows$FST, recomb.gm16$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm16.windows$Coords, spa.gm16.windows$Coords, recomb.gm16$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm16.windows$FST)), rep("SPA", length(spa.gm16.windows$FST)), rep("cM/Mb", length(recomb.gm16$cM_per_Mb))),
    p_start=7330000,
    p_end=27800000
)
fst.gm17.total <- data.frame(
    FST=c(fsts.gm17.windows$FST, spa.gm17.windows$FST, recomb.gm17$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm17.windows$Coords, spa.gm17.windows$Coords, recomb.gm17$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm17.windows$FST)), rep("SPA", length(spa.gm17.windows$FST)), rep("cM/Mb", length(recomb.gm17$cM_per_Mb))),
    p_start=14100000,
    p_end=37300000
)
fst.gm18.total <- data.frame(
    FST=c(fsts.gm18.windows$FST, spa.gm18.windows$FST, recomb.gm18$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm18.windows$Coords, spa.gm18.windows$Coords, recomb.gm18$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm18.windows$FST)), rep("SPA", length(spa.gm18.windows$FST)), rep("cM/Mb", length(recomb.gm18$cM_per_Mb))),
    p_start=17200000,
    p_end=42700000
)
fst.gm19.total <- data.frame(
    FST=c(fsts.gm19.windows$FST, spa.gm19.windows$FST, recomb.gm19$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm19.windows$Coords, spa.gm19.windows$Coords, recomb.gm19$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm19.windows$FST)), rep("SPA", length(spa.gm19.windows$FST)), rep("cM/Mb", length(recomb.gm19$cM_per_Mb))),
    p_start=4560000,
    p_end=34300000
)
fst.gm20.total <- data.frame(
    FST=c(fsts.gm20.windows$FST, spa.gm20.windows$FST, recomb.gm20$cM_per_Mb/maxrecomb),
    pos=c(fsts.gm20.windows$Coords, spa.gm20.windows$Coords, recomb.gm20$PhysicalPosition),
    lab=c(rep("FST", length(fsts.gm20.windows$FST)), rep("SPA", length(spa.gm20.windows$FST)), rep("cM/Mb", length(recomb.gm20$cM_per_Mb))),
    p_start=5390000,
    p_end=32400000
)

#   Start building plots!
p1 <- ggplot(
    fst.gm01.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm01.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm01.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p2 <- ggplot(
    fst.gm02.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm02.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm02.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p3 <- ggplot(
    fst.gm03.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm03.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm03.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p4 <- ggplot(
    fst.gm04.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm04.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm04.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p5 <- ggplot(
    fst.gm05.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm05.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm05.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p6 <- ggplot(
    fst.gm06.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm06.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm06.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p7 <- ggplot(
    fst.gm07.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm07.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm07.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p8 <- ggplot(
    fst.gm08.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm08.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm08.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p9 <- ggplot(
    fst.gm09.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm09.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm09.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p10 <- ggplot(
    fst.gm10.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm10.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm10.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p11 <- ggplot(
    fst.gm11.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm11.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm11.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p12 <- ggplot(
    fst.gm12.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm12.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm12.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p13 <- ggplot(
    fst.gm13.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm13.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm13.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p14 <- ggplot(
    fst.gm14.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm14.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm14.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p15 <- ggplot(
    fst.gm15.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm15.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm15.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p16 <- ggplot(
    fst.gm16.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm16.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm16.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p17 <- ggplot(
    fst.gm17.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm17.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm17.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p18 <- ggplot(
    fst.gm18.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm18.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm18.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p19 <- ggplot(
    fst.gm19.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm19.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm19.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))
p20 <- ggplot(
    fst.gm20.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm20.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm20.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
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
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))

#   Arrange and plot them!

legend <- g_legend(p1)
lwidth <- sum(legend$width)
pdf(file="WholeGenome_FST_SPA.pdf", width=16, height=23)
grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position="none"),
        p2 + theme(legend.position="none"),
        p3 + theme(legend.position="none"),
        p4 + theme(legend.position="none"),
        p5 + theme(legend.position="none"),
        p6 + theme(legend.position="none"),
        p7 + theme(legend.position="none"),
        p8 + theme(legend.position="none"),
        p9 + theme(legend.position="none"),
        p10 + theme(legend.position="none"),
        p11 + theme(legend.position="none"),
        p12 + theme(legend.position="none"),
        p13 + theme(legend.position="none"),
        p14 + theme(legend.position="none"),
        p15 + theme(legend.position="none"),
        p16 + theme(legend.position="none"),
        p17 + theme(legend.position="none"),
        p18 + theme(legend.position="none"),
        p19 + theme(legend.position="none"),
        p20 + theme(legend.position="none"),
        main="",
        left="",
        sub="",
        ncol=2
        ),
    legend,
    widths=unit.c(unit(1, "npc") - lwidth, lwidth),
    nrow=1)
dev.off()
# pdf(file="GM01-05_FST_SPA.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p1 + theme(legend.position="none"),
#         p2 + theme(legend.position="none"),
#         p3 + theme(legend.position="none"),
#         p4 + theme(legend.position="none"),
#         p5 + theme(legend.position="none"),
#         main="FST, SPA Score, and cM/Mb",
#         left="Scaled Value",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
# pdf(file="GM06-10_FST_SPA.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p6 + theme(legend.position="none"),
#         p7 + theme(legend.position="none"),
#         p8 + theme(legend.position="none"),
#         p9 + theme(legend.position="none"),
#         p10 + theme(legend.position="none"),
#         main="FST, SPA Score, and cM/Mb",
#         left="Scaled Value",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
# pdf(file="GM11-15_FST_SPA.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p11 + theme(legend.position="none"),
#         p12 + theme(legend.position="none"),
#         p13 + theme(legend.position="none"),
#         p14 + theme(legend.position="none"),
#         p15 + theme(legend.position="none"),
#         main="FST, SPA Score, and cM/Mb",
#         left="Scaled Value",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()
# pdf(file="GM16-20_FST_SPA.pdf", width=8, height=10.5)
# grid.arrange(
#     arrangeGrob(
#         p16 + theme(legend.position="none"),
#         p17 + theme(legend.position="none"),
#         p18 + theme(legend.position="none"),
#         p19 + theme(legend.position="none"),
#         p20 + theme(legend.position="none"),
#         main="FST, SPA Score, and cM/Mb",
#         left="Scaled Value",
#         sub="Chromosomal Position",
#         ncol=1
#         ),
#     legend,
#     widths=unit.c(unit(1, "npc") - lwidth, lwidth),
#     nrow=1)
# dev.off()

#   Plot a chr15 only view, scaled to the end of chr15
p15_scaled <- ggplot(
    fst.gm15.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    geom_vline(xintercept=fst.gm15.total$p_start, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_vline(xintercept=fst.gm15.total$p_end, colour="grey25", linetype=1, size=0.8, alpha=0.7) +
    geom_line(size=0.5, alpha=0.7) +
    theme(
        axis.text=element_text(size=9, colour="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        legend.title=element_blank(),
        axis.ticks=element_line(colour="black"),
        legend.background=element_blank(),
        legend.key=element_blank(),
        axis.text.x=element_blank()) +
    scale_color_brewer(palette="Set1") +
    labs(title="", x="", y="") +
    coord_cartesian(xlim=c(0,55000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,55000000, 5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.2))

pdf(file="GM15_SPA_FST_cMMb.pdf", 12, 4)
p15_scaled
dev.off()
#   Plot the FST-SPA correlation
spa_fst <- data.frame(
    FST=fsts$per.loc$FST,
    SPA=spa$Selction_score)
sf <- ggplot(
    spa_fst,
    aes(x=FST, y=SPA)) +
    geom_point(size=1.25, alpha=0.6) +
    geom_smooth(method=lm, se=FALSE, size=2) +
    theme(
        axis.text=element_text(size=9, colour="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        legend.title=element_blank(),
        axis.ticks=element_line(colour="black"),
        legend.background=element_blank(),
        legend.key=element_blank()) +
    scale_color_brewer(palette="Set1") +
    theme(axis.text=element_text(size=10)) +
    labs(title="Correlation Between SPA and FST", x="FST", y="SPA Selection Score") +
    coord_cartesian(xlim=c(0,1), ylim=c(0, 9)) +
    scale_x_continuous(breaks=seq(0,1, 0.1)) +
    scale_y_continuous(breaks=seq(0, 9, 1))+
pdf(file="SPA_FST_Corr.pdf", 8, 8)
sf
dev.off()
