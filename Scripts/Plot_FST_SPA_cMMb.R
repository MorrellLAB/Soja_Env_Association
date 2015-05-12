#   Calculate pairwise SNP FST, and plot that in sliding windows across the
#   genome, along with SPA score and cM/Mb.
#   This script will just plot a single chromosome, but it should easily be
#   extended to work for the entire genome.

#   load up the hierfstat library
library(hierfstat)
#   and ggplot2 for making pretty pictures
library(ggplot2)
library(grid)
library(gridExtra)
#   How big to make the windows?
windowsize <<- 5
#   and how many to slide over?
step <<- 3


#   Function yanked from
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
#   This will let us plot a whole page figure with only one legend on the side.
g_legend <- function(p) {
tmp <- ggplotGrob(p)
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
return(legend)
}


#   Create a new function to calculate sliding averages and generate the midpoints
sliding_window <- function(fsts, coords, winsize, winstep) {
    #   Create a vector of window start positions
    points <- seq(from=1, to=length(fsts) - winsize, by=winstep)
    #   Two new vectors to grow that contain the average FST and position
    fst_averages <- c()
    fst_midpoints <- c()
    #   For each window:
    for(i in 1:length(points)) {
        #   If the window is completely within the vector of FST values
        if( (i + winsize) <= length(fsts)) {
            #   We just take the mean of the FST values in the window
            fst_averages[i] <- mean(
                fsts[points[i]: (points[i] + winsize)],
                na.rm=T)
            #   A little risky, but I am going to assume that the FST values
            #   and the coords vector are the same size
            fst_midpoints[i] <- mean(coords[points[i]: (points[i] + winsize)])
        }
        else {
            #   If we are pushed outside the window, we just take until the end
            fst_averages[i] <- mean(
                fsts[points[i]:length(fsts)],
                na.rm=T)
            fst_midpoints[i] <- mean(
                coords[points[i]:length(fsts)])
        }
    }
    #   Put the window means and midpoints into a dataframe to return
    sliding_windows <- data.frame(FST=fst_averages, Coords=fst_midpoints)
    return(sliding_windows)
}

#   Read in all the separate matrices
#   Markers are columns, individuals are rows
#   The first column has an integer that describes population of origin
soja <- read.table("32k_Numeric.txt", header=TRUE, row.names=2)
#   Read in the SPA scores
spa <- read.table("SPA_Scores.txt", header=TRUE)
#   Read the cM/Mb data
recomb <- read.table("6K_cM_per_Mb.txt", header=TRUE)
#   And the pericentromere data
peric <- read.table("Pericentromeres.txt", header=TRUE)
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
marker_names <- names(soja)[2:ncol(soja)]

#   Next, we get the per-chromosome FST values and plot those
#   get the SNPs on a certain chromosome
#       The names
chr1_names <- grep("Chr01", marker_names, value=TRUE)
#       And the indices
chr1_index <- grep("Chr01", marker_names, value=FALSE)
#       Strip out the coordinates
chr1_coords <- as.numeric(gsub(pattern="Chr01_", replacement="", chr1_names))

#   Next, we split up the fst values by the chromosome they are on
fsts_gm01 <- fsts$per.loc$FST[chr1_index]

#   We get the SPA scores, too
spa_gm01 <- spa[spa$CHR == "Chr01",]
#   And the recombination rates. Note the different nomenclature for this file
recomb.gm01 <- recomb[recomb$Chrom == "Gm01",]

#   calculate a correlation between the two, just out of curiosity
cor_coeff <- cor(fsts$per.loc$FST, spa$Selction_score, method="pearson")
print(cor_coeff)
print(cor_coeff ** 2)

#   Get the sliding windows!
fsts_gm01_windows <- sliding_window(
    fsts_gm01,
    chr1_coords,
    windowsize, step)

spa_gm01_windows <- sliding_window(
    spa_gm01$Selction_score / maxspa,
    spa_gm01$BP,
    windowsize,
    step)

#   Build the data frames to plot
fst_gm01.total <- data.frame(
    FST=c(
        fsts_gm01_windows$FST,
        spa_gm01_windows$FST,
        recomb.gm01$cM_per_Mb / maxrecomb),
    pos=c(
        fsts_gm01_windows$Coords,
        spa_gm01_windows$Coords,
        recomb.gm01$PhysicalPosition),
    lab=c(
        rep("FST", length(fsts_gm01_windows$FST)),
        rep("SPA", length(spa_gm01_windows$FST)),
        rep("cM/Mb", length(recomb.gm01$cM_per_Mb))),
    p_start=peric[1, 2],
    p_end=peric[1, 3]
)

#   Start building plots!
p1 <- ggplot(
    fst_gm01.total,
    aes(x=pos, y=FST, color=lab, linetype=lab)) +
    #   Plot the left boundary of the pericentormeric region
    geom_vline(
        xintercept=fst_gm01.total$p_start,
        colour="grey25",
        linetype=1,
        size=0.8,
        alpha=0.7) +
    #   And the right boundary
    geom_vline(
        xintercept=fst_gm01.total$p_end,
        colour="grey25",
        linetype=1,
        size=0.8,
        alpha=0.7) +
    #   Plot the tracks as line plots
    geom_line(size=0.5, alpha=0.7) +
    #   Set the theme. Blank background with black axes
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
    #   Set the palette away from the gross defaults
    scale_color_brewer(palette="Set1") +
    #   Blank axis labels and title
    labs(title="", x="", y="") +
    #   Plot dimensions - biggest chromosome is ~60Mb, we scale all plots to
    #   that value.
    coord_cartesian(xlim=c(0,60000000), ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(0,60000000,5000000), labels=NULL) +
    scale_y_continuous(breaks=seq(0, 1, 0.5))

#   Arrange and plot them!

legend <- g_legend(p1)
lwidth <- sum(legend$width)
pdf(file="WholeGenome_FST_spa_pdf", width=16, height=23)
grid.arrange(
    arrangeGrob(
        p1 + theme(legend.position="none"),
        #   Uncomment these if you have created the other plot objects and
        #   would like to include them
        # p2 + theme(legend.position="none"),
        # p3 + theme(legend.position="none"),
        # p4 + theme(legend.position="none"),
        # p5 + theme(legend.position="none"),
        # p6 + theme(legend.position="none"),
        # p7 + theme(legend.position="none"),
        # p8 + theme(legend.position="none"),
        # p9 + theme(legend.position="none"),
        # p10 + theme(legend.position="none"),
        # p11 + theme(legend.position="none"),
        # p12 + theme(legend.position="none"),
        # p13 + theme(legend.position="none"),
        # p14 + theme(legend.position="none"),
        # p15 + theme(legend.position="none"),
        # p16 + theme(legend.position="none"),
        # p17 + theme(legend.position="none"),
        # p18 + theme(legend.position="none"),
        # p19 + theme(legend.position="none"),
        # p20 + theme(legend.position="none"),
        main="",
        left="",
        sub="",
        ncol=2
        ),
    legend,
    widths=unit.c(unit(1, "npc") - lwidth, lwidth),
    nrow=1)
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
    labs(
        title="Correlation Between SPA and FST",
        x="FST",
        y="SPA Selection Score") +
    coord_cartesian(xlim=c(0,1), ylim=c(0, 9)) +
    scale_x_continuous(breaks=seq(0,1, 0.1)) +
    scale_y_continuous(breaks=seq(0, 9, 1)) +
pdf(file="SPA_FST_Corr.pdf", 8, 8)
sf
dev.off()
