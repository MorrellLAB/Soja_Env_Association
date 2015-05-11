#   load up the hierfstat library
library(hierfstat)

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
combined_12 <- rbind(pop1, pop2)
combined_13 <- rbind(pop1, pop3)
combined_23 <- rbind(pop2, pop3)

#   calculate those FSTs! Since we are working with inbred samples and we have
#   discared heterozygous sites, we use the haploid calculations.
fsts_12 <- wc(combined_12, diploid=F)
fsts_13 <- wc(combined_13, diploid=F)
fsts_23 <- wc(combined_23, diploid=F)
#   And the total FST
fst_total <- wc(dat, diploid=F)

#   What is overall FST?
print(fsts_12$FST)
print(fsts_13$FST)
print(fsts_23$FST)

#   Save all the FST comparisons in a table
fst_table <- data.frame(
    SNPName=marker.names,
    MS_V_I=as.vector(unlist(fsts_12$per.loc)),
    I_V_MN=as.vector(unlist(fsts_23$per.loc)),
    MS_V_MN=as.vector(unlist(fsts_13$per.loc)),
    Total=as.vector(unlist(fst_total$per.loc))
    )

#   Write these out for plotting later
write.table(
    fst_table,
    file="Pairwise_WC_FST_PerSNP.txt",
    quote=FALSE,
    sep="\t",
    row.names=FALSE)
