#   A script to calculate and plot LD as D' and r-squared for pericentromeric
#   and euchromatic regions
#       Written by Thomas JY Kono

#   Calculate LD statistics
library(LDheatmap)
#   To color things nicely
library(RColorBrewer)
#   To work with genetic data
library(genetics)
#   Take arguments
args <- commandArgs(TRUE)

#   Define functions here
#       A function to get the numeric part of a marker name, which correpsonds
#       to the physical position
marker_phys_pos <- function(x) {
    #   Marker names are of the form
    #       CHROM_POSITION
    name_pieces <- unlist(strsplit(x, "_"))
    pos <- as.numeric(name_pieces[2])
    return(pos)
}
#       A function to replace our 1/0 genotype data with A/A T/T allele calls
#       This is because LDheatmap expects "genotype" objects, which in turn
#       expect the data to be of this form.
convert_geno <- function(x) {
    geno <- x
    geno[geno == 1] <- "A/A"
    geno[geno == 0] <- "T/T"
    return(geno)
}
#       A function to ask whether or not a marker is pericentromeric
centromeric <- function(x) {
    #   Split the chromosome name and physical position
    name_pieces <- unlist(strsplit(x, "_"))
    chrom <- name_pieces[1]
    pos <- as.numeric(name_pieces[2])
    #   Which chromosome are we looking at?
    regions <- centromeres[centromeres$Chromosome == chrom,]
    peri_start <- as.numeric(regions[2])
    peri_end <- as.numeric(regions[3])
    if( peri_start < pos & pos < peri_end ) {
        return(TRUE)
    }
    else {
        return(FALSE)
    }
}

#   Which chromosomes do we want to calcuate on?
chrom <- args[1]

#   Read in the centromeric regions
centromeres <- read.table("Pericentromeres.txt", header=TRUE)

#   Read in the genotype data
genotypes <- read.table("32k_Numeric_NewStruct.txt", header=TRUE)
#   We want to get the marker names that correspond to the chromosome that
#   we are interested in. The first two columns are population of origin
#   and the PI number.
marker_names <- colnames(genotypes)[3:ncol(genotypes)]
#   Which markers are on that chromosome?
#   With value=TRUE, we get the actual marker names
chrom_markers <- grep(chrom, marker_names, value=TRUE)
#   Then, we get the physical positions
physical_positions <- sapply(chrom_markers, marker_phys_pos)
#   Then get the marker data from the genotyping matrix
marker_data <- genotypes[,chrom_markers]
#   And convert it into the proper format
genotypes <- apply(marker_data, 2, convert_geno)
#   Which markers are pericentromeric?
peri_markers <- sapply(chrom_markers, centromeric)
#   Split the genotype matrices by which ones are pericentromeric and
#   euchromatic and convert them into Genotype objects
peri_genotypes <- makeGenotypes(genotypes[,peri_markers], tol=0)
peri_distances <- physical_positions[peri_markers]
euch_genotypes <- makeGenotypes(genotypes[,!peri_markers], tol=0)
euch_distances <- physical_positions[!peri_markers]
#   Open a new PDF file for the heatmap
pdf(file=paste(chrom, "_Pericentromere.pdf", sep=""), 8, 8)
peri_ldh <- LDheatmap(
    peri_genotypes,
    genetic.distances=peri_distances,
    LDmeasure="D'",
    add.map=TRUE,
    color=gray.colors(25),
    flip=TRUE,
    title=paste("LD as D' in", chrom, "Pericentromere")
    )
dev.off()
#   Open a new PDF file for the heatmap
pdf(file=paste(chrom, "_Euchromatin.pdf", sep=""), 8, 8)
euch_ldh <- LDheatmap(
    euch_genotypes,
    genetic.distances=euch_distances,
    LDmeasure="D'",
    add.map=TRUE,
    color=gray.colors(25),
    flip=TRUE,
    title=paste("LD as D' in", chrom, "Euchromatin")
    )
dev.off()

#   Save the LD matrices and the makes of markers for each
write.table(
    peri_ldh$LDmatrix,
    file=paste(chrom, "_Pericentromere_LDmatrix.txt", sep=""),
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE,
    sep="\t")
#   And the distances used. This will aid in regression later.
write.table(
    peri_ldh$genetic.distances,
    file=paste(chrom, "_Pericentromere_PhysDistance.txt", sep=""),
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE,
    sep="\n")
#   Do the same for euchromatic
write.table(
    euch_ldh$LDmatrix,
    file=paste(chrom, "_Euchromatin_LDmatrix.txt", sep=""),
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE,
    sep="\t")
#   And the distances used. This will aid in regression later.
write.table(
    euch_ldh$genetic.distances,
    file=paste(chrom, "_Euchromatin_PhysDistance.txt", sep=""),
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE,
    sep="\n")
