#   Simple script to calculate the median physical distance between markers in
#   euchromatic and pericentromeric regions
#   Written by Thomas JY Kono

#   Function to get the gap between markers
spacing <- function(x) {
    nmarkers <- length(x)
    spaces <- c()
    for(i in 2:length(x)) {
        spaces[i-1] <- x[i] - x[i-1]
    }
    return(spaces)
}

chr01_p <- read.table("Chr01_Pericentromere_PhysDistance.txt", header=FALSE)
chr02_p <- read.table("Chr02_Pericentromere_PhysDistance.txt", header=FALSE)
chr03_p <- read.table("Chr03_Pericentromere_PhysDistance.txt", header=FALSE)
chr04_p <- read.table("Chr04_Pericentromere_PhysDistance.txt", header=FALSE)
chr05_p <- read.table("Chr05_Pericentromere_PhysDistance.txt", header=FALSE)
chr06_p <- read.table("Chr06_Pericentromere_PhysDistance.txt", header=FALSE)
chr07_p <- read.table("Chr07_Pericentromere_PhysDistance.txt", header=FALSE)
chr08_p <- read.table("Chr08_Pericentromere_PhysDistance.txt", header=FALSE)
chr09_p <- read.table("Chr09_Pericentromere_PhysDistance.txt", header=FALSE)
chr10_p <- read.table("Chr10_Pericentromere_PhysDistance.txt", header=FALSE)
chr11_p <- read.table("Chr11_Pericentromere_PhysDistance.txt", header=FALSE)
chr12_p <- read.table("Chr12_Pericentromere_PhysDistance.txt", header=FALSE)
chr13_p <- read.table("Chr13_Pericentromere_PhysDistance.txt", header=FALSE)
chr14_p <- read.table("Chr14_Pericentromere_PhysDistance.txt", header=FALSE)
chr15_p <- read.table("Chr15_Pericentromere_PhysDistance.txt", header=FALSE)
chr16_p <- read.table("Chr16_Pericentromere_PhysDistance.txt", header=FALSE)
chr17_p <- read.table("Chr17_Pericentromere_PhysDistance.txt", header=FALSE)
chr18_p <- read.table("Chr18_Pericentromere_PhysDistance.txt", header=FALSE)
chr19_p <- read.table("Chr19_Pericentromere_PhysDistance.txt", header=FALSE)
chr20_p <- read.table("Chr20_Pericentromere_PhysDistance.txt", header=FALSE)

chr01_p <- as.numeric(chr01_p$V1)
chr02_p <- as.numeric(chr02_p$V1)
chr03_p <- as.numeric(chr03_p$V1)
chr04_p <- as.numeric(chr04_p$V1)
chr05_p <- as.numeric(chr05_p$V1)
chr06_p <- as.numeric(chr06_p$V1)
chr07_p <- as.numeric(chr07_p$V1)
chr08_p <- as.numeric(chr08_p$V1)
chr09_p <- as.numeric(chr09_p$V1)
chr10_p <- as.numeric(chr10_p$V1)
chr11_p <- as.numeric(chr11_p$V1)
chr12_p <- as.numeric(chr12_p$V1)
chr13_p <- as.numeric(chr13_p$V1)
chr14_p <- as.numeric(chr14_p$V1)
chr15_p <- as.numeric(chr15_p$V1)
chr16_p <- as.numeric(chr16_p$V1)
chr17_p <- as.numeric(chr17_p$V1)
chr18_p <- as.numeric(chr18_p$V1)
chr19_p <- as.numeric(chr19_p$V1)
chr20_p <- as.numeric(chr20_p$V1)

chr01_e <- read.table("Chr01_Euchromatin_PhysDistance.txt", header=FALSE)
chr02_e <- read.table("Chr02_Euchromatin_PhysDistance.txt", header=FALSE)
chr03_e <- read.table("Chr03_Euchromatin_PhysDistance.txt", header=FALSE)
chr04_e <- read.table("Chr04_Euchromatin_PhysDistance.txt", header=FALSE)
chr05_e <- read.table("Chr05_Euchromatin_PhysDistance.txt", header=FALSE)
chr06_e <- read.table("Chr06_Euchromatin_PhysDistance.txt", header=FALSE)
chr07_e <- read.table("Chr07_Euchromatin_PhysDistance.txt", header=FALSE)
chr08_e <- read.table("Chr08_Euchromatin_PhysDistance.txt", header=FALSE)
chr09_e <- read.table("Chr09_Euchromatin_PhysDistance.txt", header=FALSE)
chr10_e <- read.table("Chr10_Euchromatin_PhysDistance.txt", header=FALSE)
chr11_e <- read.table("Chr11_Euchromatin_PhysDistance.txt", header=FALSE)
chr12_e <- read.table("Chr12_Euchromatin_PhysDistance.txt", header=FALSE)
chr13_e <- read.table("Chr13_Euchromatin_PhysDistance.txt", header=FALSE)
chr14_e <- read.table("Chr14_Euchromatin_PhysDistance.txt", header=FALSE)
chr15_e <- read.table("Chr15_Euchromatin_PhysDistance.txt", header=FALSE)
chr16_e <- read.table("Chr16_Euchromatin_PhysDistance.txt", header=FALSE)
chr17_e <- read.table("Chr17_Euchromatin_PhysDistance.txt", header=FALSE)
chr18_e <- read.table("Chr18_Euchromatin_PhysDistance.txt", header=FALSE)
chr19_e <- read.table("Chr19_Euchromatin_PhysDistance.txt", header=FALSE)
chr20_e <- read.table("Chr20_Euchromatin_PhysDistance.txt", header=FALSE)

chr01_e <- as.numeric(chr01_e$V1)
chr02_e <- as.numeric(chr02_e$V1)
chr03_e <- as.numeric(chr03_e$V1)
chr04_e <- as.numeric(chr04_e$V1)
chr05_e <- as.numeric(chr05_e$V1)
chr06_e <- as.numeric(chr06_e$V1)
chr07_e <- as.numeric(chr07_e$V1)
chr08_e <- as.numeric(chr08_e$V1)
chr09_e <- as.numeric(chr09_e$V1)
chr10_e <- as.numeric(chr10_e$V1)
chr11_e <- as.numeric(chr11_e$V1)
chr12_e <- as.numeric(chr12_e$V1)
chr13_e <- as.numeric(chr13_e$V1)
chr14_e <- as.numeric(chr14_e$V1)
chr15_e <- as.numeric(chr15_e$V1)
chr16_e <- as.numeric(chr16_e$V1)
chr17_e <- as.numeric(chr17_e$V1)
chr18_e <- as.numeric(chr18_e$V1)
chr19_e <- as.numeric(chr19_e$V1)
chr20_e <- as.numeric(chr20_e$V1)

peri <- c(
    spacing(chr01_p),
    spacing(chr02_p),
    spacing(chr03_p),
    spacing(chr04_p),
    spacing(chr05_p),
    spacing(chr06_p),
    spacing(chr07_p),
    spacing(chr08_p),
    spacing(chr09_p),
    spacing(chr10_p),
    spacing(chr11_p),
    spacing(chr12_p),
    spacing(chr13_p),
    spacing(chr14_p),
    spacing(chr15_p),
    spacing(chr16_p),
    spacing(chr17_p),
    spacing(chr18_p),
    spacing(chr19_p),
    spacing(chr20_p)
    )
euch <- c(
    spacing(chr01_e),
    spacing(chr02_e),
    spacing(chr03_e),
    spacing(chr04_e),
    spacing(chr05_e),
    spacing(chr06_e),
    spacing(chr07_e),
    spacing(chr08_e),
    spacing(chr09_e),
    spacing(chr10_e),
    spacing(chr11_e),
    spacing(chr12_e),
    spacing(chr13_e),
    spacing(chr14_e),
    spacing(chr15_e),
    spacing(chr16_e),
    spacing(chr17_e),
    spacing(chr18_e),
    spacing(chr19_e),
    spacing(chr20_e)
    )

print(median(peri))
print(median(euch))
