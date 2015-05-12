#   Apply quality control criteria (missingness, heterozygosity) to SNP
#   genotyping data.
#   Written by Thomas JY Kono

#   Do this so we can pass command line arguments to the script
args <- commandArgs(TRUE)
#   Define the missing data thresholds
#   These are porportions of markers and individuals to filter
missing_cutoff <<- 0.1
het_cutoff <<- 0.1

#   The missing data values
missing_val <<- "U"
het_val <<- "H"

#   This function calculates the proportion of missing data
calc_missing <- function(x) {
    total <- length(x)
    missing <- sum(x == missing_val) / total
    #   we"re going to be using these to keep, so if it is above the cutoff
    #   then we get rid of it
    if(missing > missing_cutoff) {
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}
#   This function calculates the proportion of heterozygosity observed
calc_het <- function(x) {
    total <- length(x)
    het <- sum(x == het_val) / total
    if(het > het_cutoff) {
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}

#   Read the SNP matrix
#       Rows are markers
#       Columns are individuals
snp_table <- read.table(as.character(args[1]), header=T, row.names=1)
dim(snp_table)

#   Which individuals should be dropped?
individuals.missing <- apply(snp_table, 2, calc_missing)
individuals.het <- apply(snp_table, 2, calc_het)

#   And which markers?
markers.missing <- apply(snp_table, 1, calc_missing)
markers.het <- apply(snp_table, 1, calc_het)

#   Now to make a vector of columns to keep
#   If the value fails for either missingness or heterozygosity, the value
#   in the vector will become FALSE, and it will be omitted from the final table
individuals.keep <- apply(cbind(individuals.missing, individuals.het), 1, all)
#   And rows
markers.keep <- apply(cbind(markers.missing, markers.het), 1, all)

#   Who are we keeping?
individuals.keep
markers.keep

#   FALSE rows and columns, i.e., didn"t pass QC,  will be omitted
snp_table <- snp_table[markers.keep, individuals.keep]

#   Double check the dimensions to see that they make sense
dim(snp_table)

#   Write it out!
write.table(
    snp_table,
    file=sub(
        ".txt",
        "_QC.txt",
        as.character(args[1])),
    sep="\t",
    row.names=T,
    quote=F)
