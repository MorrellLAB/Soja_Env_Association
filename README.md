# Soja Environmental Association

This repository contains analysis scripts and small data files used for the *Glycine soja* environmental association manuscript. Some of these data files are available in the manuscript supplement, but will also be made available here for ease of use. Large data files are zipped with [bzip2](http://www.bzip.org).

**Manuscript status: In Preparation**

Contents:
- Data/
    - **6K_Genetic_Physical_Map**: Physical and genetic map positions for the SoySNP6K, as reported in [Lee et al. (2015)](http://link.springer.com/article/10.1007%2Fs11032-015-0209-5)
    - **Pericentromeres**: Chromosome name, pericentromere start, and pericentromere end. Fetched from [SoyBase (Grant et al. 2010)](http://soybase.org).
    - **Soja_533x32416.clst**: Cluster file for the final dataset. ([PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) format)
    - **Soja_533x32416.map**: Physical location for each SNP in the final dataset. ([PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) format)
    - **Soja_533x32416.ped**: Genotyping matrix for the final dataset. ([PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) format)
    - **SoySNP50K_V1-V2_Positions**: Translation table for SoySNP50K ([Song et al. 2013](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0054985)) from Wm82-a1 to Wm82-a2.
    - **Structure_extraparams**: Extra parameters file for [STRUCTURE](http://pritchardlab.stanford.edu/structure.html)
    - **Structure_mainparams**: Main parameters file for [STRUCTURE](http://pritchardlab.stanford.edu/structure.html)
- Scripts/
    - **CreateMap.Rmd**: Plot the samples on a geographic map, colored according to STRUCTURE assignment.
    - **Create_Manhattanplots.rmd**: Create manhattan plots of the GWAS results for each environmental variable found in Table S4.
    - **Folded_SFS.R**: Create a folded site frequency spectrum, separated by cluster.
    - **Genotype_Matrix_to_Fasta.py**: Create population-specific FASTA files for input into [libsequence](http://molpopgen.github.io/libsequence/) tools. Removes monomorphic markers in each population.
    - **LD_Decay.R**: Calculate decay of LD over physical distance in both euchromatic and pericentromeric regions using the minimum RMSE exponential decay curve in [Abecasis et al. 2001](http://www.cell.com/ajhg/abstract/S0002-9297(07)62483-5). 
    - **Make_Nicholson_FST_Inputs.py**: Create input matrices for an *F*<sub>ST</sub> estimator developed by [Nicholson et al. (2002)](http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00357/abstract) from PLINK .ped and .clst files.
    - **Matrix_to_Column.py** Convert from the square matrix output from [LDheatmap](http://cran.r-project.org/web/packages/LDheatmap/index.html) into a two-column format that is easy to parse for calculating an exponential least-squares curve.
    - **Nicholson_FST.R**: Actually run the Nicholson *F*<sub>ST</sub> estimator, as implemented in the R package '[popgen](http://cran.r-project.org/web/packages/popgen/index.html)'
    - **PCA.R**: Perform principle components analysis (PCA) on genetic and environmental data.
    - **Pairwise_LD.R**: Calculate pairwise LD as *D'* ([Lewontin 1964](http://www.genetics.org/content/49/1/49)) between each pair of SNPs in the dataset. Separates euchromatic and pericentromeric markers.
    - **PerSNP_Pairwise_FSTs.R**: Calculate per-SNP *F*<sub>ST</sub> with the [Weir and Cockerham (1984)](http://www.jstor.org/stable/2408641) estimator as implemented in the R package '[hierfstat](http://cran.r-project.org/web/packages/hierfstat/index.html)'
    - **Pericentromeric_FSTs.R**: Compare the distributions of pericentromeric and euchromatic per-SNP *F*<sub>ST</sub>.
    - **Plot_FST_SPA_cMMb.R**: Plot W-C *F*<sub>ST</sub>, SPA score ([Yang et al. 2012](http://www.nature.com/ng/journal/v44/n6/abs/ng.2285.html)), and cM/Mb across the genome 
    - **SNP_QC.R**: Quality control of genotyping matrix. Remove individuals and markers above a missingness threshold.
    - **biophysical_bioclimatic_maps.R**: Generate maps of biophysical and bioclimatic variation across geographic space. "Species Distribution Modeling With R" by Robert J. Hijmans and Jane Elith, 2014.
    - **extract_bioclimatic_biophysical_data.R**: Extract bioclimatic and biophysical data from database files with a set of specified coordinates. Code adapted from "Species Distribution Modeling With R" by Robert J. Hijmans and Jane Elith, 2014.
