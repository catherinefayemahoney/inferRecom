#' Example PLINK Dataset - simCEU
#'
#' A simulated CEU (Utah residents with Northern and Western European ancestry)
#' dataset for demonstrating crossover detection and ROH analysis.
#'
#' @format PLINK binary format files (.bed, .bim, .fam) containing:
#' \describe{
#'   \item{Samples}{Simulated family trios and larger pedigrees with 2-3 children}
#'   \item{SNPs}{Chromosome 4 markers}
#'   \item{Genotypes}{Simulated genotype data with realistic LD structure}
#'   \item{Families}{Mix of 2-child and 3-child families for testing}
#' }
#'
#' @details
#' The dataset includes three PLINK binary files:
#' \itemize{
#'   \item \code{simCEU.bed} - Binary genotype data
#'   \item \code{simCEU.bim} - Variant information (SNP IDs, positions, alleles)
#'   \item \code{simCEU.fam} - Sample information (family IDs, relationships)
#' }
#'
#' Access the files using:
#' \preformatted{
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#' }
#'
#' @source Simulated data based on CEU population structure from 1000 Genomes
#' Project and R package \code{sim1000G}
#'
#' @return No R object is returned. Access files via
#' \code{system.file("extdata", package = "inferRecom")}.
#'
#' @references
#' Siva, N. (2008). 1000 Genomes project.
#'
#' Dimitromanolakis, A., Xu, J., Krol, A., & Briollais, L. (2019). sim1000G:
#' a user-friendly genetic variant simulator in R for unrelated individuals and
#' family-based designs. BMC bioinformatics, 20(1), 26.
#' @name simCEU-data
NULL
#'
#' Genetic Maps for Chromosome 4
#'
#' Sex-specific and sex-averaged genetic maps for chromosome 4.
#'
#' @format Tab-delimited text file with columns:
#' \describe{
#'   \item{chr}{Chromosome identifier (chr4)}
#'   \item{pos}{Physical position in base pairs (hg19/GRCh37)}
#'   \item{rate}{Recombination rate (cM/Mb) at this position}
#'   \item{cM}{Cumulative genetic distance in centiMorgans from chromosome start}
#' }
#'
#' Access the files using:
#' \preformatted{
#' femaleMap <- read.delim(file.path(dataPath, "female_chr4.txt"))
#' }
#'
#' @source Derived from European sex-specific maps from Bherer, et al. (2014).
#'
#' @return No R object is returned. Access files via
#' \code{system.file("extdata", package = "inferRecom")}.
#'
#' @references
#' Bhérer, C., Campbell, C. L., & Auton, A. (2017). Refined genetic maps reveal
#' sexual dimorphism in human meiotic recombination at multiple scales. Nature
#' communications, 8(1), 14994.
#' @name geneticMaps-data
NULL
#' Example Crossover Data
#'
#' Pre-computed maternal and paternal crossover events from the simCEU dataset.
#' Used for testing and demonstrating haplotype phasing functionality.
#'
#' @format \link[GenomicRanges:GRanges-class]{GRanges} objects stored as RDS
#' files with metadata columns:
#' \describe{
#'   \item{childId}{Child identifier}
#'   \item{familyId}{Family identifier}
#'   \item{startSnp}{SNP name at crossover interval start}
#'   \item{finishSnp}{SNP name at crossover interval end}
#'   \item{startCm}{Genetic position at interval start (cM)}
#'   \item{finishCm}{Genetic position at interval end (cM)}
#' }
#'
#' @details
#' Two crossover datasets are provided:
#' \itemize{
#'   \item \code{xoMat.rds} - Maternal crossover events
#'   \item \code{xoPat.rds} - Paternal crossover events
#' }
#'
#' These objects contain detected crossover events from 3-child families
#' in the simCEU dataset. Crossovers were identified using:
#'
#' \preformatted{
#' # Maternal crossovers
#' xoMat <- xoDetect(
#'   plinkFile = "simCEU",
#'   mapFile = "female_chr4.txt",
#'   familySize = 3,
#'   parent = "mother",
#'   snpFilter = 5,
#'   cmFilter = 1
#' )
#'
#' # Paternal crossovers
#' xoPat <- xoDetect(
#'   plinkFile = "simCEU",
#'   mapFile = "male_chr4.txt",
#'   familySize = 3,
#'   parent = "father",
#'   snpFilter = 5,
#'   cmFilter = 1
#' )
#' }
#'
#' Access the crossover data using:
#' \preformatted{
#' dataPath <- system.file("extdata", package = "inferRecom")
#' xoMat <- readRDS(file.path(dataPath, "xoMat.rds"))
#' xoPat <- readRDS(file.path(dataPath, "xoPat.rds"))
#' }
#'
#' @name crossoverData
#'
#' @return No R object is returned. Access files via
#' \code{system.file("extdata", package = "inferRecom")}.
NULL
