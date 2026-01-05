#' Phase Haplotypes Using Crossover Information
#'
#' This function phases parental and child haplotypes based on detected
#' crossover events from maternal and paternal meioses.
#'
#' @param plinkFile Character string; path to the PLINK prefix (without
#'     extensions).
#' @param xoDetectPaternal \link[GenomicRanges:GRanges-class]{GRanges}; output
#'     from \code{xoDetect()} for paternal crossovers (parent = "father").
#' @param xoDetectMaternal \link[GenomicRanges:GRanges-class]{GRanges}; output from
#'     \code{xoDetect()} for maternal crossovers (parent = "mother").
#' @param famIds Character vector or NULL; specific family IDs to phase. If
#'     NULL, all families with 5+ members are analyzed. Default is NULL.
#' @param rsOnly Logical; whether to restrict to SNPs with names beginning
#'     with "rs". Default is TRUE.
#' @param outputFormat Character string; format for output. Options are "list"
#'     (default, returns named list of DataFrames), "summarizedExperiment"
#'     (returns \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} object),
#'     or "vcf" (writes VCF files).
#' @param vcfOutput Character string; path prefix for VCF output files. Only
#'     used when outputFormat = "vcf". Each family will be written to a
#'     separate VCF file with suffix "_<familyID>.vcf". Default is "phased".
#' @param BPPARAM A BiocParallelParam object specifying parallel execution.
#'     Default is \code{SerialParam()} for serial execution. Use
#'     \code{MulticoreParam(workers = n)} for parallel processing across n cores.
#'
#' @return Depends on outputFormat:
#'     \itemize{
#'       \item "list": A named list of \link[S4Vectors:DataFrame-class]{DataFrame}
#'             objects, one per family. Each DataFrame contains phased haplotypes
#'             with columns for SNP ID, physical position, and phased alleles for each
#'             family member.
#'       \item "summarizedExperiment": A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}
#'             object containing phased genotypes across all families with
#'             rowData (SNP information) and colData (sample information).
#'       \item "vcf": Writes VCF files and returns paths to created files.
#'     }
#'     Returns a message if no phase information is available for a family.
#'
#' @details
#' This function performs haplotype phasing using family-based genetic data
#' and detected crossover events. The phasing process:
#' \enumerate{
#'   \item Identifies informative SNPs where parents are heterozygous/homozygous
#'   \item Assigns alleles (1 or 2) to children based on inheritance
#'   \item Converts allele codes to nucleotides (A, C, G, T)
#'   \item Imputes homozygous sites from the other parent
#'   \item Phases parental haplotypes using crossover breakpoints
#'   \item Adds back non-informative SNPs with inferred genotypes where possible
#'   \item Uses neighboring informative SNPs to infer phase for ambiguous cases
#'   \item Outputs diploid genotypes for all family members across all SNPs
#' }
#'
#' Non-informative SNPs are handled as follows:
#' \itemize{
#'   \item Both parents homozygous (same allele): Genotypes confidently inferred
#'   \item Both parents heterozygous: Phase uncertain, set to NA
#'   \item One parent heterozygous, one homozygous: Phase inferred from
#'         neighboring informative SNPs using haplotype continuity
#' }
#'
#' The function requires crossover detection results from both parents as
#' GRanges objects. Run \code{xoDetect()} separately for maternal and paternal
#' crossovers before using this function.
#'
#' Parallel processing is performed by family, improving efficiency for datasets
#' with many families.
#'
#' Output format: Each phased family contains columns:
#' \itemize{
#'   \item \code{rsID} - SNP identifier
#'   \item \code{location} - Physical position in base pairs
#'   \item For each child: \code{<childID>Pat}, \code{<childID>Mat} - paternal
#'         and maternal haplotypes
#'   \item For parents: \code{<parentID>_1}, \code{<parentID>_2} - two
#'         haplotypes
#' }
#'
#' @import methods
#' @importFrom snpStats read.plink
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges mcols seqnames
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom utils read.delim write.table
#' @importFrom BiocParallel SerialParam MulticoreParam bplapply
#' @importFrom stats setNames
#'
#' @examples
#' # Load example data
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#' mapFemale <- file.path(dataPath, "female_chr4.txt")
#' mapMale <- file.path(dataPath, "male_chr4.txt")
#'
#' # Detect crossovers for both parents (returns GRanges objects)
#' xoPat <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapMale,
#'     familySize = 3,
#'     parent = "father"
#' )
#' xoMat <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 3,
#'     parent = "mother"
#' )
#'
#' # Phase haplotypes (serial execution, default)
#' phased <- xoPhase(
#'     plinkFile = plinkFile,
#'     xoDetectPaternal = xoPat,
#'     xoDetectMaternal = xoMat
#' )
#'
#' # Access phased haplotypes for a specific family
#' family1Phase <- phased[[1]]
#' family1Phase
#'
#' # Phase haplotypes for subset of families
#' phasedSubset <- xoPhase(
#'     plinkFile = plinkFile,
#'     xoDetectPaternal = xoPat,
#'     xoDetectMaternal = xoMat,
#'     famIds = c("F4", "F6")
#' )
#'
#' phasedSubset
#'
#' @export
xoPhase <- function(plinkFile,
                    xoDetectPaternal,
                    xoDetectMaternal,
                    famIds = NULL,
                    rsOnly = TRUE,
                    outputFormat = c("list", "summarizedExperiment", "vcf"),
                    vcfOutput = "phased",
                    BPPARAM = SerialParam()) {
    ## Argument validation
    if (!is.character(plinkFile) || length(plinkFile) != 1L) {
        stop("'plinkFile' must be a single character string")
    }

    if (!is(xoDetectPaternal, "GRanges")) {
        stop("'xoDetectPaternal' must be a GRanges object")
    }

    if (!is(xoDetectMaternal, "GRanges")) {
        stop("'xoDetectMaternal' must be a GRanges object")
    }

    if (!is.null(famIds) && !is.character(famIds)) {
        stop("'famIds' must be a character vector or NULL")
    }

    if (!is.logical(rsOnly) || length(rsOnly) != 1L) {
        stop("'rsOnly' must be a single logical value")
    }

    outputFormat <- match.arg(outputFormat)

    if (!is.character(vcfOutput) || length(vcfOutput) != 1L) {
        stop("'vcfOutput' must be a single character string")
    }

    ## Check required packages
    pkgs <- c("snpStats", "GenomicRanges", "BiocParallel")
    miss <- pkgs[!vapply(pkgs, requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )]
    if (length(miss) > 0L) {
        stop("Missing required packages: ", paste(miss, collapse = ", "))
    }

    if (outputFormat == "summarizedExperiment" &&
        !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("Package 'SummarizedExperiment' is required for this output format")
    }

    if (!.checkSameChromosome(xoDetectMaternal) | !.checkSameChromosome(xoDetectPaternal)) {
        stop("GRanges object can only contain ranges from a single chromosome")
    }

    ## Convert GRanges to data frames for internal processing
    xoDetectPaternal <- .granges2DataFrame(xoDetectPaternal)
    xoDetectMaternal <- .granges2DataFrame(xoDetectMaternal)

    ## Strip extension from plinkFile if present
    plinkFile <- sub("\\.(bed|bim|fam)$", "", plinkFile, ignore.case = TRUE)

    ## Read PLINK data
    ped <- snpStats::read.plink(plinkFile)

    ## Convert SnpMatrix to numeric matrix
    genoSnpmat <- ped$genotypes
    if (inherits(genoSnpmat, "SnpMatrix")) {
        geno <- t(methods::as(genoSnpmat, "numeric"))
    } else {
        geno <- t(as.matrix(genoSnpmat))
    }

    fam <- ped$fam

    ## Read BIM file for allele information
    bimFile <- paste0(plinkFile, ".bim")
    if (!file.exists(bimFile)) {
        stop("BIM file not found: ", bimFile)
    }
    bim <- utils::read.delim(bimFile,
        header = FALSE, sep = "\t",
        stringsAsFactors = FALSE, row.names = NULL
    )
    colnames(bim) <- c("chr", "snp", "cm", "pos", "a1", "a2")

    ## Filter for rs SNPs if requested
    if (isTRUE(rsOnly)) {
        rsIdx <- grep("^rs", rownames(geno), ignore.case = TRUE)
        if (length(rsIdx) == 0L) {
            stop("No SNPs with 'rs' prefix found")
        }
        geno <- geno[rsIdx, , drop = FALSE]
    }

    ## Determine families to analyze
    if (!is.null(famIds)) {
        familyIds <- famIds
    } else {
        famCounts <- table(fam[, 1])
        famKeep <- names(famCounts)[famCounts >= 5]
        fam <- fam[fam[, 1] %in% famKeep, , drop = FALSE]
        familyIds <- unique(fam[, 1])
    }

    ## Crossover data is already converted from GRanges to data frames
    xoPat <- xoDetectPaternal
    xoMat <- xoDetectMaternal

    ## Phase each family using parallel processing
    phaseList <- BiocParallel::bplapply(familyIds, function(id) {
        .phaseSingleFamily(
            familyId = id,
            geno = geno,
            fam = fam,
            bim = bim,
            xoPat = xoPat,
            xoMat = xoMat
        )
    }, BPPARAM = BPPARAM)

    ## Set names
    names(phaseList) <- familyIds

    ## Return based on output format
    if (outputFormat == "list") {
        return(phaseList)
    } else if (outputFormat == "summarizedExperiment") {
        return(.convertToSummarizedExperiment(phaseList, bim))
    } else if (outputFormat == "vcf") {
        return(.writeVcfFiles(phaseList, bim, vcfOutput))
    }
}


## ========================== INTERNAL HELPERS ===========================

#' Check if all chromosomes in GRanges are the same
#'
#' @param gr GRanges object
#'
#' @return Logical; TRUE if all chromosomes are identical
#' @keywords internal
.checkSameChromosome <- function(gr) {
    if (length(gr) == 0) {
        return(TRUE)
    }

    chroms <- as.character(GenomicRanges::seqnames(gr))
    length(unique(chroms)) == 1
}

#' Phase a Single Family
#'
#' @param familyId Family identifier
#' @param geno Numeric genotype matrix
#' @param fam FAM data frame
#' @param bim BIM data frame
#' @param xoPat Paternal crossover data frame
#' @param xoMat Maternal crossover data frame
#'
#' @return DataFrame with phased haplotypes or character message
#' @keywords internal
.phaseSingleFamily <- function(familyId, geno, fam, bim, xoPat, xoMat) {
    ## Extract family subset
    famSubset <- fam[fam[, 1] == familyId, ]

    ## Identify parents
    motherId <- .firstNonNA(famSubset$mother[!is.na(famSubset$mother)])
    fatherId <- .firstNonNA(famSubset$father[!is.na(famSubset$father)])

    ## Extract genotypes
    genoSubset <- geno[, which(colnames(geno) %in% famSubset[, 2]),
        drop = FALSE
    ]

    if (is.vector(genoSubset) || ncol(genoSubset) < 3) {
        return("No phase information")
    }

    ## Get parent column indices
    motherCol <- which(colnames(genoSubset) == motherId)
    fatherCol <- which(colnames(genoSubset) == fatherId)

    if (length(motherCol) == 0L || length(fatherCol) == 0L) {
        return("No phase information")
    }

    ## Reorder: father, mother, children
    childCols <- setdiff(
        seq_len(ncol(genoSubset)),
        c(fatherCol, motherCol)
    )
    genoSubset <- genoSubset[, c(fatherCol, motherCol, childCols),
        drop = FALSE
    ]
    famSubset <- famSubset[c(fatherCol, motherCol, childCols), ]

    ## Get informative genotypes
    genoAll <- .getInformativeGenotypes(genoSubset)

    if (is.null(genoAll)) {
        return("No phase information")
    }

    ## Assign allele codes (1 or 2) for each parent
    patAllele <- .assignAlleleCodes(genoAll,
        parent = "father",
        fatherCol = 1L, motherCol = 2L
    )
    matAllele <- .assignAlleleCodes(genoAll,
        parent = "mother",
        fatherCol = 1L, motherCol = 2L
    )

    ## Convert to nucleotides
    patActg <- .convertToNucleotides(patAllele, bim)
    matActg <- .convertToNucleotides(matAllele, bim)

    ## Add homozygous SNPs from other parent
    actgCombPat <- .addHomozygousSNPs(patActg, matActg, genoAll,
        bim,
        parent = "father"
    )
    actgCombMat <- .addHomozygousSNPs(matActg, patActg, genoAll,
        bim,
        parent = "mother"
    )

    ## Phase children
    phaseChildren <- .phaseChildren(actgCombPat, actgCombMat)

    ## Phase parents using crossover information
    patPhase <- .phaseParent(actgCombPat, xoPat, bim, fatherId)
    matPhase <- .phaseParent(actgCombMat, xoMat, bim, motherId)

    ## Align parent phases with children (they should have same SNPs)
    ## Get SNP order from phaseChildren
    childSnps <- as.character(phaseChildren$rsID)

    ## Align parent phases to match children SNP order
    patPhaseAligned <- patPhase[match(childSnps, rownames(patPhase)), , drop = FALSE]
    matPhaseAligned <- matPhase[match(childSnps, rownames(matPhase)), , drop = FALSE]

    ## Reset row names to match children
    rownames(patPhaseAligned) <- NULL
    rownames(matPhaseAligned) <- NULL

    ## Combine phased informative SNPs
    famPhaseInformative <- cbind(phaseChildren, patPhaseAligned, matPhaseAligned)

    ## Add non-informative SNPs back to the output
    famPhaseFull <- .addNonInformativeSNPs(
        phasedData = famPhaseInformative,
        genoSubset = genoSubset,
        bim = bim,
        fatherId = fatherId,
        motherId = motherId,
        fatherCol = 1L,
        motherCol = 2L
    )

    DataFrame(famPhaseFull)
}

#' Add Non-Informative SNPs to Phased Output
#'
#' Adds SNPs that were not informative for phasing back into the output with
#' genotypes inferred from parental homozygous states where possible. For cases
#' where one parent is homozygous and one is heterozygous, uses haplotype
#' information from neighboring informative SNPs to infer phase.
#'
#' @param phasedData Data frame with phased informative SNPs
#' @param genoSubset Full genotype matrix for the family
#' @param bim BIM data frame with SNP information
#' @param fatherId Father identifier
#' @param motherId Mother identifier
#' @param fatherCol Column index for father in genoSubset
#' @param motherCol Column index for mother in genoSubset
#'
#' @return Data frame with all SNPs (informative and non-informative)
#' @keywords internal
.addNonInformativeSNPs <- function(phasedData, genoSubset, bim,
                                   fatherId, motherId,
                                   fatherCol, motherCol) {
    ## Get all SNPs from genotype data
    allSnps <- rownames(genoSubset)
    phasedSnps <- as.character(phasedData$rsID)

    ## Identify non-informative SNPs
    nonInfSnps <- setdiff(allSnps, phasedSnps)

    if (length(nonInfSnps) == 0) {
        ## All SNPs were informative, return as is
        return(phasedData)
    }

    ## Get child column names from phased data
    childCols <- setdiff(
        colnames(phasedData),
        c(
            "rsID", "location",
            paste0(fatherId, "_1"), paste0(fatherId, "_2"),
            paste0(motherId, "_1"), paste0(motherId, "_2")
        )
    )

    ## Extract child IDs (remove Pat/Mat suffix)
    childIds <- unique(gsub("(Pat|Mat)$", "", childCols))

    ## Create template for non-informative SNPs
    nNonInf <- length(nonInfSnps)
    nCols <- ncol(phasedData)

    nonInfMatrix <- matrix(NA_character_, nrow = nNonInf, ncol = nCols)
    colnames(nonInfMatrix) <- colnames(phasedData)

    nonInfDf <- as.data.frame(nonInfMatrix, stringsAsFactors = FALSE)

    ## Create position-indexed map of phased data for neighbor lookup
    phasedPosMap <- setNames(
        seq_len(nrow(phasedData)),
        as.character(phasedData$location)
    )

    ## Fill in SNP IDs and locations
    for (i in seq_len(nNonInf)) {
        snpId <- nonInfSnps[i]
        bimIdx <- which(bim$snp == snpId)

        if (length(bimIdx) > 0) {
            nonInfDf$rsID[i] <- snpId
            snpPos <- bim$pos[bimIdx[1]]
            nonInfDf$location[i] <- as.character(snpPos)

            ## Get genotypes for this SNP
            fatherGeno <- genoSubset[snpId, fatherCol]
            motherGeno <- genoSubset[snpId, motherCol]

            ## Get alleles
            allele1 <- bim$a1[bimIdx[1]]
            allele2 <- bim$a2[bimIdx[1]]

            ## Infer genotypes based on parental genotypes
            if (fatherGeno == 0 && motherGeno == 0) {
                ## Both 0/0 - all children must be 0/0
                for (childId in childIds) {
                    nonInfDf[[paste0(childId, "Pat")]][i] <- allele1
                    nonInfDf[[paste0(childId, "Mat")]][i] <- allele1
                }
                nonInfDf[[paste0(fatherId, "_1")]][i] <- allele1
                nonInfDf[[paste0(fatherId, "_2")]][i] <- allele1
                nonInfDf[[paste0(motherId, "_1")]][i] <- allele1
                nonInfDf[[paste0(motherId, "_2")]][i] <- allele1
            } else if (fatherGeno == 2 && motherGeno == 2) {
                ## Both 2/2 - all children must be 2/2
                for (childId in childIds) {
                    nonInfDf[[paste0(childId, "Pat")]][i] <- allele2
                    nonInfDf[[paste0(childId, "Mat")]][i] <- allele2
                }
                nonInfDf[[paste0(fatherId, "_1")]][i] <- allele2
                nonInfDf[[paste0(fatherId, "_2")]][i] <- allele2
                nonInfDf[[paste0(motherId, "_1")]][i] <- allele2
                nonInfDf[[paste0(motherId, "_2")]][i] <- allele2
            } else if (fatherGeno == 1 && motherGeno == 1) {
                ## Both heterozygous - cannot phase without more info
                ## Leave as NA
            } else {
                ## One parent homozygous, one heterozygous
                ## Use haplotype information from neighboring informative SNPs

                ## Determine which parent is heterozygous
                hetParent <- if (fatherGeno == 1) "father" else "mother"
                homParent <- if (fatherGeno == 1) "mother" else "father"

                hetParentGeno <- if (hetParent == "father") fatherGeno else motherGeno
                homParentGeno <- if (homParent == "father") fatherGeno else motherGeno

                ## Find neighboring informative SNPs
                neighbors <- .findNeighboringSNPs(snpPos, phasedData, phasedPosMap)

                if (!is.null(neighbors$left) || !is.null(neighbors$right)) {
                    ## Infer phase from neighbors for each child
                    for (childId in childIds) {
                        childPatCol <- paste0(childId, "Pat")
                        childMatCol <- paste0(childId, "Mat")

                        ## Get child's genotype at this SNP
                        childColIdx <- which(colnames(genoSubset) == childId)
                        if (length(childColIdx) > 0) {
                            childGeno <- genoSubset[snpId, childColIdx]

                            ## Determine which allele child inherited from het parent
                            inferredAllele <- .inferHetParentAllele(
                                childGeno = childGeno,
                                hetParentGeno = hetParentGeno,
                                homParentGeno = homParentGeno,
                                allele1 = allele1,
                                allele2 = allele2,
                                childId = childId,
                                hetParent = hetParent,
                                neighbors = neighbors,
                                phasedData = phasedData
                            )

                            if (!is.na(inferredAllele)) {
                                ## Assign alleles based on which parent is heterozygous
                                if (hetParent == "father") {
                                    nonInfDf[[childPatCol]][i] <- inferredAllele
                                    ## Mother is homozygous, determine her allele
                                    nonInfDf[[childMatCol]][i] <- if (homParentGeno == 0) allele1 else allele2
                                } else {
                                    nonInfDf[[childMatCol]][i] <- inferredAllele
                                    ## Father is homozygous, determine his allele
                                    nonInfDf[[childPatCol]][i] <- if (homParentGeno == 0) allele1 else allele2
                                }
                            }
                        }
                    }

                    ## Phase the heterozygous parent
                    hetPhase <- .inferHetParentPhase(
                        neighbors = neighbors,
                        phasedData = phasedData,
                        hetParent = hetParent,
                        fatherId = fatherId,
                        motherId = motherId,
                        allele1 = allele1,
                        allele2 = allele2
                    )

                    if (!is.null(hetPhase)) {
                        if (hetParent == "father") {
                            nonInfDf[[paste0(fatherId, "_1")]][i] <- hetPhase$hap1
                            nonInfDf[[paste0(fatherId, "_2")]][i] <- hetPhase$hap2
                            ## Mother is homozygous
                            homAllele <- if (homParentGeno == 0) allele1 else allele2
                            nonInfDf[[paste0(motherId, "_1")]][i] <- homAllele
                            nonInfDf[[paste0(motherId, "_2")]][i] <- homAllele
                        } else {
                            nonInfDf[[paste0(motherId, "_1")]][i] <- hetPhase$hap1
                            nonInfDf[[paste0(motherId, "_2")]][i] <- hetPhase$hap2
                            ## Father is homozygous
                            homAllele <- if (homParentGeno == 0) allele1 else allele2
                            nonInfDf[[paste0(fatherId, "_1")]][i] <- homAllele
                            nonInfDf[[paste0(fatherId, "_2")]][i] <- homAllele
                        }
                    }
                }
            }
        }
    }

    ## Combine informative and non-informative SNPs
    combinedDf <- rbind(phasedData, nonInfDf)

    ## Sort by chromosome and position
    combinedDf$location <- as.numeric(combinedDf$location)

    ## Get chromosome information for sorting
    combinedDf$chr <- NA_integer_
    for (i in seq_len(nrow(combinedDf))) {
        snpId <- as.character(combinedDf$rsID[i])
        bimIdx <- which(bim$snp == snpId)
        if (length(bimIdx) > 0) {
            combinedDf$chr[i] <- bim$chr[bimIdx[1]]
        }
    }

    ## Sort by chromosome then position
    combinedDf <- combinedDf[order(combinedDf$chr, combinedDf$location), ]

    ## Remove temporary chr column
    combinedDf$chr <- NULL

    ## Reset row names
    rownames(combinedDf) <- NULL

    combinedDf
}

#' Find Neighboring Informative SNPs
#'
#' Finds the closest informative SNPs on either side of a target position
#'
#' @param targetPos Physical position of target SNP
#' @param phasedData Data frame with phased informative SNPs
#' @param phasedPosMap Named vector mapping positions to row indices
#'
#' @return List with 'left' and 'right' neighbor indices, or NULL if not found
#' @keywords internal
.findNeighboringSNPs <- function(targetPos, phasedData, phasedPosMap) {
    if (nrow(phasedData) == 0) {
        return(list(left = NULL, right = NULL))
    }

    phasedPositions <- as.numeric(phasedData$location)

    ## Find closest SNP to the left (smaller position)
    leftPositions <- phasedPositions[phasedPositions < targetPos]
    leftIdx <- NULL
    if (length(leftPositions) > 0) {
        leftPos <- max(leftPositions)
        leftIdx <- which(phasedPositions == leftPos)[1]
    }

    ## Find closest SNP to the right (larger position)
    rightPositions <- phasedPositions[phasedPositions > targetPos]
    rightIdx <- NULL
    if (length(rightPositions) > 0) {
        rightPos <- min(rightPositions)
        rightIdx <- which(phasedPositions == rightPos)[1]
    }

    list(left = leftIdx, right = rightIdx)
}

#' Infer Heterozygous Parent Allele for Child
#'
#' Infers which allele a child inherited from the heterozygous parent
#' using neighboring informative SNPs
#'
#' @param childGeno Child's genotype (0, 1, or 2)
#' @param hetParentGeno Heterozygous parent's genotype (always 1)
#' @param homParentGeno Homozygous parent's genotype (0 or 2)
#' @param allele1 Reference allele
#' @param allele2 Alternate allele
#' @param childId Child identifier
#' @param hetParent "father" or "mother" - which parent is heterozygous
#' @param neighbors List with left and right neighbor indices
#' @param phasedData Data frame with phased informative SNPs
#'
#' @return Inferred allele (allele1 or allele2), or NA if cannot infer
#' @keywords internal
.inferHetParentAllele <- function(childGeno, hetParentGeno, homParentGeno,
                                  allele1, allele2, childId, hetParent,
                                  neighbors, phasedData) {
    ## Determine which haplotype column to examine based on parent
    hapCol <- if (hetParent == "father") {
        paste0(childId, "Pat")
    } else {
        paste0(childId, "Mat")
    }

    ## Try to infer from nearest neighbor (prefer closer one)
    neighborIdxList <- c(neighbors$right, neighbors$left)
    neighborIdxList <- neighborIdxList[!vapply(neighborIdxList, is.null, logical(1))]

    if (length(neighborIdxList) == 0) {
        return(NA_character_)
    }

    ## Check neighbors in order of proximity
    for (neighborIdx in neighborIdxList) {
        neighborAllele <- phasedData[[hapCol]][neighborIdx]

        if (!is.na(neighborAllele)) {
            ## Assume haplotype continuity: same allele at nearby SNP
            ## means likely same haplotype

            ## If child is heterozygous (genotype 1), both alleles present
            ## Use neighbor to determine which came from het parent
            if (childGeno == 1) {
                ## Child has one copy of each allele
                ## The one matching the neighbor is likely from het parent
                return(neighborAllele)
            } else if (childGeno == 0) {
                ## Child is 0/0, so inherited allele1 from both parents
                ## If neighbor shows allele1, consistent with het parent giving allele1
                if (neighborAllele == allele1) {
                    return(allele1)
                }
            } else if (childGeno == 2) {
                ## Child is 2/2, so inherited allele2 from both parents
                ## If neighbor shows allele2, consistent with het parent giving allele2
                if (neighborAllele == allele2) {
                    return(allele2)
                }
            }
        }
    }

    return(NA_character_)
}

#' Infer Heterozygous Parent Phase
#'
#' Infers the phase of the heterozygous parent using neighboring informative SNPs
#'
#' @param neighbors List with left and right neighbor indices
#' @param phasedData Data frame with phased informative SNPs
#' @param hetParent "father" or "mother" - which parent is heterozygous
#' @param fatherId Father identifier
#' @param motherId Mother identifier
#' @param allele1 Reference allele
#' @param allele2 Alternate allele
#'
#' @return List with hap1 and hap2, or NULL if cannot infer
#' @keywords internal
.inferHetParentPhase <- function(neighbors, phasedData, hetParent,
                                 fatherId, motherId, allele1, allele2) {
    parentId <- if (hetParent == "father") fatherId else motherId
    hap1Col <- paste0(parentId, "_1")
    hap2Col <- paste0(parentId, "_2")

    ## Try to use nearest neighbor
    neighborIdxList <- c(neighbors$right, neighbors$left)
    neighborIdxList <- neighborIdxList[!vapply(neighborIdxList, is.null, logical(1))]

    if (length(neighborIdxList) == 0) {
        return(NULL)
    }

    ## Use first available neighbor
    for (neighborIdx in neighborIdxList) {
        hap1 <- phasedData[[hap1Col]][neighborIdx]
        hap2 <- phasedData[[hap2Col]][neighborIdx]

        if (!is.na(hap1) && !is.na(hap2)) {
            ## Assume haplotype continuity
            ## If neighbor has allele1 on haplotype 1, current SNP likely does too
            return(list(hap1 = hap1, hap2 = hap2))
        }
    }

    ## If no neighbor information, use simple assignment
    ## Arbitrarily assign allele1 to haplotype 1, allele2 to haplotype 2
    return(NULL)
}

#' Get Informative Genotypes for Phasing
#'
#' Filters for informative SNPs and checks Mendelian consistency
#'
#' @param genoSubset Genotype matrix for one family
#'
#' @return Matrix of informative genotypes or NULL
#' @keywords internal
.getInformativeGenotypes <- function(genoSubset) {
    ## Remove missing data (genotype == 0 or NA)
    genoNoMissing <- genoSubset[apply(genoSubset, 1, function(row) {
        all(row != 0 & !is.na(row))
    }), , drop = FALSE]

    if (nrow(genoNoMissing) == 0L) {
        return(NULL)
    }

    ## Informative SNP configurations: (0,1), (1,0), (2,1), (1,2)
    ## In numeric coding: heterozygous/homozygous parent combinations
    infSnpTypes <- list(c(0, 1), c(1, 0), c(2, 1), c(1, 2))

    infIndex <- vapply(seq_len(nrow(genoNoMissing)), function(i) {
        parentGeno <- as.vector(genoNoMissing[i, seq_len(2)])
        any(vapply(infSnpTypes, function(x) {
            identical(x, parentGeno)
        }, logical(1)))
    }, logical(1))

    if (!any(infIndex)) {
        return(NULL)
    }

    genoInf <- genoNoMissing[infIndex, , drop = FALSE]

    if (nrow(genoInf) == 0L) {
        return(NULL)
    }

    ## Filter for Mendelian consistency
    mendel <- vapply(seq_len(nrow(genoInf)), function(i) {
        all(genoInf[i, 3:ncol(genoInf)] %in% genoInf[i, seq_len(2)])
    }, logical(1))

    if (!any(mendel)) {
        return(NULL)
    }

    genoAll <- genoInf[mendel, , drop = FALSE]

    if (nrow(genoAll) == 0L) {
        return(NULL)
    }

    genoAll
}

#' Assign Allele Codes for One Parent
#'
#' Determines which allele (1 or 2) was inherited from specified parent
#'
#' @param genoDf Informative genotype matrix
#' @param parent "father" or "mother"
#' @param fatherCol Column index for father
#' @param motherCol Column index for mother
#'
#' @return Data frame with allele codes
#' @keywords internal
.assignAlleleCodes <- function(genoDf, parent, fatherCol, motherCol) {
    parentCol <- if (parent == "father") fatherCol else motherCol

    ## Select SNPs where parent is heterozygous (genotype = 1)
    genoParent <- genoDf[genoDf[, parentCol] == 1, , drop = FALSE]

    if (nrow(genoParent) == 0L) {
        return(NULL)
    }

    ## Determine allele codes for children
    genoChildren <- genoParent[, -c(fatherCol, motherCol), drop = FALSE]
    nChildren <- ncol(genoChildren)

    codeMatrix <- matrix(0L, nrow = nrow(genoParent), ncol = nChildren)

    for (i in seq_len(nrow(genoParent))) {
        fatherGt <- genoParent[i, fatherCol]
        motherGt <- genoParent[i, motherCol]

        for (j in seq_len(nChildren)) {
            childGt <- genoChildren[i, j]

            codeMatrix[i, j] <- .matchAlleleCode(
                parent, fatherGt, motherGt, childGt
            )
        }
    }

    ## Create data frame with SNP IDs
    result <- data.frame(
        rsID = rownames(genoParent),
        codeMatrix,
        stringsAsFactors = FALSE
    )
    colnames(result) <- c("rsID", colnames(genoChildren))
    rownames(result) <- rownames(genoParent)

    result
}

#' Match Allele Code for Single Genotype
#'
#' @param parent "father" or "mother"
#' @param fatherGt Father genotype (0, 1, or 2)
#' @param motherGt Mother genotype (0, 1, or 2)
#' @param childGt Child genotype (0, 1, or 2)
#'
#' @return Allele code (1 or 2)
#' @keywords internal
.matchAlleleCode <- function(parent, fatherGt, motherGt, childGt) {
    if (parent == "father") {
        ## Father is heterozygous (1), determine which allele passed to child
        if (fatherGt == 1 && motherGt == 2 && childGt == 2) {
            return(2L)
        } else if (fatherGt == 1 && motherGt == 2 && childGt == 1) {
            return(1L)
        } else if (fatherGt == 1 && motherGt == 0 && childGt == 0) {
            return(1L)
        } else if (fatherGt == 1 && motherGt == 0 && childGt == 1) {
            return(2L)
        }
    } else {
        ## Mother is heterozygous (1)
        if (fatherGt == 2 && motherGt == 1 && childGt == 2) {
            return(2L)
        } else if (fatherGt == 2 && motherGt == 1 && childGt == 1) {
            return(1L)
        } else if (fatherGt == 0 && motherGt == 1 && childGt == 0) {
            return(1L)
        } else if (fatherGt == 0 && motherGt == 1 && childGt == 1) {
            return(2L)
        }
    }

    ## If no matches
    NA_integer_
}

#' Convert Allele Codes to Nucleotides
#'
#' @param alleleCodeDf Data frame with allele codes (1 or 2)
#' @param bim BIM data frame with allele information
#'
#' @return Data frame with nucleotide alleles
#' @keywords internal
.convertToNucleotides <- function(alleleCodeDf, bim) {
    if (is.null(alleleCodeDf)) {
        return(NULL)
    }

    nRows <- nrow(alleleCodeDf)
    nCols <- ncol(alleleCodeDf) + 1 ## Add location column

    actgMatrix <- matrix(character(nRows * nCols),
        nrow = nRows, ncol = nCols
    )

    ## Add SNP IDs and locations
    actgMatrix[, 1] <- as.character(alleleCodeDf$rsID)

    for (i in seq_len(nRows)) {
        snpId <- alleleCodeDf$rsID[i]
        bimIdx <- which(bim$snp == snpId)

        if (length(bimIdx) > 0L) {
            actgMatrix[i, 2] <- as.character(bim$pos[bimIdx[1]])

            for (j in 2:ncol(alleleCodeDf)) {
                code <- alleleCodeDf[i, j]
                if (code == 1) {
                    actgMatrix[i, j + 1] <- bim$a1[bimIdx[1]]
                } else if (code == 2) {
                    actgMatrix[i, j + 1] <- bim$a2[bimIdx[1]]
                }
            }
        }
    }

    result <- as.data.frame(actgMatrix, stringsAsFactors = FALSE)
    colnames(result) <- c("rsID", "location", colnames(alleleCodeDf)[-1])
    rownames(result) <- rownames(alleleCodeDf)

    result
}

#' Add Homozygous SNPs from Other Parent
#'
#' @param actgMatrixHom Matrix for parent being imputed
#' @param actgMatrixHet Matrix from other parent (het SNPs)
#' @param genotypeDf Original genotype data frame
#' @param bim BIM data frame
#' @param parent "father" or "mother"
#'
#' @return Combined data frame with het and hom SNPs
#' @keywords internal
.addHomozygousSNPs <- function(actgMatrixHom, actgMatrixHet,
                               genotypeDf, bim, parent) {
    if (is.null(actgMatrixHom) || is.null(actgMatrixHet)) {
        return(actgMatrixHom)
    }

    pCol <- if (parent == "father") 1L else 2L
    nSamples <- ncol(actgMatrixHom) - 2L

    ## Extract homozygous alleles from heterozygous parent's SNPs
    homList <- lapply(rownames(actgMatrixHet), function(snp) {
        bimIdx <- which(bim$snp == snp)
        genoIdx <- which(rownames(genotypeDf) == snp)

        if (length(bimIdx) == 0L || length(genoIdx) == 0L) {
            return(NULL)
        }

        ## Determine homozygous allele
        if (genotypeDf[genoIdx, pCol] == 0) {
            homAllele <- bim$a1[bimIdx[1]]
        } else {
            homAllele <- bim$a2[bimIdx[1]]
        }

        ## Create row with homozygous allele for all samples
        data.frame(
            rsID = snp,
            location = as.character(bim$pos[bimIdx[1]]),
            matrix(homAllele, nrow = 1, ncol = nSamples),
            stringsAsFactors = FALSE
        )
    })

    homDf <- do.call(rbind, homList[!vapply(homList, is.null, logical(1))])

    if (is.null(homDf) || nrow(homDf) == 0L) {
        actgMatrixHom$hom <- 1
        return(actgMatrixHom)
    }

    colnames(homDf) <- colnames(actgMatrixHom)

    ## Mark het vs hom SNPs
    actgMatrixHom$hom <- 1
    homDf$hom <- 0

    ## Combine and sort by location
    combinedDf <- rbind(actgMatrixHom, homDf)
    combinedDf$location <- as.numeric(combinedDf$location)
    combinedDf <- combinedDf[order(combinedDf$location), ]

    combinedDf
}

#' Phase Child Haplotypes
#'
#' @param actgCombPat Paternal allele data frame
#' @param actgCombMat Maternal allele data frame
#'
#' @return Data frame with phased children
#' @keywords internal
.phaseChildren <- function(actgCombPat, actgCombMat) {
    ## Ensure both data frames have the same SNPs
    ## Union of all SNPs
    allSnps <- unique(c(
        as.character(actgCombPat$rsID),
        as.character(actgCombMat$rsID)
    ))

    ## Create aligned data frames with all SNPs
    actgCombPatAligned <- .alignSnpData(actgCombPat, allSnps)
    actgCombMatAligned <- .alignSnpData(actgCombMat, allSnps)

    nRows <- nrow(actgCombPatAligned)
    nChildren <- (ncol(actgCombPatAligned) - 3) ## Subtract rsID, location, hom
    nCols <- 2 + 2 * nChildren ## rsID, location, 2 haplotypes per child

    phaseChildren <- matrix(character(nRows * nCols),
        nrow = nRows, ncol = nCols
    )

    ## Add SNP IDs and locations
    phaseChildren[, 1] <- as.character(actgCombPatAligned$rsID)
    phaseChildren[, 2] <- as.character(actgCombPatAligned$location)

    ## Add paternal and maternal haplotypes for each child
    nameVec <- c("rsID", "location")

    for (j in seq_len(nChildren)) {
        colIdx <- j + 2 ## Skip rsID and location
        phaseChildren[, 2 + (j - 1) * 2 + 1] <- as.character(actgCombPatAligned[, colIdx])
        phaseChildren[, 2 + (j - 1) * 2 + 2] <- as.character(actgCombMatAligned[, colIdx])

        childName <- colnames(actgCombPatAligned)[colIdx]
        nameVec <- c(
            nameVec, paste0(childName, "Pat"),
            paste0(childName, "Mat")
        )
    }

    result <- as.data.frame(phaseChildren, stringsAsFactors = FALSE)
    colnames(result) <- nameVec

    result
}

#' Align SNP Data to Common Set
#'
#' Ensures a data frame contains all specified SNPs in order, filling missing
#' SNPs with NA values
#'
#' @param snpDf Data frame with SNP information (must have rsID column)
#' @param allSnps Character vector of all SNPs to include
#'
#' @return Data frame with all SNPs aligned
#' @keywords internal
.alignSnpData <- function(snpDf, allSnps) {
    if (nrow(snpDf) == 0 || length(allSnps) == 0) {
        return(snpDf)
    }

    ## Identify missing SNPs
    existingSnps <- as.character(snpDf$rsID)
    missingSnps <- setdiff(allSnps, existingSnps)

    if (length(missingSnps) == 0) {
        ## All SNPs present; reorder
        snpDf <- snpDf[match(allSnps, snpDf$rsID), , drop = FALSE]
        return(snpDf)
    }

    ## Template for missing SNPs
    nCols <- ncol(snpDf)
    nMissing <- length(missingSnps)

    missingDf <- as.data.frame(
        matrix(NA_character_,
            nrow = nMissing,
            ncol = nCols
        ),
        stringsAsFactors = FALSE
    )
    colnames(missingDf) <- colnames(snpDf)
    missingDf$rsID <- missingSnps

    ## Combine existing and missing data
    combinedDf <- rbind(snpDf, missingDf)

    ## Reorder to match allSnps
    combinedDf <- combinedDf[match(allSnps, combinedDf$rsID), , drop = FALSE]

    ## Ensure location is numeric for sorting
    if ("location" %in% colnames(combinedDf)) {
        combinedDf$location <- as.numeric(combinedDf$location)
    }

    combinedDf
}

#' Phase Parent Haplotypes Using Crossover Information
#'
#' @param actgMatrix Allele matrix for parent
#' @param xoDetectMatrix Crossover detection results
#' @param bim BIM data frame
#' @param parentId Parent identifier
#'
#' @return Data frame with two phased haplotypes
#' @keywords internal
.phaseParent <- function(actgMatrix, xoDetectMatrix, bim, parentId) {
    nSnps <- nrow(actgMatrix)

    ## Find a child without crossovers (if any)
    childCols <- setdiff(
        colnames(actgMatrix),
        c("rsID", "location", "hom")
    )

    ## Check if any children have no crossovers
    childrenNoXo <- setdiff(childCols, xoDetectMatrix$CHILD)

    if (length(childrenNoXo) > 0) {
        ## Use first child without crossovers
        parentL <- as.character(actgMatrix[[childrenNoXo[1]]])
        names(parentL) <- actgMatrix$rsID

        ## Other haplotype is complementary allele
        parentR <- vapply(seq_along(parentL), function(i) {
            .getOtherAllele(names(parentL)[i], parentL[i], bim)
        }, character(1))

        ## For homozygous sites, both haplotypes are the same
        homSites <- which(actgMatrix$hom == 0)
        parentR[homSites] <- parentL[homSites]

        names(parentR) <- names(parentL)
    } else {
        ## Use crossover information to phase
        childName <- childCols[1]
        parentL <- as.character(actgMatrix[[childName]])
        names(parentL) <- actgMatrix$rsID

        ## Get crossover positions for this child
        xoChild <- xoDetectMatrix[xoDetectMatrix$CHILD == childName, ]

        if (nrow(xoChild) > 0) {
            ## Get indices of crossover SNPs
            xoSnps <- xoChild$START
            xoIndices <- match(xoSnps, names(parentL))
            xoIndices <- xoIndices[!is.na(xoIndices)]

            ## Switch haplotypes at crossover points
            if (length(xoIndices) > 0) {
                for (k in seq_along(xoIndices)) {
                    if (k %% 2 == 1) { ## Odd crossovers - switch
                        endIdx <- if (k < length(xoIndices)) {
                            xoIndices[k + 1] - 1
                        } else {
                            length(parentL)
                        }

                        ## Get complementary alleles for this segment
                        for (j in xoIndices[k]:endIdx) {
                            parentL[j] <- .getOtherAllele(
                                names(parentL)[j],
                                parentL[j],
                                bim
                            )
                        }
                    }
                }
            }
        }

        ## Other haplotype is complementary
        parentR <- vapply(seq_along(parentL), function(i) {
            .getOtherAllele(names(parentL)[i], parentL[i], bim)
        }, character(1))

        ## For homozygous sites, both are the same
        homSites <- which(actgMatrix$hom == 0)
        parentR[homSites] <- parentL[homSites]

        names(parentR) <- names(parentL)
    }

    ## Create output with parent ID column names
    result <- data.frame(
        parentL,
        parentR,
        stringsAsFactors = FALSE
    )
    colnames(result) <- c(paste0(parentId, "_1"), paste0(parentId, "_2"))
    rownames(result) <- names(parentL)

    result
}

#' Get Complementary Allele
#'
#' @param snpId SNP identifier
#' @param allele Current allele (A, C, G, or T)
#' @param bim BIM data frame
#'
#' @return Complementary allele
#' @keywords internal
.getOtherAllele <- function(snpId, allele, bim) {
    bimIdx <- which(bim$snp == snpId)

    if (length(bimIdx) == 0L) {
        return(NA_character_)
    }

    if (bim$a1[bimIdx[1]] == allele) {
        return(bim$a2[bimIdx[1]])
    } else {
        return(bim$a1[bimIdx[1]])
    }
}

#' Get First Non-NA Value
#'
#' @param x Vector
#'
#' @return First non-NA value
#' @keywords internal
.firstNonNA <- function(x) {
    x[!is.na(x)][1]
}

#' Convert Phased List to SummarizedExperiment
#'
#' @param phaseList Named list of phased DataFrames
#' @param bim BIM data frame
#'
#' @return SummarizedExperiment object
#' @keywords internal
.convertToSummarizedExperiment <- function(phaseList, bim) {
    ## Keep only valid phased DataFrames
    validFamilies <- phaseList[vapply(
        phaseList,
        function(x) is(x, "DataFrame") || is.data.frame(x),
        logical(1)
    )]

    if (!length(validFamilies)) {
        stop("No valid phased data to convert to SummarizedExperiment")
    }

    ## Collect all SNPs and samples
    allSnps <- unique(unlist(lapply(validFamilies, function(fam) {
        fam$rsID
    })))

    allSamples <- unique(unlist(lapply(validFamilies, function(fam) {
        setdiff(colnames(fam), c("rsID", "location"))
    })))

    ## Initialize genotype matrix
    genoMatrix <- matrix(
        NA_character_,
        nrow = length(allSnps),
        ncol = length(allSamples),
        dimnames = list(allSnps, allSamples)
    )

    ## Fill genotype matrix family-wise
    for (famData in validFamilies) {
        sampleCols <- setdiff(colnames(famData), c("rsID", "location"))
        if (!length(sampleCols)) next

        rowIdx <- match(famData$rsID, allSnps)
        colIdx <- match(sampleCols, allSamples)

        genoMatrix[rowIdx, colIdx] <- as.matrix(famData[, sampleCols, drop = FALSE])
    }

    ## Build rowData via BIM lookup
    bimIdx <- match(allSnps, bim$snp)

    rowDataDf <- data.frame(
        snpId = allSnps,
        chromosome = bim$chr[bimIdx],
        position = bim$pos[bimIdx],
        allele1 = bim$a1[bimIdx],
        allele2 = bim$a2[bimIdx],
        stringsAsFactors = FALSE
    )

    ## Build colData
    colDataDf <- data.frame(
        sampleId = allSamples,
        stringsAsFactors = FALSE
    )

    colDataDf$sampleType <- ifelse(
        grepl("Pat$|Mat$", allSamples),
        "child_haplotype",
        ifelse(grepl("_1$|_2$", allSamples),
            "parent_haplotype",
            "unknown"
        )
    )

    colDataDf$individualId <- gsub("(Pat|Mat|_1|_2)$", "", allSamples)

    ## Map samples to families
    sampleToFamily <- setNames(
        rep(NA_character_, length(allSamples)),
        allSamples
    )

    for (famId in names(validFamilies)) {
        famSamples <- setdiff(
            colnames(validFamilies[[famId]]),
            c("rsID", "location")
        )
        sampleToFamily[famSamples] <- famId
    }

    colDataDf$familyId <- unname(sampleToFamily[allSamples])

    ## Create SummarizedExperiment
    SummarizedExperiment::SummarizedExperiment(
        assays  = list(phased = genoMatrix),
        rowData = rowDataDf,
        colData = colDataDf
    )
}

#' Write Phased Genotypes to VCF Files
#'
#' @param phaseList Named list of phased DataFrames
#' @param bim BIM data frame
#' @param vcfOutput Path prefix for VCF output
#'
#' @return Character vector of created VCF file paths
#' @keywords internal
.writeVcfFiles <- function(phaseList, bim, vcfOutput) {
    ## Keep only valid phased DataFrames
    validFamilies <- phaseList[vapply(
        phaseList,
        function(x) is(x, "DataFrame") || is.data.frame(x),
        logical(1)
    )]

    if (!length(validFamilies)) {
        stop("No valid phased data to write to VCF")
    }

    ## Precompute BIM lookup
    bimIdx <- match(bim$snp, bim$snp)
    bimMap <- data.frame(
        snp = bim$snp,
        chr = bim$chr,
        pos = bim$pos,
        ref = bim$a1,
        alt = bim$a2,
        stringsAsFactors = FALSE
    )
    rownames(bimMap) <- bimMap$snp

    vcfFiles <- character(length(validFamilies))

    for (i in seq_along(validFamilies)) {
        famId <- names(validFamilies)[i]
        famData <- validFamilies[[i]]

        vcfFile <- paste0(vcfOutput, "_", famId, ".vcf")
        vcfFiles[i] <- vcfFile

        sampleCols <- setdiff(colnames(famData), c("rsID", "location"))

        ## Header
        header <- c(
            "##fileformat=VCFv4.2",
            paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
            "##source=xoPhase",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased Genotype\">",
            paste0(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
                paste(sampleCols[grepl("Pat$|_1$", sampleCols)], collapse = "\t")
            )
        )

        writeLines(header, vcfFile)

        ## SNP lookup
        snpIds <- as.character(famData$rsID)
        keep <- snpIds %in% rownames(bimMap)

        famData <- famData[keep, , drop = FALSE]
        snpIds <- snpIds[keep]

        bimSub <- bimMap[snpIds, ]

        ## Genotype matrix
        geno <- as.matrix(famData[, sampleCols, drop = FALSE])

        ## Allele to numeric conversion
        ref <- bimSub$ref
        alt <- bimSub$alt

        genoNum <- matrix(
            ".",
            nrow = nrow(geno),
            ncol = ncol(geno)
        )

        genoNum[geno == ref] <- "0"
        genoNum[geno == alt] <- "1"

        ## Collapse haplotypes into phased genotypes
        patCols <- grepl("Pat$|_1$", sampleCols)
        matCols <- grepl("Mat$|_2$", sampleCols)

        phased <- matrix(
            paste0(
                genoNum[, patCols, drop = FALSE],
                "|",
                genoNum[, matCols, drop = FALSE]
            ),
            nrow = nrow(genoNum)
        )

        genoStr <- do.call(
            paste,
            c(as.data.frame(phased), sep = "\t")
        )

        vcfLines <- paste(
            bimSub$chr,
            bimSub$pos,
            snpIds,
            ref,
            alt,
            ".",
            "PASS",
            ".",
            "GT",
            genoStr,
            sep = "\t"
        )

        ## Append data lines to VCF
        write.table(vcfLines, vcfFile,
            append = TRUE, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t"
        )
    }

    message("VCF files written:")
    message(paste("  ", vcfFiles, collapse = "\n"))

    invisible(vcfFiles)
}

#' Convert GRanges to DataFrame
#'
#' Converts a GRanges object to a data.frame/DataFrame for internal processing
#'
#' @param granges GRanges object
#'
#' @return data.frame with genomic ranges information and metadata columns
#' @keywords internal
.granges2DataFrame <- function(granges) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        stop("Package 'GenomicRanges' is required")
    }

    if (!is(granges, "GRanges")) {
        stop("Input must be a GRanges object")
    }

    if (length(granges) == 0) {
        warning("Empty GRanges object provided")
        return(data.frame())
    }

    ## Extract metadata columns first
    mcolsData <- GenomicRanges::mcols(granges)

    if (!is.null(mcolsData) && ncol(mcolsData) > 0) {
        ## Convert mcols to data.frame
        df <- as.data.frame(mcolsData, stringsAsFactors = FALSE)

        ## Add genomic location columns (prepend to maintain order)
        locationDf <- data.frame(
            seqnames = as.character(GenomicRanges::seqnames(granges)),
            start = GenomicRanges::start(granges),
            end = GenomicRanges::end(granges),
            width = GenomicRanges::width(granges),
            strand = as.character(GenomicRanges::strand(granges)),
            stringsAsFactors = FALSE
        )

        ## Combine location and metadata
        df <- cbind(locationDf, df)
    } else {
        ## Basic data frame with genomic coordinates (no metadata)
        df <- data.frame(
            seqnames = as.character(GenomicRanges::seqnames(granges)),
            start = GenomicRanges::start(granges),
            end = GenomicRanges::end(granges),
            width = GenomicRanges::width(granges),
            strand = as.character(GenomicRanges::strand(granges)),
            stringsAsFactors = FALSE
        )
    }

    ## Set row names if the GRanges object has names
    if (!is.null(names(granges)) && length(names(granges)) == nrow(df)) {
        rownames(df) <- names(granges)
    }

    df
}
