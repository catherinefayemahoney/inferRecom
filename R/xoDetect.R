#' Detect Meiotic Recombinations from PLINK Genotype Data
#'
#' This function identifies recombination events in family-based genotype data,
#' using PLINK-formatted input files and a genetic map. It supports 2- or
#' 3-child family structures, and can optionally write results to disk if an
#' output path is provided.
#'
#' @param plinkFile Character string; path to the PLINK prefix (without
#'     extensions).
#' @param mapFile Character string; path to a tab-delimited map file containing
#'     columns: \code{chr}, \code{pos}, \code{rate}, \code{cM}.
#' @param familySize Integer; either 2 or 3, specifying the family pedigree size.
#' @param parent Character; either \code{"mother"} or \code{"father"},
#'     specifying the parent to trace crossovers from.
#' @param rsOnly Logical; whether to restrict to SNPs with names beginning
#'     with "rs". Default is \code{FALSE}.
#' @param snpFilter Integer; minimum number of SNPs separating putative
#'     crossovers. Default is 5.
#' @param cmFilter Numeric; minimum genetic distance (cM) separating putative
#'     crossovers. Default is 1.
#' @param caseControl Logical; if TRUE, function returns a
#'    \link[GenomicRanges:GRangesList-class]{GRangesList} of case and control
#'    crossovers for 3+ child families. Default is \code{FALSE}.
#' @param out Character or NULL; if provided, output is written as a CSV at this
#'     path.
#' @param BPPARAM A BiocParallelParam object specifying parallel execution.
#'     Default is \code{SerialParam()} for serial execution. Use
#'     \code{MulticoreParam(workers = n)} for parallel processing across n cores.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} listing detected
#'     crossover intervals with metadata columns: \code{childId} (in 3-child
#'     families), \code{familyId}, \code{startSnp}, \code{finishSnp},
#'     \code{startPos}, \code{finishPos}, \code{startCm}, \code{finishCm}.
#'
#' @details
#' This function reads PLINK genotype data (via \pkg{snpStats}) and infers
#' loci where meiotic recombination has occurred, resolved to the nearest
#' informative SNPs. The algorithm identifies informative marker configurations
#' where parents are heterozygous/homozygous and filters for Mendelian
#' consistency before detecting inheritance state changes that indicate
#' crossovers.
#'
#' For 3-child families, children are analyzed in triples to allow
#' identification of the specific recombinant child. For 2-child families, state
#' changes between the two siblings are used to infer crossover locations, but
#' the individual cannot be determined.
#'
#' If \code{out} is provided, the resulting table is written to disk as
#' UTF-8 encoded CSV without row names.
#'
#' Parallel processing is performed by family, improving efficiency for datasets
#' with many families.
#'
#' @import methods
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges GRangesList seqnames start makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom snpStats read.plink
#' @importFrom utils read.delim write.csv
#' @importFrom dplyr distinct
#' @importFrom BiocParallel SerialParam MulticoreParam bplapply
#' @importFrom stats var
#'
#' @examples
#' # Serial execution (default)
#'
#' # Load example data
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#'
#' # Basic maternal crossover detection
#' mapFemale <- file.path(dataPath, "female_chr4.txt")
#' xoMat3 <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 3,
#'     parent = "mother"
#' )
#'
#' # View results
#' xoMat3
#'
#' # Basic 2-child family
#' xoMat2 <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 2,
#'     parent = "mother"
#' )
#'
#' # View results
#' xoMat2
#'
#' # Maternal crossover detection with case/control separation
#' mapFemale <- file.path(dataPath, "female_chr4.txt")
#' xoMat3CC <- xoDetect(
#'     plinkFile = plinkFile,
#'     mapFile = mapFemale,
#'     familySize = 3,
#'     parent = "mother",
#'     caseControl = TRUE
#' )
#'
#' # View results
#' xoMat3CC$case
#' xoMat3CC$control
#'
#' @export
xoDetect <- function(plinkFile,
                     mapFile,
                     familySize = c(2, 3),
                     parent = c("mother", "father"),
                     rsOnly = FALSE,
                     snpFilter = 5,
                     cmFilter = 1,
                     caseControl = FALSE,
                     out = NULL,
                     BPPARAM = SerialParam()) {
    ## Argument validation
    parent <- match.arg(parent)
    familySize <- as.integer(match.arg(as.character(familySize), c("2", "3")))

    bedFile <- paste0(plinkFile, ".bed")
    bimFile <- paste0(plinkFile, ".bim")
    famFile <- paste0(plinkFile, ".fam")

    missingFiles <- character()
    if (!file.exists(bedFile)) {
        missingFiles <- c(missingFiles, bedFile)
    }
    if (!file.exists(bimFile)) {
        missingFiles <- c(missingFiles, bimFile)
    }
    if (!file.exists(famFile)) {
        missingFiles <- c(missingFiles, famFile)
    }

    if (length(missingFiles) > 0L) {
        stop("PLINK file(s) not found: ", paste(missingFiles, collapse = ", "))
    }

    if (!file.exists(mapFile)) {
        stop("Genetic map file not found: ", mapFile)
    }

    if (!is.null(out) && !is.character(out)) {
        stop("'out' must be a character string or NULL")
    }

    if (!is.logical(rsOnly)) {
        stop("'rsOnly' must be logical")
    }

    if (!is.numeric(snpFilter) || snpFilter < 0) {
        stop("'snpFilter' must be a non-negative integer")
    }

    if (!is.numeric(cmFilter) || cmFilter < 0) {
        stop("'cmFilter' must be a non-negative numeric value")
    }

    if (familySize == 2 && caseControl == TRUE) {
        stop("No case-control separation for 2-child families")
    }

    ## Dependency checks
    pkgs <- c("snpStats", "S4Vectors", "GenomicRanges", "IRanges", "dplyr", "BiocParallel")
    miss <- pkgs[!vapply(pkgs, requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )]
    if (length(miss) > 0L) {
        stop("Missing required packages: ", paste(miss, collapse = ", "))
    }

    ## Strip extension from plinkFile if present
    plinkFile <- sub("\\.(bed|bim|fam)$", "", plinkFile, ignore.case = TRUE)

    ## Read PLINK data
    pedData <- snpStats::read.plink(plinkFile)
    mapDf <- utils::read.delim(mapFile, stringsAsFactors = FALSE)

    if (!all(c("chr", "pos", "cM") %in% colnames(mapDf))) {
        stop("mapFile must have columns: chr, pos, cM")
    }

    ## Extract genotype matrix
    genoMatrix <- as.data.frame(t(pedData$genotypes))

    ## Filter for rs SNPs if requested
    if (isTRUE(rsOnly)) {
        rsIdx <- grep("^rs", rownames(genoMatrix), ignore.case = TRUE)
        if (length(rsIdx) == 0L) {
            stop("No SNPs with 'rs' prefix found")
        }
        genoMatrix <- genoMatrix[rsIdx, , drop = FALSE]
    }

    famData <- pedData$fam

    ## Dispatch to appropriate helper function
    if (familySize == 3L) {
        resDf <- .xoDetectThreeChild(
            pedData, genoMatrix, famData, mapDf, parent,
            snpFilter, cmFilter, BPPARAM
        )
    } else {
        resDf <- .xoDetectTwoChild(
            pedData, genoMatrix, famData, mapDf, parent,
            snpFilter, cmFilter, BPPARAM
        )
    }

    ## Write output if requested
    if (!is.null(out) && nrow(resDf) > 0L) {
        utils::write.csv(as.data.frame(resDf),
            file = out,
            row.names = FALSE, fileEncoding = "UTF-8"
        )
        message("Results written to: ", out)
    }

    ## Convert to appropriate GRanges object or list
    if (nrow(resDf) > 0L) {
        # Add chromosome information
        resDf$chromosome <- mapDf$chr[1]

        if (caseControl == FALSE) {
            # Convert to GRanges
            grResult <- GenomicRanges::makeGRangesFromDataFrame(
                df = resDf,
                seqnames.field = "chromosome",
                start.field = "startPos",
                end.field = "finishPos",
                keep.extra.columns = TRUE
            )
        } else {
            # Add case/control status
            affectedIds <- famData$member[famData$affected == 1]
            resDf$caseStatus <- ifelse(resDf$childId %in% affectedIds,
                "case",
                "control"
            )

            # Convert to GRanges first
            grTemp <- GenomicRanges::makeGRangesFromDataFrame(
                df = resDf,
                seqnames.field = "chromosome",
                start.field = "startPos",
                end.field = "finishPos",
                keep.extra.columns = TRUE
            )

            # Split into GRangesList by case/control status
            grResult <- split(grTemp, grTemp$caseStatus)
        }
    } else {
        # Return empty GRanges or GRangesList
        if (caseControl == FALSE) {
            grResult <- GRanges()
        } else {
            grResult <- GRangesList(case = GRanges(), control = GRanges())
        }
    }

    grResult
}


## ========================== INTERNAL HELPERS ===========================

## Extract first non-NA value from vector
.firstNonNA <- function(x) {
    idx <- which(!is.na(x))
    if (length(idx) == 0L) {
        return(NA)
    }
    x[idx[1L]]
}

## Identify informative SNP types (heterozygous/homozygous parent pairs)
.getInformativeSNPs <- function(genoSubset) {
    ## Remove missing data (genotype == 0)
    genoNoMissing <- genoSubset[apply(genoSubset, 1, function(row) {
        all(row != 0)
    }), , drop = FALSE]

    if (nrow(genoNoMissing) == 0L) {
        return(NULL)
    }

    ## Informative SNP configurations: (1,2), (2,1), (3,2), (2,3)
    infSnpTypes <- list(c(1, 2), c(2, 1), c(3, 2), c(2, 3))

    infIndex <- vapply(seq_len(nrow(genoNoMissing)), function(i) {
        parentGeno <- as.vector(genoNoMissing[i, seq_len(2)])
        any(vapply(infSnpTypes, function(x) {
            identical(x, parentGeno)
        }, logical(1)))
    }, logical(1))

    if (!any(infIndex)) {
        return(NULL)
    }

    genoInformative <- genoNoMissing[infIndex, , drop = FALSE]

    if (nrow(genoInformative) == 0L || is.vector(genoInformative)) {
        return(NULL)
    }

    ## Filter for Mendelian consistency
    mendelConsistent <- vapply(seq_len(nrow(genoInformative)), function(i) {
        all(genoInformative[i, 3:ncol(genoInformative)] %in%
            as.vector(genoInformative[i, seq_len(2)]))
    }, logical(1))

    if (!any(mendelConsistent)) {
        return(NULL)
    }

    genoFiltered <- genoInformative[mendelConsistent, , drop = FALSE]

    if (nrow(genoFiltered) == 0L || is.vector(genoFiltered)) {
        return(NULL)
    }

    genoFiltered
}

## Create binary inheritance matrix for parent-specific haplotypes
.makeBinaryInheritance <- function(parentSubset) {
    if (nrow(parentSubset) <= 1L) {
        return(NULL)
    }

    binaryMatrix <- matrix(
        nrow = nrow(parentSubset),
        ncol = (ncol(parentSubset) - 2)
    )

    for (i in seq_len(nrow(parentSubset))) {
        for (j in 3:ncol(parentSubset)) {
            if (parentSubset[i, 2] == parentSubset[i, j]) {
                binaryMatrix[i, j - 2] <- 1
            } else {
                binaryMatrix[i, j - 2] <- 0
            }
        }
    }

    rownames(binaryMatrix) <- rownames(parentSubset)
    binaryMatrix
}

## Apply SNP filter to remove closely-spaced state changes
.applySNPFilter <- function(stateChange, snpFilter) {
    if (snpFilter == 0) {
        return(character())
    }

    snpFilterDrop <- character()
    changeLength <- length(stateChange)

    if (changeLength <= snpFilter) {
        return(character())
    }

    changeIndex <- which(stateChange != 0)
    changeIndex <- changeIndex[changeIndex < (changeLength - snpFilter)]

    for (k in changeIndex) {
        window <- stateChange[k:(k + snpFilter)]
        if (sum(window) == 0 && var(window) != 0) {
            snpFilterDrop <- c(snpFilterDrop, rownames(stateChange)[k])
        }
    }

    snpFilterDrop
}

## Apply cM filter to remove closely-spaced crossovers
.applyCMFilter <- function(xoData, cmFilter) {
    if (is.na(cmFilter) || nrow(xoData) <= 1L) {
        return(xoData)
    }

    ## Sort by child and position
    if ("childId" %in% colnames(xoData)) {
        xoSorted <- xoData[order(xoData$childId, xoData$startPos), ]

        dropIndex <- numeric()
        for (i in 2:nrow(xoSorted)) {
            if (xoSorted$childId[i - 1] == xoSorted$childId[i] &&
                (xoSorted$startCm[i] - xoSorted$finishCm[i - 1] < cmFilter)) {
                dropIndex <- c(dropIndex, i - 1, i)
            }
        }
    } else {
        xoSorted <- xoData[order(xoData$familyId, xoData$startPos), ]

        dropIndex <- numeric()
        for (i in 2:nrow(xoSorted)) {
            if (xoSorted$familyId[i - 1] == xoSorted$familyId[i] &&
                (xoSorted$startCm[i] - xoSorted$finishCm[i - 1] < cmFilter)) {
                dropIndex <- c(dropIndex, i - 1, i)
            }
        }
    }

    dropIndex <- unique(dropIndex)
    if (length(dropIndex) > 0L) {
        xoSorted <- xoSorted[-dropIndex, ]
    }

    xoSorted
}

## Add physical and genetic positions to crossover results
.addPositions <- function(xoData, pedMap, geneticMap) {
    if (nrow(xoData) == 0L) {
        return(xoData)
    }

    startPosList <- numeric(nrow(xoData))
    finishPosList <- numeric(nrow(xoData))
    startCmList <- numeric(nrow(xoData))
    finishCmList <- numeric(nrow(xoData))

    for (rowIdx in seq_len(nrow(xoData))) {
        ## Physical positions
        startSnpIdx <- which(pedMap$snp.name == xoData$startSnp[rowIdx])
        finishSnpIdx <- which(pedMap$snp.name == xoData$finishSnp[rowIdx])

        if (length(startSnpIdx) > 0L) {
            startPosList[rowIdx] <- pedMap$position[startSnpIdx[1]]
        } else {
            startPosList[rowIdx] <- NA_real_
        }

        if (length(finishSnpIdx) > 0L) {
            finishPosList[rowIdx] <- pedMap$position[finishSnpIdx[1]]
        } else {
            finishPosList[rowIdx] <- NA_real_
        }

        ## Genetic positions (cM)
        if (!is.na(startPosList[rowIdx])) {
            startCmList[rowIdx] <- geneticMap$cM[which.min(abs(startPosList[rowIdx] -
                geneticMap$pos))]
        } else {
            startCmList[rowIdx] <- NA_real_
        }

        if (!is.na(finishPosList[rowIdx])) {
            finishCmList[rowIdx] <- geneticMap$cM[which.min(abs(finishPosList[rowIdx] -
                geneticMap$pos))]
        } else {
            finishCmList[rowIdx] <- NA_real_
        }
    }

    xoData$startPos <- startPosList
    xoData$finishPos <- finishPosList
    xoData$startCm <- startCmList
    xoData$finishCm <- finishCmList

    xoData
}


## ---------------- Process single family (for parallel execution) ------------

.processSingleFamily <- function(currentFamId, famSubset, genoMatrix, pedMap,
                                 geneticMap, parent, snpFilter, familySize) {
    currentFamily <- famSubset[famSubset[, 1] == currentFamId, ]

    ## Identify parents
    motherId <- .firstNonNA(currentFamily$mother[!is.na(currentFamily$mother)])
    fatherId <- .firstNonNA(currentFamily$father[!is.na(currentFamily$father)])

    ## Extract genotypes for this family
    familyGeno <- genoMatrix[, which(colnames(genoMatrix) %in% currentFamily[, 2]),
        drop = FALSE
    ]

    minCols <- if (familySize == 2L) 4L else 5L
    if (is.vector(familyGeno) || ncol(familyGeno) < minCols) {
        return(NULL)
    }

    motherColIdx <- which(colnames(familyGeno) == motherId)
    fatherColIdx <- which(colnames(familyGeno) == fatherId)

    if (length(motherColIdx) == 0L || length(fatherColIdx) == 0L) {
        return(NULL)
    }

    ## Convert to numeric and reorder (father, mother, children)
    genoRowNames <- rownames(familyGeno)
    familyGenoNumeric <- apply(familyGeno, 2, as.numeric)
    rownames(familyGenoNumeric) <- genoRowNames

    childColIndices <- setdiff(
        seq_len(ncol(familyGenoNumeric)),
        c(fatherColIdx, motherColIdx)
    )
    familyGenoOrdered <- familyGenoNumeric[, c(fatherColIdx, motherColIdx, childColIndices),
        drop = FALSE
    ]

    ## Get informative SNPs
    informativeGeno <- .getInformativeSNPs(familyGenoOrdered)
    if (is.null(informativeGeno)) {
        return(NULL)
    }

    ## Select parent-specific informative SNPs
    parentColIdx <- if (parent == "father") 1L else 2L
    parentSpecificGeno <- informativeGeno[informativeGeno[, parentColIdx] == 2, , drop = FALSE]

    ## Create binary inheritance matrix
    binaryInheritance <- .makeBinaryInheritance(parentSpecificGeno)
    if (is.null(binaryInheritance)) {
        return(NULL)
    }

    if (familySize == 2L) {
        return(.processTwoChildFamily(currentFamId, binaryInheritance, snpFilter))
    } else {
        currentFamilyOrdered <- currentFamily[c(fatherColIdx, motherColIdx, childColIndices), ]
        return(.processThreeChildFamily(
            currentFamId, binaryInheritance,
            currentFamilyOrdered, snpFilter
        ))
    }
}


## ---------------- Two-child family processing ----------------

.processTwoChildFamily <- function(currentFamId, binaryInheritance, snpFilter) {
    ## Assign inheritance states (1 = same, 2 = different)
    inheritanceState <- numeric(nrow(binaryInheritance))
    names(inheritanceState) <- rownames(binaryInheritance)

    for (i in seq_len(nrow(binaryInheritance))) {
        if ((binaryInheritance[i, 1] == 0 && binaryInheritance[i, 2] == 0) ||
            (binaryInheritance[i, 1] == 1 && binaryInheritance[i, 2] == 1)) {
            inheritanceState[i] <- 1
        } else {
            inheritanceState[i] <- 2
        }
    }

    if (length(inheritanceState) < 2L) {
        return(NULL)
    }

    ## Calculate state changes
    stateChange <- c(0, diff(inheritanceState))
    names(stateChange) <- names(inheritanceState)

    ## Apply SNP filter
    snpsToDrop <- .applySNPFilter(stateChange, snpFilter)

    ## Detect crossovers
    xoEvents <- NULL
    for (changeIdx in seq(2, length(stateChange))) {
        if (stateChange[changeIdx] != 0 &&
            !(names(stateChange)[changeIdx - 1] %in% snpsToDrop)) {
            eventData <- c(
                familyId = as.character(currentFamId),
                startSnp = names(stateChange)[changeIdx - 1],
                finishSnp = names(stateChange)[changeIdx]
            )
            xoEvents <- rbind(xoEvents, eventData)
        }
    }

    xoEvents
}


## ---------------- Three-child family processing ----------------

.processThreeChildFamily <- function(currentFamId, binaryInheritance,
                                     currentFamilyOrdered, snpFilter) {
    numChildren <- ncol(binaryInheritance)
    numTriples <- ceiling(numChildren / 3)

    allSnpsToDrop <- character()
    familyXoEvents <- NULL

    for (tripleIdx in seq_len(numTriples)) {
        ## Define triple indices
        tripleStartIdx <- (tripleIdx - 1) * 3 + 1
        tripleEndIdx <- min(tripleIdx * 3, numChildren)
        tripleColIndices <- tripleStartIdx:tripleEndIdx

        if (length(tripleColIndices) < 3) {
            tripleColIndices <- (numChildren - 2):numChildren
        }

        tripleGeno <- binaryInheritance[, tripleColIndices, drop = FALSE]
        tripleChildIds <- currentFamilyOrdered$member[(tripleColIndices + 2)]

        ## Assign states based on inheritance pattern
        inheritanceState <- numeric(nrow(tripleGeno))
        names(inheritanceState) <- rownames(tripleGeno)

        for (i in seq_len(nrow(tripleGeno))) {
            pattern <- tripleGeno[i, ]

            if (length(pattern) == 3) {
                if (all(pattern == 0) || all(pattern == 1)) {
                    inheritanceState[i] <- 1
                } else if ((pattern[1] == 0 && pattern[2] == 0 &&
                    pattern[3] == 1) ||
                    (pattern[1] == 1 && pattern[2] == 1 &&
                        pattern[3] == 0)) {
                    inheritanceState[i] <- 3
                } else if ((pattern[1] == 0 && pattern[2] == 1 &&
                    pattern[3] == 1) ||
                    (pattern[1] == 1 && pattern[2] == 0 &&
                        pattern[3] == 0)) {
                    inheritanceState[i] <- 7
                } else if ((pattern[1] == 1 && pattern[2] == 0 &&
                    pattern[3] == 1) ||
                    (pattern[1] == 0 && pattern[2] == 1 &&
                        pattern[3] == 0)) {
                    inheritanceState[i] <- 14
                }
            }
        }

        if (length(inheritanceState) < 10L) {
            next
        }

        ## Calculate state changes
        stateChange <- c(0, diff(inheritanceState))
        names(stateChange) <- names(inheritanceState)

        ## Apply SNP filter
        snpsToDrop <- .applySNPFilter(stateChange, snpFilter)
        allSnpsToDrop <- c(allSnpsToDrop, snpsToDrop)

        ## Detect crossovers based on state transitions
        for (changeIdx in seq(2, length(stateChange))) {
            stateChangeValue <- stateChange[changeIdx]

            affectedChildId <- NULL
            if (abs(stateChangeValue) == 2) {
                affectedChildId <- tripleChildIds[3]
            } else if (abs(stateChangeValue) == 4) {
                affectedChildId <- tripleChildIds[2]
            } else if (abs(stateChangeValue) %in% c(6, 11)) {
                affectedChildId <- tripleChildIds[1]
            } else if (abs(stateChangeValue) %in% c(7, 13)) {
                if (abs(stateChangeValue) == 7) {
                    affectedChildId <- tripleChildIds[3]
                } else {
                    affectedChildId <- tripleChildIds[2]
                }
            }

            if (!is.null(affectedChildId)) {
                eventData <- c(
                    childId = as.character(affectedChildId),
                    familyId = as.character(currentFamId),
                    startSnp = names(stateChange)[changeIdx - 1],
                    finishSnp = names(stateChange)[changeIdx]
                )
                familyXoEvents <- rbind(familyXoEvents, eventData)
            }
        }
    }

    ## Remove events flagged by SNP filter
    if (!is.null(familyXoEvents) && length(allSnpsToDrop) > 0L) {
        keepIndices <- !(familyXoEvents[, "startSnp"] %in% allSnpsToDrop)
        if (any(keepIndices)) {
            familyXoEvents <- familyXoEvents[keepIndices, , drop = FALSE]
        } else {
            familyXoEvents <- NULL
        }
    }

    familyXoEvents
}


## ---------------- Two-child families (with parallel support) ----------------

.xoDetectTwoChild <- function(pedData, genoMatrix, famData, geneticMap, parent,
                              snpFilter, cmFilter, BPPARAM) {
    ## Identify families with exactly 4 members (2 parents + 2 children)
    famCounts <- table(famData[, 1])
    famKeep <- names(famCounts)[famCounts == 4]
    famSubset <- famData[famData[, 1] %in% famKeep, , drop = FALSE]
    famIds <- unique(famSubset[, 1])
    pedMap <- pedData$map

    if (length(famIds) == 0L) {
        return(DataFrame(
            familyId = character(),
            startSnp = character(),
            finishSnp = character(),
            startPos = numeric(),
            finishPos = numeric(),
            startCm = numeric(),
            finishCm = numeric()
        ))
    }

    ## Process families in parallel
    xoList <- BiocParallel::bplapply(famIds, function(currentFamId) {
        .processSingleFamily(currentFamId, famSubset, genoMatrix, pedMap,
            geneticMap, parent, snpFilter,
            familySize = 2L
        )
    }, BPPARAM = BPPARAM)

    ## Remove NULL results
    xoList <- xoList[!vapply(xoList, is.null, logical(1))]

    ## Combine results
    if (length(xoList) == 0L) {
        return(DataFrame(
            familyId = character(),
            startSnp = character(),
            finishSnp = character(),
            startPos = numeric(),
            finishPos = numeric(),
            startCm = numeric(),
            finishCm = numeric()
        ))
    }

    xoAllEvents <- do.call(rbind, xoList)
    xoAllEvents <- as.data.frame(xoAllEvents, stringsAsFactors = FALSE)

    ## Add positions
    xoAllEvents <- .addPositions(xoAllEvents, pedMap, geneticMap)

    ## Remove duplicates
    xoAllEvents <- dplyr::distinct(xoAllEvents)

    ## Apply cM filter
    xoAllEvents <- .applyCMFilter(xoAllEvents, cmFilter)

    ## Convert to DataFrame
    resultDf <- DataFrame(
        familyId = xoAllEvents$familyId,
        startSnp = xoAllEvents$startSnp,
        finishSnp = xoAllEvents$finishSnp,
        startPos = as.numeric(xoAllEvents$startPos),
        finishPos = as.numeric(xoAllEvents$finishPos),
        startCm = as.numeric(xoAllEvents$startCm),
        finishCm = as.numeric(xoAllEvents$finishCm)
    )

    resultDf
}


## ---------------- Three-child families (with parallel support) ----------------

.xoDetectThreeChild <- function(pedData, genoMatrix, famData, geneticMap, parent,
                                snpFilter, cmFilter, BPPARAM) {
    ## Identify families with 5+ members (2 parents + 3+ children)
    famCounts <- table(famData[, 1])
    famKeep <- names(famCounts)[famCounts >= 5]
    famSubset <- famData[famData[, 1] %in% famKeep, , drop = FALSE]
    famIds <- unique(famSubset[, 1])
    pedMap <- pedData$map

    if (length(famIds) == 0L) {
        return(DataFrame(
            childId = character(),
            familyId = character(),
            startSnp = character(),
            finishSnp = character(),
            startPos = numeric(),
            finishPos = numeric(),
            startCm = numeric(),
            finishCm = numeric()
        ))
    }

    ## Process families in parallel
    xoList <- BiocParallel::bplapply(famIds, function(currentFamId) {
        .processSingleFamily(currentFamId, famSubset, genoMatrix, pedMap,
            geneticMap, parent, snpFilter,
            familySize = 3L
        )
    }, BPPARAM = BPPARAM)

    ## Remove NULL results
    xoList <- xoList[!vapply(xoList, is.null, logical(1))]

    ## Combine results
    if (length(xoList) == 0L) {
        return(DataFrame(
            childId = character(),
            familyId = character(),
            startSnp = character(),
            finishSnp = character(),
            startPos = numeric(),
            finishPos = numeric(),
            startCm = numeric(),
            finishCm = numeric()
        ))
    }

    xoAllEvents <- do.call(rbind, xoList)
    xoAllEvents <- as.data.frame(xoAllEvents, stringsAsFactors = FALSE)

    ## Add positions
    xoAllEvents <- .addPositions(xoAllEvents, pedMap, geneticMap)

    ## Remove duplicates
    xoAllEvents <- dplyr::distinct(xoAllEvents)

    ## Apply cM filter
    xoAllEvents <- .applyCMFilter(xoAllEvents, cmFilter)

    ## Convert to DataFrame
    resultDf <- DataFrame(
        childId = xoAllEvents$childId,
        familyId = xoAllEvents$familyId,
        startSnp = xoAllEvents$startSnp,
        finishSnp = xoAllEvents$finishSnp,
        startPos = as.numeric(xoAllEvents$startPos),
        finishPos = as.numeric(xoAllEvents$finishPos),
        startCm = as.numeric(xoAllEvents$startCm),
        finishCm = as.numeric(xoAllEvents$finishCm)
    )

    resultDf
}
