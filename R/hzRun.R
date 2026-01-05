#' Detect Runs of Homozygosity in PLINK Genotype Data
#'
#' This function identifies runs of homozygosity (ROH) in PLINK genotype data
#' based on minimum physical length and SNP count thresholds. It can optionally
#' filter for specific samples, apply a minimum genetic distance criterion, and
#' separate results by case/control status.
#'
#' @param plinkFile Character string; path to the PLINK prefix (without
#'     extensions).
#' @param mapFile Character string or NULL; path to a tab-delimited genetic map
#'     file containing columns: \code{chr}, \code{pos}, \code{cM}. If NULL, cM
#'     filtering is not available. Default is NULL.
#' @param minMb Numeric; minimum physical length of ROH in megabases.
#'     Default is 1.
#' @param minSnps Integer; minimum number of consecutive homozygous SNPs
#'     required for an ROH. Default is 5.
#' @param minCm Numeric or NULL; minimum genetic length of ROH in centiMorgans.
#'     If NULL, no genetic length filtering is applied. Requires mapFile to be
#'     specified. Default is NULL.
#' @param rsOnly Logical; whether to restrict analysis to SNPs with names
#'     beginning with "rs". Default is FALSE.
#' @param sampleIds Character vector or NULL; specific sample IDs to analyze. If
#'     NULL, all samples are analyzed. Default is NULL.
#' @param caseControl Logical; if TRUE, function returns a
#'     \link[GenomicRanges:GRanges-class]{GRanges} of case and control ROH.
#'     Requires affected status in FAM file. Default is FALSE.
#' @param BPPARAM A BiocParallelParam object specifying parallel execution.
#'     Default is \code{SerialParam()} for serial execution. Use
#'     \code{MulticoreParam(workers = n)} for parallel processing across n cores.
#'
#' @return A \link[GenomicRanges:GRanges-class]{GRanges} object containing detected ROH with metadata
#'     columns: \code{sampleId}, \code{startSnp}, \code{finishSnp}, \code{numSnps}.
#'     If genetic map is provided, also includes \code{startCm}, \code{finishCm}.
#'     If \code{caseControl = TRUE}, returns a \link[GenomicRanges:GRangesList-class]{GRanges} with
#'     two elements: \code{case} and \code{control}. Returns an empty
#'     \link[GenomicRanges:GRanges-class]{GRanges} if no ROH are detected.
#'
#' @details
#' Runs of homozygosity (ROH) are continuous genomic segments where an
#' individual is homozygous at all marker positions. ROH can indicate
#' autozygosity (inheritance of identical haplotypes from a common ancestor)
#' and are used to estimate inbreeding, identify disease-associated loci,
#' and study population history.
#'
#' This function identifies ROH by:
#' \enumerate{
#'   \item Converting genotypes to binary (1 = homozygous, 0 = heterozygous/missing)
#'   \item Identifying runs using run-length encoding
#'   \item Filtering runs by minimum SNP count and physical length
#'   \item Optionally filtering by minimum genetic length (cM)
#' }
#'
#' Genotypes from \pkg{snpStats} are coded as:
#' \itemize{
#'   \item 0 = homozygous reference (AA) - counted as homozygous
#'   \item 1 = heterozygous (AB) - breaks ROH
#'   \item 2 = homozygous alternate (BB) - counted as homozygous
#'   \item NA = missing - breaks ROH
#' }
#'
#' Filtering strategy:
#' \itemize{
#'   \item Physical length (minMb): Always applied
#'   \item SNP count (minSnps): Always applied
#'   \item Genetic length (minCm): Applied only if mapFile is provided and minCm is not NULL
#' }
#'
#' Parallel processing is performed by individual, improving efficiency for datasets
#' with many samples.
#'
#' @import methods
#' @importFrom snpStats read.plink
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomicRanges GRanges GRangesList makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom utils read.delim
#' @importFrom BiocParallel SerialParam MulticoreParam bplapply
#' @importFrom stats setNames
#'
#' @examples
#' # Serial execution (default)
#'
#' # Get path to example data
#' dataPath <- system.file("extdata", package = "inferRecom")
#' plinkFile <- file.path(dataPath, "simCEU")
#'
#' # Detect ROH
#' rohData <- hzRun(
#'     plinkFile = plinkFile,
#'     minMb = 1,
#'     minSnps = 50
#' )
#'
#' # View results
#' head(rohData)
#'
#' # Analyze specific individuals
#' sampleIds <- c("F3C1", "F4C3", "F5P2")
#' rohSubset <- hzRun(
#'     plinkFile = plinkFile,
#'     minMb = .5,
#'     minSnps = 50,
#'     sampleIds = sampleIds
#' )
#'
#' rohSubset
#'
#' # Separate by case/control status
#' rohCc <- hzRun(
#'     plinkFile = plinkFile,
#'     minMb = .5,
#'     minSnps = 50,
#'     caseControl = TRUE
#' )
#'
#' rohCases <- rohCc$case
#' rohControls <- rohCc$control
#'
#' rohCases
#' rohControls
#'
#' @export
hzRun <- function(plinkFile,
                  mapFile = NULL,
                  minMb = 1,
                  minSnps = 5,
                  minCm = NULL,
                  rsOnly = FALSE,
                  sampleIds = NULL,
                  caseControl = FALSE,
                  BPPARAM = SerialParam()) {
    ## Argument validation
    if (!is.character(plinkFile) || length(plinkFile) != 1L) {
        stop("'plinkFile' must be a single character string")
    }

    if (!is.null(mapFile) && (!is.character(mapFile) || length(mapFile) != 1L)) {
        stop("'mapFile' must be a single character string or NULL")
    }

    if (!is.numeric(minMb) || length(minMb) != 1L || minMb <= 0) {
        stop("'minMb' must be a positive numeric value")
    }

    if (!is.numeric(minSnps) || length(minSnps) != 1L || minSnps <= 0) {
        stop("'minSnps' must be a positive integer")
    }
    minSnps <- as.integer(minSnps)

    if (!is.null(minCm)) {
        if (!is.numeric(minCm) || length(minCm) != 1L || minCm <= 0) {
            stop("'minCm' must be a positive numeric value or NULL")
        }
        if (is.null(mapFile)) {
            stop("'mapFile' must be provided when 'minCm' is specified")
        }
    }

    if (!is.logical(rsOnly) || length(rsOnly) != 1L) {
        stop("'rsOnly' must be a single logical value")
    }

    if (!is.null(sampleIds) && !is.character(sampleIds)) {
        stop("'sampleIds' must be a character vector or NULL")
    }

    if (!is.logical(caseControl) || length(caseControl) != 1L) {
        stop("'caseControl' must be a single logical value")
    }

    ## Check required packages
    pkgs <- c("snpStats", "BiocParallel")
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

    ## Convert SnpMatrix to numeric matrix
    ## as(object, "numeric") gives 0, 1, 2 for genotype counts, NA for missing
    genoSnpMat <- pedData$genotypes
    if (inherits(genoSnpMat, "SnpMatrix")) {
        ## Convert to numeric: 0=AA, 1=AB, 2=BB, NA=missing
        genoNumeric <- methods::as(genoSnpMat, "numeric")
        ## Transpose so SNPs are rows, samples are columns
        genoMatrix <- t(genoNumeric)
    } else {
        genoMatrix <- t(as.matrix(genoSnpMat))
    }

    ## Ensure matrix is numeric
    if (!is.numeric(genoMatrix)) {
        storage.mode(genoMatrix) <- "numeric"
    }

    mapData <- pedData$map
    famData <- pedData$fam

    ## Read genetic map if provided
    geneticMap <- NULL
    if (!is.null(mapFile)) {
        if (!file.exists(mapFile)) {
            stop("Genetic map file not found: ", mapFile)
        }

        geneticMapDf <- utils::read.delim(mapFile, stringsAsFactors = FALSE)

        if (!all(c("chr", "pos", "cM") %in% colnames(geneticMapDf))) {
            stop("mapFile must have columns: chr, pos, cM")
        }

        geneticMap <- data.frame(
            chr = geneticMapDf$chr,
            pos = geneticMapDf$pos,
            cM = geneticMapDf$cM,
            stringsAsFactors = FALSE
        )
    }

    ## Filter for rs SNPs if requested
    if (isTRUE(rsOnly)) {
        rsIdx <- grep("^rs", rownames(genoMatrix), ignore.case = TRUE)
        if (length(rsIdx) == 0L) {
            stop("No SNPs with 'rs' prefix found")
        }
        genoMatrix <- genoMatrix[rsIdx, , drop = FALSE]
    }

    ## Filter for specific IDs if provided
    if (!is.null(sampleIds)) {
        sampleIdx <- which(colnames(genoMatrix) %in% sampleIds)
        if (length(sampleIdx) == 0L) {
            stop("None of the specified IDs found in genotype data")
        }
        genoMatrix <- genoMatrix[, sampleIdx, drop = FALSE]
    }

    ## Get physical positions for SNPs
    snpNames <- rownames(genoMatrix)
    locationVec <- .getSnpLocations(snpNames, mapData)

    ## Convert threshold to base pairs
    minBp <- minMb * 1e6

    ## Detect ROH for all samples (with parallel support)
    rohList <- .detectRohAllSamples(
        genoMatrix, locationVec, snpNames,
        minSnps, minBp, geneticMap, minCm,
        BPPARAM
    )

    ## Convert to DataFrame
    if (length(rohList) == 0L) {
        ## Create empty DataFrame with appropriate columns
        if (!is.null(geneticMap)) {
            resultDf <- DataFrame(
                sampleId = character(),
                startSnp = character(),
                finishSnp = character(),
                numSnps = integer(),
                startPos = numeric(),
                finishPos = numeric(),
                startCm = numeric(),
                finishCm = numeric(),
                chromosome = character()
            )
        } else {
            resultDf <- DataFrame(
                sampleId = character(),
                startSnp = character(),
                finishSnp = character(),
                numSnps = integer(),
                startPos = numeric(),
                finishPos = numeric(),
                chromosome = character()
            )
        }
    } else {
        rohMat <- do.call(rbind, rohList)

        if (!is.null(geneticMap)) {
            resultDf <- DataFrame(
                sampleId = rohMat[, 1],
                startSnp = rohMat[, 2],
                finishSnp = rohMat[, 3],
                numSnps = as.integer(rohMat[, 4]),
                startPos = as.numeric(rohMat[, 5]),
                finishPos = as.numeric(rohMat[, 6]),
                startCm = as.numeric(rohMat[, 7]),
                finishCm = as.numeric(rohMat[, 8]),
                chromosome = mapData$chromosome[1]
            )
        } else {
            resultDf <- DataFrame(
                sampleId = rohMat[, 1],
                startSnp = rohMat[, 2],
                finishSnp = rohMat[, 3],
                numSnps = as.integer(rohMat[, 4]),
                startPos = as.numeric(rohMat[, 5]),
                finishPos = as.numeric(rohMat[, 6]),
                chromosome = mapData$chromosome[1]
            )
        }
    }

    ## Convert to GRanges or GRangesList
    if (nrow(resultDf) > 0L) {
        if (caseControl == FALSE) {
            # Convert to GRanges
            grResult <- makeGRangesFromDataFrame(
                df = resultDf,
                seqnames.field = "chromosome",
                start.field = "startPos",
                end.field = "finishPos",
                keep.extra.columns = TRUE
            )
        } else {
            # Add case/control status
            affectedIds <- famData$member[famData$affected == 1]
            resultDf$caseStatus <- ifelse(resultDf$sampleId %in% affectedIds,
                "case",
                "control"
            )

            # Convert to GRanges initially
            grTemp <- makeGRangesFromDataFrame(
                df = resultDf,
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

#' Get Physical Positions for SNPs
#'
#' Map SNP names to physical positions
#'
#' @param snpNames Character vector of SNP names
#' @param mapData Data frame with SNP map information
#'
#' @return Numeric vector of positions
#' @keywords internal
.getSnpLocations <- function(snpNames, mapData) {
    ## Create lookup table for matching
    lookupTable <- setNames(mapData$position, mapData$snp.name)

    positions <- lookupTable[snpNames]

    positions[is.na(positions)] <- 0

    as.numeric(positions)
}

#' Classify Genotype as Homozygous or Heterozygous
#'
#' Converts numeric genotype codes to binary (1 = homozygous, 0 = heterozygous)
#'
#' @param genoVec Numeric vector of genotype codes (0=AA, 1=AB, 2=BB, NA=missing)
#'
#' @return Integer vector: 1 for homozygous (0 or 2), 0 for heterozygous (1) or missing (NA)
#' @keywords internal
.isHomozygous <- function(genoVec) {
    ## Homozygous: 0 (AA) or 2 (BB)
    ## Heterozygous: 1 (AB)
    ## Missing: NA (treated as 0 to break runs)
    resultVec <- as.integer(genoVec == 0 | genoVec == 2)
    resultVec[is.na(resultVec)] <- 0L
    resultVec
}

#' Detect ROH for a Single Sample
#'
#' Identifies runs of homozygosity in one individual
#'
#' @param genoVec Numeric vector of genotypes for one sample (0, 1, 2, or NA)
#' @param locationVec Physical positions of SNPs in base pairs
#' @param snpNames SNP identifiers
#' @param sampleId Sample identifier
#' @param minSnps Minimum number of SNPs
#' @param minBp Minimum physical length in base pairs
#' @param geneticMap Data frame with genetic map (chr, pos, cM), or NULL
#' @param minCm Minimum genetic length in cM, or NULL
#'
#' @return Matrix with ROH information, or NULL if no ROH detected
#' @keywords internal
.detectRohSample <- function(genoVec, locationVec, snpNames, sampleId,
                             minSnps, minBp, geneticMap = NULL, minCm = NULL) {
    ## Convert to binary (homozygous = 1, heterozygous/missing = 0)
    homozygousVec <- .isHomozygous(genoVec)

    ## Run-length encoding to find consecutive homozygous stretches
    homozygousRuns <- rle(homozygousVec)

    ## Calculate start/end positions of runs
    endIdx <- cumsum(homozygousRuns$lengths)
    startIdx <- c(1L, endIdx[-length(endIdx)] + 1L)

    ## Get physical positions
    startBp <- locationVec[startIdx]
    endBp <- locationVec[endIdx]

    ## Filter for homozygous runs meeting SNP threshold
    isHomozygousRun <- homozygousRuns$values == 1L
    meetsSnpThreshold <- homozygousRuns$lengths >= minSnps
    passesSnpFilter <- isHomozygousRun & meetsSnpThreshold

    if (!any(passesSnpFilter)) {
        return(NULL)
    }

    ## Apply SNP filter
    startFiltered <- startIdx[passesSnpFilter]
    endFiltered <- endIdx[passesSnpFilter]
    startBpFiltered <- startBp[passesSnpFilter]
    endBpFiltered <- endBp[passesSnpFilter]

    ## Filter by physical length
    runLengthBp <- endBpFiltered - startBpFiltered
    meetsLengthThreshold <- runLengthBp >= minBp

    if (!any(meetsLengthThreshold)) {
        return(NULL)
    }

    ## Apply physical length filter
    startFinal <- startFiltered[meetsLengthThreshold]
    endFinal <- endFiltered[meetsLengthThreshold]
    startBpFinal <- startBpFiltered[meetsLengthThreshold]
    endBpFinal <- endBpFiltered[meetsLengthThreshold]

    ## Get genetic positions if available
    if (!is.null(geneticMap)) {
        startCm <- .getGeneticPosition(startBpFinal, geneticMap)
        endCm <- .getGeneticPosition(endBpFinal, geneticMap)

        ## Apply genetic length filter if specified
        if (!is.null(minCm)) {
            runLengthCm <- endCm - startCm
            meetsCmThreshold <- !is.na(runLengthCm) & runLengthCm >= minCm

            if (!any(meetsCmThreshold)) {
                return(NULL)
            }

            startFinal <- startFinal[meetsCmThreshold]
            endFinal <- endFinal[meetsCmThreshold]
            startBpFinal <- startBpFinal[meetsCmThreshold]
            endBpFinal <- endBpFinal[meetsCmThreshold]
            startCm <- startCm[meetsCmThreshold]
            endCm <- endCm[meetsCmThreshold]
        }
    } else {
        startCm <- rep(NA_real_, length(startFinal))
        endCm <- rep(NA_real_, length(startFinal))
    }

    ## Calculate final number of ROH
    numRoh <- length(startFinal)

    if (numRoh == 0L) {
        return(NULL)
    }

    ## Calculate number of SNPs in each run
    numSnpsPerRun <- endFinal - startFinal + 1L

    ## Build result matrix
    if (!is.null(geneticMap)) {
        ## 8 columns with cM information
        resultMat <- matrix(character(numRoh * 8L), nrow = numRoh, ncol = 8L)
        resultMat[, 1] <- sampleId
        resultMat[, 2] <- snpNames[startFinal]
        resultMat[, 3] <- snpNames[endFinal]
        resultMat[, 4] <- as.character(numSnpsPerRun)
        resultMat[, 5] <- as.character(startBpFinal)
        resultMat[, 6] <- as.character(endBpFinal)
        resultMat[, 7] <- as.character(startCm)
        resultMat[, 8] <- as.character(endCm)
    } else {
        ## 6 columns without cM information
        resultMat <- matrix(character(numRoh * 6L), nrow = numRoh, ncol = 6L)
        resultMat[, 1] <- sampleId
        resultMat[, 2] <- snpNames[startFinal]
        resultMat[, 3] <- snpNames[endFinal]
        resultMat[, 4] <- as.character(numSnpsPerRun)
        resultMat[, 5] <- as.character(startBpFinal)
        resultMat[, 6] <- as.character(endBpFinal)
    }

    resultMat
}

#' Process Single Individual for ROH Detection
#'
#' Wrapper function for parallel processing of a single individual
#'
#' @param sampleIdx Integer index of sample in genotype matrix
#' @param genoMatrix Matrix of genotypes (SNPs x samples)
#' @param locationVec Physical positions of SNPs
#' @param snpNames SNP identifiers
#' @param sampleIds Sample identifiers
#' @param minSnps Minimum number of SNPs
#' @param minBp Minimum physical length
#' @param geneticMap Genetic map data frame or NULL
#' @param minCm Minimum genetic length or NULL
#'
#' @return Matrix with ROH information, or NULL
#' @keywords internal
.processSingleIndividual <- function(sampleIdx, genoMatrix, locationVec,
                                     snpNames, sampleIds, minSnps, minBp,
                                     geneticMap, minCm) {
    .detectRohSample(
        genoVec = genoMatrix[, sampleIdx],
        locationVec = locationVec,
        snpNames = snpNames,
        sampleId = sampleIds[sampleIdx],
        minSnps = minSnps,
        minBp = minBp,
        geneticMap = geneticMap,
        minCm = minCm
    )
}

#' Detect ROH Across All Samples
#'
#' Efficiently processes all samples to detect runs of homozygosity
#' with optional parallel processing
#'
#' @param genoMatrix Matrix of genotypes (SNPs x samples)
#' @param locationVec Physical positions of SNPs in base pairs
#' @param snpNames SNP identifiers
#' @param minSnps Minimum number of SNPs
#' @param minBp Minimum physical length in base pairs
#' @param geneticMap Data frame with genetic map, or NULL
#' @param minCm Minimum genetic length in cM, or NULL
#' @param BPPARAM BiocParallelParam object for parallel execution
#'
#' @return List of matrices, one per sample with detected ROH
#' @keywords internal
.detectRohAllSamples <- function(genoMatrix, locationVec, snpNames,
                                 minSnps, minBp, geneticMap = NULL,
                                 minCm = NULL, BPPARAM = SerialParam()) {
    numSamples <- ncol(genoMatrix)
    sampleIds <- colnames(genoMatrix)

    ## Process all samples in parallel
    rohList <- BiocParallel::bplapply(seq_len(numSamples), function(i) {
        .processSingleIndividual(
            sampleIdx = i,
            genoMatrix = genoMatrix,
            locationVec = locationVec,
            snpNames = snpNames,
            sampleIds = sampleIds,
            minSnps = minSnps,
            minBp = minBp,
            geneticMap = geneticMap,
            minCm = minCm
        )
    }, BPPARAM = BPPARAM)

    ## Remove NULL results
    rohList <- rohList[!vapply(rohList, is.null, logical(1))]

    rohList
}

#' Get Genetic Position from Physical Position
#'
#' Maps physical position (bp) to genetic position (cM)
#'
#' @param posBp Numeric vector of physical positions in base pairs
#' @param geneticMap Data frame with chr, pos and cM columns
#'
#' @return Numeric vector of genetic positions in centiMorgans
#' @keywords internal
.getGeneticPosition <- function(posBp, geneticMap) {
    if (is.null(geneticMap) || nrow(geneticMap) == 0L) {
        return(rep(NA_real_, length(posBp)))
    }

    ## For each position, find nearest marker in genetic map
    cmPositions <- vapply(posBp, function(bp) {
        if (is.na(bp) || bp == 0) {
            return(NA_real_)
        }
        nearestIdx <- which.min(abs(geneticMap$pos - bp))
        if (length(nearestIdx) == 0L) {
            return(NA_real_)
        }
        geneticMap$cM[nearestIdx]
    }, numeric(1))

    cmPositions
}
