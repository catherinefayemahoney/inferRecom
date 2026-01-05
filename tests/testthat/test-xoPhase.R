## ============================================================================
## Setup and Helper Functions
## ============================================================================

test_that("xoPhase helper functions work correctly", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Test .firstNonNA
    expect_equal(inferRecom:::.firstNonNA(c(NA, NA, 5, 10)), 5)
    expect_equal(inferRecom:::.firstNonNA(c(1, 2, 3)), 1)
    expect_true(is.na(inferRecom:::.firstNonNA(c(NA, NA, NA))))

    # Test .checkSameChromosome with empty GRanges
    emptyGr <- GRanges()
    expect_true(inferRecom:::.checkSameChromosome(emptyGr))

    # Test .checkSameChromosome with single chromosome
    singleChrGr <- GRanges(
        seqnames = rep("chr4", 3),
        ranges = IRanges(start = c(1000, 2000, 3000), width = 100)
    )
    expect_true(inferRecom:::.checkSameChromosome(singleChrGr))

    # Test .checkSameChromosome with multiple chromosomes
    multiChrGr <- GRanges(
        seqnames = c("chr4", "chr5", "chr4"),
        ranges = IRanges(start = c(1000, 2000, 3000), width = 100)
    )
    expect_false(inferRecom:::.checkSameChromosome(multiChrGr))
})

test_that(".granges2DataFrame converts GRanges correctly", {
    # Test with metadata
    gr <- GRanges(
        seqnames = c("chr4", "chr4"),
        ranges = IRanges(start = c(1000, 2000), end = c(1500, 2500)),
        childId = c("child1", "child2"),
        familyId = c("fam1", "fam1")
    )

    df <- inferRecom:::.granges2DataFrame(gr)

    expect_true(is.data.frame(df))
    expect_equal(nrow(df), 2)
    expect_true("childId" %in% colnames(df))
    expect_true("familyId" %in% colnames(df))
    expect_true("start" %in% colnames(df))
    expect_true("end" %in% colnames(df))
    expect_equal(df$childId, c("child1", "child2"))

    # Test with empty GRanges
    emptyGr <- GRanges()
    emptyDf <- suppressWarnings(inferRecom:::.granges2DataFrame(emptyGr))
    expect_true(is.data.frame(emptyDf))
    expect_equal(nrow(emptyDf), 0)
})

## ============================================================================
## Input Validation Tests
## ============================================================================

test_that("xoPhase validates input arguments", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Load crossover data
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Test invalid plinkFile
    expect_error(
        xoPhase(
            plinkFile = c("file1", "file2"),
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat
        ),
        "'plinkFile' must be a single character string"
    )

    expect_error(
        xoPhase(
            plinkFile = 123,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat
        ),
        "'plinkFile' must be a single character string"
    )

    # Test invalid xoDetect objects
    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = data.frame(),
            xoDetectMaternal = xoMat
        ),
        "'xoDetectPaternal' must be a GRanges object"
    )

    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = list()
        ),
        "'xoDetectMaternal' must be a GRanges object"
    )

    # Test invalid famIds
    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat,
            famIds = 123
        ),
        "'famIds' must be a character vector or NULL"
    )

    # Test invalid rsOnly
    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat,
            rsOnly = "yes"
        ),
        "'rsOnly' must be a single logical value"
    )
})

test_that("xoPhase rejects mixed chromosomes", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Create GRanges with multiple chromosomes
    multiChrGr <- GRanges(
        seqnames = c("chr4", "chr5"),
        ranges = IRanges(start = c(1000, 2000), width = 100),
        childId = c("child1", "child1"),
        familyId = c("fam1", "fam1"),
        startSnp = c("rs1", "rs2"),
        finishSnp = c("rs3", "rs4")
    )

    singleChrGr <- GRanges(
        seqnames = "chr4",
        ranges = IRanges(start = 1000, width = 100),
        childId = "child1",
        familyId = "fam1",
        startSnp = "rs1",
        finishSnp = "rs2"
    )

    expect_error(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = multiChrGr,
            xoDetectMaternal = singleChrGr
        ),
        "GRanges object can only contain ranges from a single chromosome"
    )
})

## ============================================================================
## Basic Functionality Tests
## ============================================================================

test_that("xoPhase returns list output format", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Load crossover data
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        outputFormat = "list"
    )

    # Check return type
    expect_type(testResult, "list")
    expect_true(length(testResult) > 0)
    expect_true(!is.null(names(testResult)))

    # Check first family
    firstFamily <- testResult[[1]]

    if (is(firstFamily, "DataFrame") || is.data.frame(firstFamily)) {
        # Check required columns
        expect_true("rsID" %in% colnames(firstFamily))
        expect_true("location" %in% colnames(firstFamily))

        # Check for child haplotype columns (ending in Pat/Mat)
        childCols <- grep("(Pat|Mat)$", colnames(firstFamily), value = TRUE)
        expect_true(length(childCols) > 0)

        # Check for parent haplotype columns (ending in _1/_2)
        parentCols <- grep("_[12]$", colnames(firstFamily), value = TRUE)
        expect_true(length(parentCols) > 0)
    }
})

test_that("xoPhase phases multiple families", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    # Should have multiple families
    expect_true(length(testResult) >= 1)

    # Each family should have unique ID
    expect_equal(length(names(testResult)), length(unique(names(testResult))))
})

test_that("xoPhase handles specific family IDs", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Get all families first
    allFamilies <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    if (length(allFamilies) > 1) {
        # Select specific family
        targetFamily <- names(allFamilies)[1]

        testResult <- xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat,
            famIds = targetFamily
        )

        expect_equal(length(testResult), 1)
        expect_true(targetFamily %in% names(testResult))
    }
})

test_that("xoPhase filters rs SNPs correctly", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # With rsOnly = TRUE
    testResultRs <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        rsOnly = TRUE
    )

    # With rsOnly = FALSE
    testResultAll <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        rsOnly = FALSE
    )

    # Both should return lists
    expect_type(testResultRs, "list")
    expect_type(testResultAll, "list")

    # Check that rs filtering works
    if (length(testResultRs) > 0) {
        firstFam <- testResultRs[[1]]
        if (is(firstFam, "DataFrame") || is.data.frame(firstFam)) {
            rsSnps <- grep("^rs", firstFam$rsID, value = TRUE)
            expect_true(length(rsSnps) > 0)
        }
    }
})

## ============================================================================
## Output Format Tests
## ============================================================================

test_that("xoPhase returns SummarizedExperiment format", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not_installed("SummarizedExperiment")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        outputFormat = "summarizedExperiment"
    )

    expect_s4_class(testResult, "SummarizedExperiment")

    # Check assays
    expect_true("phased" %in% names(assays(testResult)))

    # Check rowData
    rData <- rowData(testResult)
    expect_true("snpId" %in% colnames(rData))

    # Check colData
    cData <- colData(testResult)
    expect_true("sampleId" %in% colnames(cData))
    expect_true("sampleType" %in% colnames(cData))
    expect_true("familyId" %in% colnames(cData))
})

test_that("xoPhase writes VCF files", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Create temporary output directory
    tempDir <- tempdir()
    vcfPrefix <- file.path(tempDir, "test_phased")

    testResult <- suppressMessages(
        xoPhase(
            plinkFile = testData$plinkFile,
            xoDetectPaternal = xoPat,
            xoDetectMaternal = xoMat,
            outputFormat = "vcf",
            vcfOutput = vcfPrefix
        )
    )

    # Check that VCF files were created
    expect_type(testResult, "character")
    expect_true(length(testResult) > 0)

    # Check that files exist
    for (vcfFile in testResult) {
        expect_true(file.exists(vcfFile))

        # Check file has content
        vcfContent <- readLines(vcfFile, n = 10)
        expect_true(any(grepl("^##fileformat=VCF", vcfContent)))
        expect_true(any(grepl("^#CHROM", vcfContent)))
    }

    # Clean up
    unlink(testResult)
})

## ============================================================================
## Phasing Quality Tests
## ============================================================================

test_that("xoPhase produces valid haplotypes", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    # Check first valid family
    for (famData in testResult) {
        if (is(famData, "DataFrame") || is.data.frame(famData)) {
            # Check that locations are numeric and sorted
            locations <- as.numeric(famData$location)
            expect_true(all(!is.na(locations)))
            expect_true(all(diff(locations) >= 0))

            # Check that SNP IDs are present
            expect_true(all(nchar(as.character(famData$rsID)) > 0))

            # Check that haplotypes contain valid nucleotides or NA
            validNuc <- c("A", "C", "G", "T", NA)
            for (col in colnames(famData)) {
                if (col %in% c("rsID", "location")) next
                expect_true(all(famData[[col]] %in% validNuc))
            }

            break # Only check first valid family
        }
    }
})

test_that("xoPhase maintains parent-child relationships", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    for (famData in testResult) {
        if (is(famData, "DataFrame") || is.data.frame(famData)) {
            # Get child columns
            childPatCols <- grep("Pat$", colnames(famData), value = TRUE)
            childMatCols <- grep("Mat$", colnames(famData), value = TRUE)

            # Each child should have both Pat and Mat columns
            childIds <- unique(gsub("(Pat|Mat)$", "", c(childPatCols, childMatCols)))

            for (childId in childIds) {
                expect_true(paste0(childId, "Pat") %in% colnames(famData))
                expect_true(paste0(childId, "Mat") %in% colnames(famData))
            }

            # Get parent columns
            parentCols <- grep("_[12]$", colnames(famData), value = TRUE)
            expect_true(length(parentCols) >= 2) # At least one parent with 2 haplotypes

            break
        }
    }
})

test_that("xoPhase handles homozygous sites correctly", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    testResult <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    for (famData in testResult) {
        if (is(famData, "DataFrame") || is.data.frame(famData)) {
            # For parent haplotypes, check that some SNPs have identical alleles
            # (homozygous sites)
            parentCols <- grep("_[12]$", colnames(famData), value = TRUE)

            if (length(parentCols) >= 2) {
                # Get pairs of parent haplotypes
                parentIds <- unique(gsub("_[12]$", "", parentCols))

                for (parentId in parentIds) {
                    hap1Col <- paste0(parentId, "_1")
                    hap2Col <- paste0(parentId, "_2")

                    if (hap1Col %in% colnames(famData) && hap2Col %in% colnames(famData)) {
                        # Some sites should be homozygous (both haplotypes same)
                        homSites <- famData[[hap1Col]] == famData[[hap2Col]]
                        homSites <- homSites[!is.na(homSites)]

                        if (length(homSites) > 0) {
                            # Expect at least some homozygous sites
                            expect_true(any(homSites))
                        }
                    }
                }
            }

            break
        }
    }
})

## ============================================================================
## Edge Case Tests
## ============================================================================

test_that("xoPhase handles empty crossover data", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Create empty GRanges
    emptyGr <- GRanges()

    expect_warning(
        expect_warning(
            testResult <- xoPhase(
                plinkFile = testData$plinkFile,
                xoDetectPaternal = emptyGr,
                xoDetectMaternal = emptyGr
            ),
            "Empty GRanges object provided"
        ),
        "Empty GRanges object provided"
    )

    # Should still return a list
    expect_type(testResult, "list")
})

test_that("xoPhase handles missing PLINK files", {
    xoPat <- GRanges(
        seqnames = "chr4",
        ranges = IRanges(start = 1000, width = 100),
        childId = "child1",
        familyId = "fam1"
    )

    expect_error(
        expect_warning(
            xoPhase(
                plinkFile = "nonexistent_file",
                xoDetectPaternal = xoPat,
                xoDetectMaternal = xoPat
            ),
            "cannot open file 'nonexistent_file.bim': No such file or directory"
        )
    )
})

## ============================================================================
## Parallel Processing Tests
## ============================================================================

test_that("xoPhase works with parallel processing", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_on_os("windows") # MulticoreParam doesn't work on Windows

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Serial execution
    resultSerial <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        BPPARAM = SerialParam()
    )

    # Parallel execution
    resultParallel <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat,
        BPPARAM = MulticoreParam(workers = 2)
    )

    # Results should be equivalent
    expect_equal(length(resultSerial), length(resultParallel))
    expect_equal(names(resultSerial), names(resultParallel))
})

## ============================================================================
## Internal Helper Function Tests
## ============================================================================

test_that(".matchAlleleCode works correctly", {
    # Father heterozygous (1), mother homozygous (2)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 2, 2), 2L)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 2, 1), 1L)

    # Father heterozygous (1), mother homozygous (0)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 0, 0), 1L)
    expect_equal(inferRecom:::.matchAlleleCode("father", 1, 0, 1), 2L)

    # Mother heterozygous (1), father homozygous (2)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 2, 1, 2), 2L)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 2, 1, 1), 1L)

    # Mother heterozygous (1), father homozygous (0)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 0, 1, 0), 1L)
    expect_equal(inferRecom:::.matchAlleleCode("mother", 0, 1, 1), 2L)

    # Unexpected combination
    expect_true(is.na(inferRecom:::.matchAlleleCode("father", 0, 0, 0)))
})

test_that(".getOtherAllele retrieves complementary alleles", {
    # Create mock BIM data
    bim <- data.frame(
        chr = c(4, 4, 4),
        snp = c("rs1", "rs2", "rs3"),
        cm = c(0, 0, 0),
        pos = c(1000, 2000, 3000),
        a1 = c("A", "C", "G"),
        a2 = c("T", "G", "A"),
        stringsAsFactors = FALSE
    )

    # Test getting other allele
    expect_equal(inferRecom:::.getOtherAllele("rs1", "A", bim), "T")
    expect_equal(inferRecom:::.getOtherAllele("rs1", "T", bim), "A")
    expect_equal(inferRecom:::.getOtherAllele("rs2", "C", bim), "G")
    expect_equal(inferRecom:::.getOtherAllele("rs3", "G", bim), "A")

    # Test missing SNP
    expect_true(is.na(inferRecom:::.getOtherAllele("rs999", "A", bim)))
})

test_that(".alignSnpData aligns SNPs correctly", {
    # Create test data
    snpDf <- data.frame(
        rsID = c("rs1", "rs3", "rs5"),
        location = c(1000, 3000, 5000),
        value = c("A", "C", "G"),
        stringsAsFactors = FALSE
    )

    allSnps <- c("rs1", "rs2", "rs3", "rs4", "rs5")

    aligned <- inferRecom:::.alignSnpData(snpDf, allSnps)

    expect_equal(nrow(aligned), 5)
    expect_equal(as.character(aligned$rsID), allSnps)
    expect_equal(aligned$value[1], "A")
    expect_true(is.na(aligned$value[2])) # rs2 was missing
    expect_equal(aligned$value[3], "C")
    expect_true(is.na(aligned$value[4])) # rs4 was missing
    expect_equal(aligned$value[5], "G")
})

test_that(".findNeighboringSNPs finds correct neighbors", {
    # Create mock phased data
    phasedData <- data.frame(
        rsID = c("rs1", "rs2", "rs3", "rs4", "rs5"),
        location = c(1000, 2000, 3000, 4000, 5000),
        stringsAsFactors = FALSE
    )

    phasedPosMap <- setNames(1:5, as.character(phasedData$location))

    # Test finding neighbors
    neighbors <- inferRecom:::.findNeighboringSNPs(2500, phasedData, phasedPosMap)

    expect_equal(neighbors$left, 2) # rs2 at position 2000
    expect_equal(neighbors$right, 3) # rs3 at position 3000

    # Test at start
    neighbors <- inferRecom:::.findNeighboringSNPs(500, phasedData, phasedPosMap)
    expect_null(neighbors$left)
    expect_equal(neighbors$right, 1)

    # Test at end
    neighbors <- inferRecom:::.findNeighboringSNPs(6000, phasedData, phasedPosMap)
    expect_equal(neighbors$left, 5)
    expect_null(neighbors$right)
})

## ============================================================================
## Integration Tests
## ============================================================================

test_that("xoPhase integrates with xoDetect output", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Run xoDetect to generate crossovers
    xoPat <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapMale,
        familySize = 3,
        parent = "father"
    )

    xoMat <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother"
    )

    # Phase using detected crossovers
    phased <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    # Verify output
    expect_type(phased, "list")
    expect_true(length(phased) > 0)
})

test_that("xoPhase output structure is consistent", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()
    xoPat <- readRDS(file.path(dirname(testData$plinkFile), "xoPat.rds"))
    xoMat <- readRDS(file.path(dirname(testData$plinkFile), "xoMat.rds"))

    # Run multiple times
    result1 <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    result2 <- xoPhase(
        plinkFile = testData$plinkFile,
        xoDetectPaternal = xoPat,
        xoDetectMaternal = xoMat
    )

    # Results should be identical
    expect_equal(length(result1), length(result2))
    expect_equal(names(result1), names(result2))

    # Check first family data is identical
    if (length(result1) > 0) {
        for (i in seq_along(result1)) {
            if (is(result1[[i]], "DataFrame") || is.data.frame(result1[[i]])) {
                expect_equal(nrow(result1[[i]]), nrow(result2[[i]]))
                expect_equal(ncol(result1[[i]]), ncol(result2[[i]]))
                expect_equal(colnames(result1[[i]]), colnames(result2[[i]]))
            }
        }
    }
})
