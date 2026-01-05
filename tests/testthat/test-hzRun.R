test_that("hzRun detects runs of homozygosity", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    testResult <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = 1,
        minSnps = 50
    )

    # Check return type - should be GRanges when caseControl = FALSE
    expect_s4_class(testResult, "GRanges")

    # Check expected metadata columns
    expectedCols <- c("sampleId", "startSnp", "finishSnp", "numSnps")
    expect_true(all(expectedCols %in% colnames(mcols(testResult))))

    # Check ROH conform to arguments
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(mcols(testResult)$numSnps >= 50))

        # Check length in Mb
        lengthBp <- end(testResult) - start(testResult)
        expect_true(all(lengthBp >= 1e6))
    }
})

test_that("hzRun filters by minimum length", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Short ROH
    resultShort <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = 0.5,
        minSnps = 10
    )

    # Long ROH
    resultLong <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = 1,
        minSnps = 100
    )

    # Should detect more short ROH than long ROH
    expect_true(length(resultShort) >= length(resultLong))
})

test_that("hzRun filters by minimum genetic distance", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Short genetic distance
    resultShort <- hzRun(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        minMb = 0.5,
        minSnps = 10,
        minCm = 0.3
    )

    # Long genetic distance
    resultLong <- hzRun(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        minMb = 1,
        minSnps = 100,
        minCm = 0.8
    )

    # Should detect more short ROH than long ROH
    expect_true(length(resultShort) >= length(resultLong))
})

test_that("hzRun handles case/control split", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Use caseControl = TRUE
    testResult <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = .5,
        minSnps = 10,
        caseControl = TRUE
    )

    # Should have "case" and "control" elements
    expect_true("case" %in% names(testResult))
    expect_true("control" %in% names(testResult))

    # Output list items are GRanges
    expect_s4_class(testResult$case, "GRanges")
    expect_s4_class(testResult$control, "GRanges")
})

test_that("hzRun filters for specific sample IDs", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Get all samples first
    allRoh <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = .1,
        minSnps = 10
    )

    if (length(allRoh) > 0) {
        # Select subset of samples
        uniqueSampleIds <- unique(mcols(allRoh)$sampleId)
        selectedSampleIds <- uniqueSampleIds[1:min(2, length(uniqueSampleIds))]

        testResult <- hzRun(
            plinkFile = testData$plinkFile,
            minMb = .5,
            minSnps = 50,
            sampleIds = selectedSampleIds
        )

        # Should only contain specified samples
        expect_true(all(mcols(testResult)$sampleId %in% selectedSampleIds))
    }
})

test_that("hzRun filters by rs SNPs only", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # With rsOnly = TRUE
    resultFiltered <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = 1,
        minSnps = 50,
        rsOnly = TRUE
    )

    # With rsOnly = FALSE
    resultAll <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = 1,
        minSnps = 50,
        rsOnly = FALSE
    )

    # Both should return GRanges
    expect_s4_class(resultFiltered, "GRanges")
    expect_s4_class(resultAll, "GRanges")
})

test_that("hzRun returns empty GRanges when no ROH detected", {
    skip_if_not(checkTestData(), "Test data not available")

    testData <- loadSimCeuData()

    # Use very strict filters to potentially get no results
    testResult <- hzRun(
        plinkFile = testData$plinkFile,
        minMb = 100,
        minSnps = 10000
    )

    # Should still return GRanges (possibly empty)
    expect_s4_class(testResult, "GRanges")
})

test_that("hzRun handles invalid input", {
    testData <- loadSimCeuData()

    # Invalid minMb
    expect_error(
        hzRun(
            plinkFile = testData$plinkFile,
            minMb = -1,
            minSnps = 50
        ),
        "'minMb' must be a positive numeric value"
    )

    # Invalid minSnps
    expect_error(
        hzRun(
            plinkFile = testData$plinkFile,
            minMb = 1,
            minSnps = -10
        ),
        "'minSnps' must be a positive integer"
    )

    # Missing file
    expect_error(
        suppressWarnings(
            hzRun(
                plinkFile = "nonexistent_file",
                minMb = 1,
                minSnps = 50
            )
        ),
        regexp = NULL
    )
})

test_that(".getGeneticPosition maps physical to genetic positions", {
    # Create mock genetic map
    geneticMap <- data.frame(
        chr = rep("chr4", 5),
        pos = c(1000, 5000, 10000, 15000, 20000),
        cM = c(0.1, 0.5, 1.0, 1.5, 2.0),
        stringsAsFactors = FALSE
    )

    # Test exact match
    result <- inferRecom:::.getGeneticPosition(c(1000, 10000), geneticMap)
    expect_equal(result, c(0.1, 1.0))

    # Test nearest neighbor (should find closest position)
    result <- inferRecom:::.getGeneticPosition(c(1100, 4900), geneticMap)
    expect_equal(result, c(0.1, 0.5))

    # Test position exactly between two markers
    result <- inferRecom:::.getGeneticPosition(7500, geneticMap)
    expect_true(result %in% c(0.5, 1.0)) # Could match either

    # Test NA input
    result <- inferRecom:::.getGeneticPosition(NA, geneticMap)
    expect_true(is.na(result))

    # Test zero position
    result <- inferRecom:::.getGeneticPosition(0, geneticMap)
    expect_true(is.na(result))

    # Test NULL genetic map
    result <- inferRecom:::.getGeneticPosition(c(1000, 5000), NULL)
    expect_equal(result, c(NA_real_, NA_real_))

    # Test empty genetic map
    emptyMap <- data.frame(chr = character(), pos = numeric(), cM = numeric())
    result <- inferRecom:::.getGeneticPosition(c(1000), emptyMap)
    expect_true(is.na(result))

    # Test multiple positions
    result <- inferRecom:::.getGeneticPosition(c(1000, 5000, 10000), geneticMap)
    expect_equal(length(result), 3)
    expect_equal(result, c(0.1, 0.5, 1.0))

    # Test position beyond range
    result <- inferRecom:::.getGeneticPosition(100000, geneticMap)
    expect_equal(result, 2.0) # Should match furthest marker
})

test_that(".getGeneticPosition handles edge cases", {
    geneticMap <- data.frame(
        chr = "chr4",
        pos = c(1000, 5000),
        cM = c(0.1, 0.5)
    )

    # Test vector input
    posBp <- c(1000, 3000, 5000, 7000)
    result <- inferRecom:::.getGeneticPosition(posBp, geneticMap)
    expect_equal(length(result), 4)
    expect_false(any(is.na(result[1:3])))

    # Test monotonicity (positions should generally increase with bp)
    largMap <- data.frame(
        chr = rep("chr4", 10),
        pos = seq(1000, 100000, length.out = 10),
        cM = seq(0.1, 10, length.out = 10)
    )
    testPos <- seq(1000, 100000, length.out = 5)
    result <- inferRecom:::.getGeneticPosition(testPos, largMap)
    expect_true(all(diff(result) >= 0))
})
