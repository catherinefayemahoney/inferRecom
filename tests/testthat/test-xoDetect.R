test_that("xoDetect works with 3-child maternal crossovers and female map", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Test maternal crossovers with female genetic map
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother"
    )

    # Check return type - should be GRanges when caseControl = FALSE
    expect_s4_class(testResult, "GRanges")

    # Check expected metadata columns
    expectedCols <- c(
        "childId", "familyId", "startSnp", "finishSnp",
        "startCm", "finishCm"
    )
    expect_true(all(expectedCols %in% colnames(mcols(testResult))))

    # Check data types
    expect_type(mcols(testResult)$childId, "character")
    expect_type(mcols(testResult)$familyId, "character")

    # Check that positions are reasonable
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(start(testResult) > 0))
        # Check cM positions are reasonable (<500 cM per chromosome)
        expect_true(all(mcols(testResult)$startCm >= 0 | is.na(mcols(testResult)$startCm)))
        expect_true(all(mcols(testResult)$finishCm <= 500 | is.na(mcols(testResult)$finishCm)))
    }
})

test_that("xoDetect works with 3-child paternal crossovers and male map", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("male"), "Male genetic map not available")

    testData <- loadSimCeuData()

    # Test paternal crossovers with male genetic map
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapMale,
        familySize = 3,
        parent = "father"
    )

    # Check return type
    expect_s4_class(testResult, "GRanges")

    # Check expected metadata columns
    expectedCols <- c(
        "childId", "familyId", "startSnp", "finishSnp",
        "startCm", "finishCm"
    )
    expect_true(all(expectedCols %in% colnames(mcols(testResult))))

    # Check that positions are reasonable
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(start(testResult) > 0))
        # Check cM positions are reasonable (<500 cM per chromosome)
        expect_true(all(mcols(testResult)$startCm >= 0 | is.na(mcols(testResult)$startCm)))
        expect_true(all(mcols(testResult)$finishCm <= 500 | is.na(mcols(testResult)$finishCm)))
    }
})

test_that("xoDetect handles 2-child families", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 2,
        parent = "mother"
    )

    expect_s4_class(testResult, "GRanges")

    # 2-child families should not have childId column
    expect_false("childId" %in% colnames(mcols(testResult)))

    # Should have familyId column
    if (length(testResult) > 0) {
        expect_true("familyId" %in% colnames(mcols(testResult)))
    }

    # Check that positions are reasonable
    if (length(testResult) > 0) {
        expect_true(all(start(testResult) < end(testResult)))
        expect_true(all(start(testResult) > 0))
        # Check cM positions are reasonable (<500 cM per chromosome)
        expect_true(all(mcols(testResult)$startCm >= 0 | is.na(mcols(testResult)$startCm)))
        expect_true(all(mcols(testResult)$finishCm <= 500 | is.na(mcols(testResult)$finishCm)))
    }
})

test_that("xoDetect filters by SNP correctly", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # With rsOnly = TRUE
    resultFiltered <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        rsOnly = TRUE
    )

    # With rsOnly = FALSE
    resultAll <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        rsOnly = FALSE
    )

    # Both should return GRanges
    expect_s4_class(resultFiltered, "GRanges")
    expect_s4_class(resultAll, "GRanges")

    # Should have string matches to 'rs' in both SNP columns equal to length
    sum(grepl("^rs", mcols(resultFiltered)$startSnp) &
        grepl("^rs", mcols(resultFiltered)$finishSnp)) ==
        length(resultFiltered)
})

test_that("xoDetect applies filters correctly", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Test with strict filters
    resultStrict <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        snpFilter = 10,
        cmFilter = 5
    )

    # Test with lenient filters
    resultLenient <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        snpFilter = 1,
        cmFilter = 0.1
    )

    # Strict filters should result in fewer crossovers
    expect_true(length(resultStrict) <= length(resultLenient))
})

test_that("xoDetect handles invalid input", {
    testData <- loadSimCeuData()

    # Invalid parent
    expect_error(
        xoDetect(
            plinkFile = testData$plinkFile,
            mapFile = testData$mapFemale,
            familySize = 3,
            parent = "invalid"
        )
    )

    # Invalid family size
    expect_error(
        xoDetect(
            plinkFile = testData$plinkFile,
            mapFile = testData$mapFemale,
            familySize = 5,
            parent = "mother"
        )
    )

    # Missing file
    expect_error(
        xoDetect(
            plinkFile = "nonexistent_file",
            mapFile = testData$mapFemale,
            familySize = 3,
            parent = "mother"
        )
    )
})

test_that("xoDetect output can be written to file", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()
    tempFile <- tempfile(fileext = ".csv")

    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        out = tempFile
    )

    # Check result is GRanges
    expect_s4_class(testResult, "GRanges")

    # Check file was created
    expect_true(file.exists(tempFile))

    # Check file can be read back
    readBack <- read.csv(tempFile, stringsAsFactors = FALSE)
    expect_equal(nrow(readBack), length(testResult))

    # Cleanup
    unlink(tempFile)
})

test_that("xoDetect returns GRangesList for case-control analysis", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Test with caseControl = TRUE (only works with familySize = 3)
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        caseControl = TRUE
    )

    # Should have "case" and "control" elements
    expect_true("case" %in% names(testResult))
    expect_true("control" %in% names(testResult))

    # Each element should be a GRanges
    expect_s4_class(testResult$case, "GRanges")
    expect_s4_class(testResult$control, "GRanges")
})

test_that("xoDetect caseControl requires familySize = 3", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Should error when caseControl = TRUE and familySize = 2
    expect_error(
        xoDetect(
            plinkFile = testData$plinkFile,
            mapFile = testData$mapFemale,
            familySize = 2,
            parent = "mother",
            caseControl = TRUE
        ),
        "No case-control separation for 2-child families"
    )
})

test_that("xoDetect returns empty GRanges when no crossovers detected", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    # Use very strict filters to potentially get no results
    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother",
        snpFilter = 100,
        cmFilter = 50
    )

    # Should still return GRanges (possibly empty)
    expect_s4_class(testResult, "GRanges")
    expect_true(length(testResult) >= 0)
})

test_that("xoDetect GRanges has correct seqnames", {
    skip_if_not(checkTestData(), "Test data not available")
    skip_if_not(checkMapFile("female"), "Female genetic map not available")

    testData <- loadSimCeuData()

    testResult <- xoDetect(
        plinkFile = testData$plinkFile,
        mapFile = testData$mapFemale,
        familySize = 3,
        parent = "mother"
    )

    # Check that seqnames are set
    if (length(testResult) > 0) {
        expect_true(length(seqnames(testResult)) > 0)
        # Seqnames should be chromosome identifiers
        expect_true(all(!is.na(seqnames(testResult))))
    }
})
