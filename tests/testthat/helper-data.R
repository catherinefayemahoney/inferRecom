# Helper functions to load test data
getTestDataPath <- function() {
    system.file("extdata", package = "inferRecom")
}

loadSimCeuData <- function() {
    dataPath <- getTestDataPath()
    list(
        plinkFile = file.path(dataPath, "simCEU"),
        mapMale = file.path(dataPath, "male_chr4.txt"),
        mapFemale = file.path(dataPath, "female_chr4.txt")
    )
}

# Check if test data exists
checkTestData <- function() {
    testData <- loadSimCeuData()
    all(file.exists(
        paste0(testData$plinkFile, ".bed"),
        paste0(testData$plinkFile, ".bim"),
        paste0(testData$plinkFile, ".fam"),
        testData$mapMale,
        testData$mapFemale
    ))
}

# Check if specific map exists
checkMapFile <- function(mapType = c("male", "female")) {
    mapType <- match.arg(mapType)
    testData <- loadSimCeuData()
    mapFile <- switch(mapType,
        male = testData$mapMale,
        female = testData$mapFemale
    )
    file.exists(mapFile)
}
