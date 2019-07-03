extdataPath <- system.file("extdata", package = "ALTBED")
filesetPath <- paste0(extdataPath, "/example")

nrows <- 50
ncols <- 1000
length <- nrows * ncols

parseRaw <- function(path) {
    lines <- strsplit(readLines(path), " ")
    header <- lines[[1]]
    data <- matrix(data = unlist(lines[2:length(lines)]), nrow = nrows, ncol = ncols + 6, byrow = TRUE)
    pheno <- data[, 1:6]
    geno <- data[, 7:ncol(data)]
    suppressWarnings(mode(geno) <- "integer")
    rownames(geno) <- paste0(pheno[, 1], "_", pheno[, 2])
    colnames(geno) <- header[7:length(header)]
    return(geno)
}

expect_error(ALTBED(path = filesetPath, n = -1))
expect_error(ALTBED(path = filesetPath, p = -1))
expect_error(ALTBED(path = filesetPath, n = -1, p = -1))

raw <- parseRaw(paste0(filesetPath, ".raw"))
bed <- suppressMessages(ALTBED(path = filesetPath))

expect_equal(raw[0], bed[0])
for (k in 1:length) {
    expect_equal(raw[k], bed[k])
}
expect_equal(raw[length + 1], bed[length + 1])

expect_equal(raw[0, ], bed[0, ])
expect_equal(raw[, 0], bed[, 0])
expect_equal(raw[0, 0], bed[0, 0])
for (i in 1:nrows) {
    for (j in 1:ncols) {
        expect_equal(raw[i, j], bed[i, j])
    }
}
expect_error(raw[nrows + 1, ])
expect_error(bed[nrows + 1, ])
expect_error(raw[, ncols + 1])
expect_error(bed[, ncols + 1])
expect_error(raw[nrows + 1, ncols + 1])
expect_error(bed[nrows + 1, ncols + 1])
