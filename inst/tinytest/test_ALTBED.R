extdataPath <- system.file("extdata", package = "ALTBED")
filesetPath <- paste0(extdataPath, "/example")

parseRaw <- function(path) {
    lines <- strsplit(readLines(path), " ")
    header <- lines[[1]]
    data <- matrix(data = unlist(lines[2:length(lines)]), nrow = 50, ncol = 1006, byrow = TRUE)
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

expect_equal(raw[1], bed[1])
