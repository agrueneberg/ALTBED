# Delimiters used in .fam and .bim files
delims <- "[ \t]"

map <- function(path, n = NULL, p = NULL, simple_names = FALSE) {
    path <- path.expand(path)
    if (!file.exists(path)) {
        # Try to add extension (common in PLINK)
        path <- paste0(path, ".bed")
        if (!file.exists(path)) {
            stop("File not found.")
        }
    }
    pathSansExt <- tools::file_path_sans_ext(path)
    filesetName <- basename(pathSansExt)
    if (is.null(n)) {
        # Check if FAM file exists
        famPath <- paste0(pathSansExt, ".fam")
        if (!file.exists(famPath)) {
            stop(filesetName, ".fam not found. Provide number of samples (n).")
        } else {
            message("Extracting number of samples and rownames from ", filesetName, ".fam...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                if (simple_names) {
                    famColumns <- c(2L)
                } else {
                    famColumns <- c(1L, 2L)
                }
                fam <- data.table::fread(famPath, select = famColumns, colClasses = list(character = famColumns), data.table = FALSE, showProgress = FALSE)
                # Determine n
                n <- nrow(fam)
                # Determine rownames
                if (simple_names) {
                    # Use within-family ID only
                    rownames <- fam[, 1L]
                } else {
                    # Concatenate family ID and within-family ID
                    rownames <- paste0(fam[, 1L], "_", fam[, 2L])
                }
            } else {
                fam <- readLines(famPath) # much faster than read.table
                # Determine n
                n <- length(fam)
                # Determine rownames
                if (simple_names) {
                    rownames <- sapply(strsplit(fam, delims), function(line) {
                        # Use within-family ID only
                        line[2L]
                    })
                } else {
                    rownames <- sapply(strsplit(fam, delims), function(line) {
                        # Concatenate family ID and within-family ID
                        paste0(line[1L], "_", line[2L])
                    })
                }
            }
        }
    } else {
        n <- as.integer(n)
        rownames <- NULL
    }
    if (is.null(p)) {
        # Check if BIM file exists
        bimPath <- paste0(pathSansExt, ".bim")
        if (!file.exists(bimPath)) {
            stop(filesetName, ".bim not found. Provide number of variants (p).")
        } else {
            message("Extracting number of variants and colnames from ", filesetName, ".bim...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                if (simple_names) {
                    bimColumns <- c(2L)
                } else {
                    bimColumns <- c(2L, 5L)
                }
                bim <- data.table::fread(bimPath, select = bimColumns, colClasses = list(character = bimColumns), data.table = FALSE, showProgress = FALSE)
                # Determine p
                p <- nrow(bim)
                # Determine colnames
                if (simple_names) {
                    # Use variant name only
                    colnames <- bim[, 1L]
                } else {
                    # Concatenate variant name and minor allele (like --recodeA)
                    colnames <- paste0(bim[, 1L], "_", bim[, 2L])
                }
            } else {
                bim <- readLines(bimPath) # much faster than read.table
                # Determine p
                p <- length(bim)
                # Determine colnames
                if (simple_names) {
                    colnames <- sapply(strsplit(bim, delims), function(line) {
                        # Use variant name only
                        line[2L]
                    })
                } else {
                    colnames <- sapply(strsplit(bim, delims), function(line) {
                        # Concatenate variant name and minor allele (like --recodeA)
                        paste0(line[2L], "_", line[5L])
                    })
                }
            }
        }
    } else {
        p <- as.integer(p)
        colnames <- NULL
    }
    obj <- .Call(C_map, path, n, p)
    dimnames(obj) <- list(rownames, colnames)
    return(obj)
}
