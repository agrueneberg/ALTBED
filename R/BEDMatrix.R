# Delimiters used in PED files
delims <- "[ \t]"

map <- function(path, n = NULL, p = NULL) {
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
                fam <- data.table::fread(famPath, select = c(1L, 2L), data.table = FALSE, showProgress = FALSE)
                # Determine n
                n <- nrow(fam)
                # Determine rownames
                rownames <- paste0(fam[, 1L], "_", fam[, 2L])
            } else {
                fam <- readLines(famPath)
                # Determine n
                n <- length(fam)
                # Determine rownames
                rownames <- sapply(strsplit(fam, delims), function(line) {
                    # Concatenate family ID and subject ID
                    return(paste0(line[1L], "_", line[2L]))
                })
            }
        }
    } else {
        n <- as.integer(n)
        # TODO: test if n is positive
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
                bim <- data.table::fread(bimPath, select = c(2L, 5L), data.table = FALSE, showProgress = FALSE)
                # Determine p
                p <- nrow(bim)
                # Determine colnames
                colnames <- paste0(bim[, 1L], "_", bim[, 2L])
            } else {
                bim <- readLines(bimPath)
                # Determine p
                p <- length(bim)
                # Determine colnames
                colnames <- sapply(strsplit(bim, delims), function(line) {
                    # Concatenate SNP name and minor allele (like --recodeA)
                    return(paste0(line[2L], "_", line[5L]))
                })
            }
        }
    } else {
        p <- as.integer(p)
        # TODO: test if p is positive
        colnames <- NULL
    }
    obj <- .Call(C_map, path, n, p)
    dimnames(obj) <- list(rownames, colnames)
    return(obj)
}

unmap <- function(x) {
    .Call(C_unmap, x)
}
