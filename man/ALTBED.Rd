\name{ALTBED}
\alias{ALTBED}
\title{Create an ALTBED Object from a PLINK .bed File}
\description{
    This function constructs a new \code{ALTBED} object by mapping the
    specified \href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK
        .bed} file into memory.
}
\usage{
ALTBED(path, n = NULL, p = NULL, simple_names = FALSE)
}
\arguments{
    \item{path}{
        Path to the
        \href{https://www.cog-genomics.org/plink2/formats#bed}{.bed} file (with
        or without extension).
    }
    \item{n}{
        The number of samples. If \code{NULL} (the default), this number will
        be determined from the accompanying
        \href{https://www.cog-genomics.org/plink2/formats#fam}{.fam} file (of
        the same name as the
        \href{https://www.cog-genomics.org/plink2/formats#bed}{.bed} file).  If
        a positive integer, the .fam file is not read and \code{rownames} will
        be set to \code{NULL} and have to be provided manually.
    }
    \item{p}{
        The number of variants. If \code{NULL} (the default) the number of
        variants will be determined from the accompanying
        \href{https://www.cog-genomics.org/plink2/formats#bim}{.bim} file (of
        the same name as the
        \href{https://www.cog-genomics.org/plink2/formats#bed}{.bed} file).  If
        a positive integer, the .bim file is not read and \code{colnames} will
        be set to \code{NULL} and have to be provided manually.
    }
    \item{simple_names}{
        Whether to simplify the format of the dimension names.  If \code{FALSE}
        (the default), row names are concatenations of family IDs, \code{_},
        and within-family IDs, while column names are concatenations of variant
        names, \code{_}, and minor alleles. If \code{TRUE}, row names are
        within-family IDs only and column names are variant names only.
    }
}
\details{
    \href{https://www.cog-genomics.org/plink2/formats#bed}{.bed} files must be
    accompanied by \href{https://www.cog-genomics.org/plink2/formats#fam}{.fam}
    and \href{https://www.cog-genomics.org/plink2/formats#bim}{.bim} files:
    .fam files contain sample information, and .bim files contain variant
    information. If the name of the .bed file is \emph{plink}.bed then the
    names of the .fam and .bim files have to be \emph{plink}.fam and
    \emph{plink}.bim, respectively. The .fam and .bim files are used to extract
    the number and names of samples and variants.

    For very large .bed files, reading the .fam and .bim files can take a long
    time. If \code{n} and \code{p} are provided, these files are not read and
    \code{dimnames} have to be provided manually.

    Currently, only the variant-major mode of .bed files is supported.
    \href{https://www.cog-genomics.org/plink2/}{PLINK2} "dropped" support for
    the sample-major mode by automatically converting files in this format to
    the variant-major mode. Therefore, it is recommended to run files in
    sample-major mode through PLINK2 first.
}
\value{
    An \code{ALTBED} object.
}
\examples{
# Get the path to the example .bed file
path <- system.file("extdata", "example.bed",
                    package = "ALTBED")

# Create an ALTBED object the example .bed file
m1 <- ALTBED(path)

# Create an ALTBED object the example .bed file without loading
# the .fam and .bim files
m2 <- ALTBED(path, n = 50, p = 1000)
}
