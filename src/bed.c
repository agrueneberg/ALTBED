#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "bed.h"

#define PLINK_BED_HEADER_LENGTH 3
#define PLINK_BED_GENOTYPES_PER_BYTE 4

int is_bed_file(uint8_t *bed) {
 // Check magic number
    if (!(bed[0] == 0x6c && bed[1] == 0x1b)) {
        errno = 1;
        return -1;
    }
 // Check mode: 00000001 indicates the default variant-major mode (i.e.
 // list all samples for first variant, all samples for second variant,
 // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
 // all variants for the first sample, list all variants for the second
 // sample, etc)
    if (bed[2] != 0x01) {
        errno = 2;
        return -1;
    }
    return 0;
}

int has_valid_dimensions(size_t length, int nrows, int ncols) {
    int retval = 0;
 // File is a sequence of V blocks of N/4 (rounded up) bytes each, where V
 // is the number of variants and N is the number of samples.
    int num_bytes_per_block = ceil((double) nrows / PLINK_BED_GENOTYPES_PER_BYTE);
 // Check if given dimensions match the file
    if (((size_t) ncols * num_bytes_per_block) != (length - PLINK_BED_HEADER_LENGTH)) {
        retval = -1;
    }
    return retval;
}

int extract_genotype(uint8_t *bed, int n, int k, int na_value) {
 // Convert index from one-dimensional to two-dimensional
    int i = k % n;
    int j = k / n;
 // Each byte contains 4 genotypes; adjust indices
    int i_bytes = i / PLINK_BED_GENOTYPES_PER_BYTE;
    int n_bytes = n / PLINK_BED_GENOTYPES_PER_BYTE + (n % PLINK_BED_GENOTYPES_PER_BYTE != 0); // fast ceil for int
    int i_genotypes = 2 * (i - i_bytes * PLINK_BED_GENOTYPES_PER_BYTE);
 // Get byte from mapping
    uint8_t genotypes = bed[PLINK_BED_HEADER_LENGTH + (j * n_bytes + i_bytes)];
 // Extract genotype from byte by shifting bit pair of interest to the
 // right, then masking with 11
    uint8_t genotype = genotypes >> i_genotypes & 3;
 // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
 // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
 // allele. In BED, the coding is different: homozygous minor allele is 0
 // (00) and homozygous major allele is 3 (11). Each byte is read backwards.
    int mapping = na_value; // missing
    if (genotype == 0) {
        mapping = 2; // homozygous AA
    } else if (genotype == 3) {
        mapping = 0; // homozygous BB
    } else if (genotype == 2) {
        mapping = 1; // heterozygous AB
    }
    return mapping;
}
