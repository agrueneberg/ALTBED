#ifndef BED_H
#define BED_H

#include <stddef.h>
#include <stdint.h>

int is_bed_file(uint8_t *bed);

int has_valid_dimensions(size_t length, int nrows, int ncols);

int extract_genotype(uint8_t *bed, int n, size_t k, int na_value);

#endif
