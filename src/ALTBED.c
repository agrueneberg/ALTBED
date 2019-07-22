#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Altrep.h>

#include <errno.h>

#include "mapping.h"
#include "bed.h"

#define ALTBED_EPTR(x) R_altrep_data1(x)
#define ALTBED_STATE(x) R_altrep_data2(x)
#define ALTBED_STATE_PATH(x) CAR(x)
#define ALTBED_STATE_LENGTH(x) CADR(x)
#define ALTBED_STATE_NROWS(x) CADDR(x)
#define ALTBED_STATE_NCOLS(x) CADDDR(x)

static R_altrep_class_t altrep_class_altbed;

static void *ALTBED_ADDR(SEXP x) {
    SEXP eptr = ALTBED_EPTR(x);
    void *addr = R_ExternalPtrAddr(eptr);
    if (addr == NULL) {
        Rf_error("object has been unmapped");
    }
    return addr;
}

static int ALTBED_NROWS(SEXP x) {
    SEXP state = ALTBED_STATE(x);
    SEXP nrows = ALTBED_STATE_NROWS(state);
    return Rf_asInteger(nrows);
}

static int ALTBED_NCOLS(SEXP x) {
    SEXP state = ALTBED_STATE(x);
    SEXP ncols = ALTBED_STATE_NCOLS(state);
    return Rf_asInteger(ncols);
}

/**
 * Object Creation
 */
static void unmake_altbed(SEXP eptr, SEXP state) {
    struct mapped_region mapped_region;
    mapped_region.addr = R_ExternalPtrAddr(eptr);
    SEXP length = ALTBED_STATE_LENGTH(state);
    mapped_region.length = Rf_asReal(length);
    if (unmap_region(&mapped_region) == -1) {
        Rf_error("could not unmap region");
    }
    R_ClearExternalPtr(eptr);
}

static void altbed_finalizer(SEXP eptr) {
    if (!R_ExternalPtrAddr(eptr)) {
        return;
    }
    SEXP state = R_ExternalPtrProtected(eptr);
    unmake_altbed(eptr, state);
}

static SEXP make_altbed_state(SEXP s_path, size_t length, SEXP s_nrows, SEXP s_ncols) {
    SEXP s_length = PROTECT(Rf_ScalarReal(length));
    SEXP s_state = PROTECT(Rf_list4(s_path, s_length, s_nrows, s_ncols));
    UNPROTECT(2); // s_length, s_state
    return s_state;
}

static SEXP make_altbed_eptr(void *p, SEXP state) {
    SEXP eptr = PROTECT(R_MakeExternalPtr(p, R_NilValue, state)); // protect state for use by the finalizer
    R_RegisterCFinalizerEx(eptr, altbed_finalizer, TRUE);
    UNPROTECT(1); // eptr
    return eptr;
}

static SEXP make_altbed(SEXP s_path, SEXP s_nrows, SEXP s_ncols) {
    const char *path = CHAR(Rf_asChar(s_path));
    int nrows = Rf_asInteger(s_nrows);
    if (nrows < 0) {
        Rf_error("nrows has to be nonnegative");
    }
    int ncols = Rf_asInteger(s_ncols);
    if (ncols < 0) {
        Rf_error("ncols has to be nonnegative");
    }
 // Map file
    struct mapped_region mapped_region;
    if (map_region(path, &mapped_region) == -1) {
        switch(errno) {
            case 1:
                Rf_error("could not find file");
                break;
            case 2:
                Rf_error("could not open file");
                break;
            case 3:
                Rf_error("could not create file mapping");
                break;
            case 4:
                Rf_error("could not map file");
                break;
            case 5:
                Rf_error("could not close file mapping");
                break;
            case 6:
                Rf_error("could not close file");
                break;
            case 7:
                Rf_error("file is not a regular file");
                break;
        }
    }
 // Test if valid .bed file
    if (is_bed_file(mapped_region.addr) == -1) {
        switch(errno) {
            case 1:
                Rf_error("file is not a binary PED file");
                break;
            case 2:
                Rf_error("sample-major mode is not supported");
                break;
        }
        unmap_region(&mapped_region); // ignore errors
    }
 // Test if n and p are correct
    if (has_valid_dimensions(mapped_region.length, nrows, ncols) == -1) {
        Rf_error("n or p does not match the dimensions of the file");
        unmap_region(&mapped_region); // ignore errors
    }
 // Create ALTBED object
    SEXP state = PROTECT(make_altbed_state(s_path, mapped_region.length, s_nrows, s_ncols));
    SEXP eptr = PROTECT(make_altbed_eptr(mapped_region.addr, state));
    SEXP altrep = PROTECT(R_new_altrep(altrep_class_altbed, eptr, state));
    MARK_NOT_MUTABLE(altrep); // always duplicate
    SEXP dim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(dim)[0] = nrows;
    INTEGER(dim)[1] = ncols;
    Rf_setAttrib(altrep, R_DimSymbol, dim);
    UNPROTECT(4); // state, eptr, altrep, dim
    return altrep;
}


/**
 * ALTREP Methods
 */
static R_xlen_t altbed_Length(SEXP x) {
    return (R_xlen_t) ALTBED_NROWS(x) * ALTBED_NCOLS(x);
}

static SEXP altbed_Duplicate(SEXP x, Rboolean deep) {
 // Not sure if this is OK...
    return x;
}

static SEXP altbed_Serialized_state(SEXP x) {
    return ALTBED_STATE(x);
}

static SEXP altbed_Unserialize(SEXP class, SEXP state) {
    SEXP path = ALTBED_STATE_PATH(state);
    SEXP nrows = ALTBED_STATE_NROWS(state);
    SEXP ncols = ALTBED_STATE_NCOLS(state);
    return make_altbed(path, nrows, ncols);
}


/**
 * ALTVEC Methods
 */
static void *altbed_Dataptr(SEXP x, Rboolean writeable) {
    Rf_error("cannot access data pointer for this ALTVEC object: contents of a .bed file cannot be cast to int or double -- this is a current limitation of the ALTREP framework");
}

static const void *altbed_Dataptr_or_null(SEXP x) {
 // No pointer available
    return NULL;
}

static SEXP altbed_Extract_subset(SEXP x, SEXP i, SEXP call) {
    uint8_t *ptr = ALTBED_ADDR(x);
    int nrows = ALTBED_NROWS(x);
    int num_bytes_per_variant = compute_num_bytes_per_variant(nrows);
    R_xlen_t nx = Rf_xlength(x);
    R_xlen_t ni = Rf_xlength(i);
    SEXP result = PROTECT(Rf_allocVector(INTSXP, ni));
    int *presult = INTEGER(result);
    R_xlen_t ci, ii;
    if (TYPEOF(i) == INTSXP) {
        int const *pi = INTEGER_RO(i);
        for (ci = 0; ci < ni; ci++) {
            ii = pi[ci];
            if (0 < ii && ii <= nx) {
                ii--;
                presult[ci] = extract_genotype(ptr, nrows, num_bytes_per_variant, ii, NA_INTEGER);
            } else {
                presult[ci] = NA_INTEGER;
            }
        }
    } else {
        double const *pi = REAL_RO(i);
        for (ci = 0; ci < ni; ci++) {
            double di = pi[ci];
            ii = (R_xlen_t) (di - 1);
            if (R_FINITE(di) && 0 <= ii && ii < nx) {
                presult[ci] = extract_genotype(ptr, nrows, num_bytes_per_variant, ii, NA_INTEGER);
            } else {
                presult[ci] = NA_INTEGER;
            }
        }
    }
    UNPROTECT(1); // result
    return result;
}


/**
 * ALTINTEGER Methods
 */
static int altbed_Elt(SEXP x, R_xlen_t i) {
    uint8_t *ptr = ALTBED_ADDR(x);
    int nrows = ALTBED_NROWS(x);
    int num_bytes_per_variant = compute_num_bytes_per_variant(nrows);
    return extract_genotype(ptr, nrows, num_bytes_per_variant, i, NA_INTEGER);
}


/**
 * Class Creation
 */
static void init_altbed_class(DllInfo *dll) {
    R_altrep_class_t class = R_make_altinteger_class("ALTBED", "ALTBED", dll);
    altrep_class_altbed = class;
 // Override ALTREP methods
    R_set_altrep_Length_method(class, altbed_Length);
    R_set_altrep_Duplicate_method(class, altbed_Duplicate);
    R_set_altrep_Serialized_state_method(class, altbed_Serialized_state);
    R_set_altrep_Unserialize_method(class, altbed_Unserialize);
 // Override ALTVEC methods
    R_set_altvec_Dataptr_method(class, altbed_Dataptr);
    R_set_altvec_Dataptr_or_null_method(class, altbed_Dataptr_or_null);
    R_set_altvec_Extract_subset_method(class, altbed_Extract_subset);
 // Override ALTINTEGER methods
    R_set_altinteger_Elt_method(class, altbed_Elt);
}


/**
 * Routines
 */
SEXP C_ALTBED(SEXP path, SEXP nrows, SEXP ncols) {
    return make_altbed(path, nrows, ncols);
}


/**
 * Shared Library Init
 */
static const R_CallMethodDef CallEntries[] = {
    {"C_ALTBED", (DL_FUNC) &C_ALTBED, 3},
    {NULL, NULL, 0}
};

void R_init_ALTBED(DllInfo *dll) {
    init_altbed_class(dll);
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}
