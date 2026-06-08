#ifndef ANEWDSC_H
#define ANEWDSC_H

#include <mpfr.h>
#include <mpfi.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------
 * Error Handling (Return Codes)
 * ------------------------------------------------------------------------- */
typedef enum {
    ANEWDSC_SUCCESS = 0,
    ANEWDSC_ERR_NOMEM = -1,
    ANEWDSC_ERR_INVALID_ARG = -2,
    ANEWDSC_ERR_MPFI = -3
} anewdsc_status_t;

/* -------------------------------------------------------------------------
 * Data Structures (Opaque structs with GMP-style init/clear)
 * ------------------------------------------------------------------------- */

// Represents a polynomial
typedef struct anewdsc_poly_struct {
    mpfi_t *coeffs;  // Array of interval coefficients
    size_t degree;
    size_t alloc;
} anewdsc_poly_struct;

typedef anewdsc_poly_struct anewdsc_poly_t[1];

// Represents a real interval (wrapper around mpfi_t)
typedef struct anewdsc_interval_struct {
    mpfi_t bounds;
} anewdsc_interval_struct;

typedef anewdsc_interval_struct anewdsc_interval_t[1];

/* -------------------------------------------------------------------------
 * Initialization and Memory Management
 * ------------------------------------------------------------------------- */
anewdsc_status_t anewdsc_poly_init(anewdsc_poly_t poly);
anewdsc_status_t anewdsc_poly_clear(anewdsc_poly_t poly);

anewdsc_status_t anewdsc_interval_init(anewdsc_interval_t interval);
anewdsc_status_t anewdsc_interval_clear(anewdsc_interval_t interval);

/* -------------------------------------------------------------------------
 * Interactive Builder
 * ------------------------------------------------------------------------- */
// Sets the coefficient for x^degree. The polynomial grows automatically.
anewdsc_status_t anewdsc_poly_set_coeff_str(anewdsc_poly_t poly, size_t degree, const char* str_val, int base);
anewdsc_status_t anewdsc_poly_set_coeff_d(anewdsc_poly_t poly, size_t degree, double val);
anewdsc_status_t anewdsc_poly_set_coeff_mpfr(anewdsc_poly_t poly, size_t degree, const mpfr_t val);

/* -------------------------------------------------------------------------
 * Polynomial Evaluation
 * ------------------------------------------------------------------------- */
// Evaluates the polynomial at the given interval x, storing the result in res.
// Uses Horner's method with interval arithmetic.
anewdsc_status_t anewdsc_poly_eval_interval(const anewdsc_poly_t poly, const anewdsc_interval_t x, anewdsc_interval_t res);

// Copies a polynomial.
anewdsc_status_t anewdsc_poly_copy(anewdsc_poly_t dest, const anewdsc_poly_t src);

// Performs a Taylor shift: res(x) = in(x + a)
anewdsc_status_t anewdsc_poly_taylor_shift(anewdsc_poly_t res, const anewdsc_poly_t in, const anewdsc_interval_t a);

// Scales the roots: res(x) = in(x * scale)
anewdsc_status_t anewdsc_poly_scale_roots(anewdsc_poly_t res, const anewdsc_poly_t in, const anewdsc_interval_t scale);

// Reverses the coefficients: res(x) = x^n * in(1/x)
anewdsc_status_t anewdsc_poly_reverse(anewdsc_poly_t res, const anewdsc_poly_t in);

// Computes the derivative of the polynomial: res(x) = in'(x)
anewdsc_status_t anewdsc_poly_derivative(anewdsc_poly_t res, const anewdsc_poly_t in);

/* -------------------------------------------------------------------------
 * Polynomial Arithmetic
 * ------------------------------------------------------------------------- */

// Adds two polynomials: res = p1 + p2
anewdsc_status_t anewdsc_poly_add(anewdsc_poly_t res, const anewdsc_poly_t p1, const anewdsc_poly_t p2);

// Subtracts two polynomials: res = p1 - p2
anewdsc_status_t anewdsc_poly_sub(anewdsc_poly_t res, const anewdsc_poly_t p1, const anewdsc_poly_t p2);

// Multiplies two polynomials: res = p1 * p2
anewdsc_status_t anewdsc_poly_mul(anewdsc_poly_t res, const anewdsc_poly_t p1, const anewdsc_poly_t p2);

// Divides two polynomials: dividend = quotient * divisor + remainder
// Returns ANEWDSC_ERR_INVALID_ARG if divisor is zero or its leading coefficient contains zero.
anewdsc_status_t anewdsc_poly_divrem(
    anewdsc_poly_t quotient, 
    anewdsc_poly_t remainder, 
    const anewdsc_poly_t dividend, 
    const anewdsc_poly_t divisor
);

/* -------------------------------------------------------------------------
 * Core Algorithm
 * ------------------------------------------------------------------------- */
typedef enum {
    ANEWDSC_VAR_ZERO = 0,
    ANEWDSC_VAR_ONE = 1,
    ANEWDSC_VAR_UNKNOWN = 2
} anewdsc_var_result_t;

// Performs the 01-Test according to Descartes.
anewdsc_status_t anewdsc_01_test(
    const anewdsc_poly_t poly, 
    const anewdsc_interval_t I, 
    anewdsc_var_result_t *result
);

// Finds a pseudo-admissible point near the midpoint of interval I.
// Stores the point as an exact interval [m_star, m_star] in m_star.
anewdsc_status_t anewdsc_find_admissible_point(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t I,
    anewdsc_interval_t m_star,
    int *found
);

// Performs the Newton-Test to accelerate convergence against root clusters.
// N_I is the convergence speed parameter.
// On success, sets success=1 and modifies I_prime.
anewdsc_status_t anewdsc_newton_test(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t I,
    size_t N_I,
    anewdsc_interval_t I_prime,
    int *success
);
// Isolates the real roots.
// output_roots: Pointer to an array of anewdsc_interval_t. Allocated internally or provided by user.
// num_roots: Number of found roots.
anewdsc_status_t anewdsc_isolate_roots(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t search_domain,
    anewdsc_interval_t **output_roots,
    size_t *num_roots
);

// Finds the first (smallest) real root in the given search domain.
// Returns ANEWDSC_SUCCESS and sets *found = 1 if a root is found, populating first_root.
// If no root exists, sets *found = 0.
anewdsc_status_t anewdsc_first_root_in_range(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t search_domain,
    anewdsc_interval_t first_root,
    int *found
);

// Helper function to free the output array
anewdsc_status_t anewdsc_free_roots_array(anewdsc_interval_t *output_roots, size_t num_roots);

#ifdef __cplusplus
}
#endif

#endif // ANEWDSC_H
