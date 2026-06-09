#include "anewdsc.h"
#include <stdlib.h>

anewdsc_status_t anewdsc_poly_init(anewdsc_poly_t poly) {
    if (!poly) return ANEWDSC_ERR_INVALID_ARG;
    poly->degree = 0;
    poly->alloc = 4;
    poly->coeffs = malloc(poly->alloc * sizeof(mpfi_t));
    if (!poly->coeffs) return ANEWDSC_ERR_NOMEM;
    
    for (size_t i = 0; i < poly->alloc; i++) {
        mpfi_init(poly->coeffs[i]);
        mpfi_set_ui(poly->coeffs[i], 0);
    }
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_clear(anewdsc_poly_t poly) {
    if (!poly) return ANEWDSC_ERR_INVALID_ARG;
    if (poly->coeffs) {
        for (size_t i = 0; i < poly->alloc; i++) {
            mpfi_clear(poly->coeffs[i]);
        }
        free(poly->coeffs);
        poly->coeffs = NULL;
    }
    poly->degree = 0;
    poly->alloc = 0;
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_interval_init(anewdsc_interval_t interval) {
    if (!interval) return ANEWDSC_ERR_INVALID_ARG;
    mpfi_init(interval->bounds);
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_interval_clear(anewdsc_interval_t interval) {
    if (!interval) return ANEWDSC_ERR_INVALID_ARG;
    mpfi_clear(interval->bounds);
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t ensure_capacity(anewdsc_poly_t poly, size_t required_degree) {
    if (required_degree < poly->alloc) return ANEWDSC_SUCCESS;
    
    size_t new_alloc = poly->alloc;
    while (new_alloc <= required_degree) {
        new_alloc *= 2;
    }
    
    mpfi_t *new_coeffs = realloc(poly->coeffs, new_alloc * sizeof(mpfi_t));
    if (!new_coeffs) return ANEWDSC_ERR_NOMEM;
    
    for (size_t i = poly->alloc; i < new_alloc; i++) {
        mpfi_init(new_coeffs[i]);
        mpfi_set_ui(new_coeffs[i], 0);
    }
    
    poly->coeffs = new_coeffs;
    poly->alloc = new_alloc;
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_set_coeff_str(anewdsc_poly_t poly, size_t degree, const char* str_val, int base) {
    if (!poly || !str_val) return ANEWDSC_ERR_INVALID_ARG;
    
    anewdsc_status_t status = ensure_capacity(poly, degree);
    if (status != ANEWDSC_SUCCESS) return status;
    
    if (mpfi_set_str(poly->coeffs[degree], str_val, base) != 0) {
        return ANEWDSC_ERR_MPFI;
    }
    
    if (degree > poly->degree) {
        poly->degree = degree;
    }
    
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_set_coeff_mpfr(anewdsc_poly_t poly, size_t degree, const mpfr_t val) {
    if (!poly) return ANEWDSC_ERR_INVALID_ARG;

    anewdsc_status_t status = ensure_capacity(poly, degree);
    if (status != ANEWDSC_SUCCESS) return status;

    mpfi_set_fr(poly->coeffs[degree], val);

    if (degree > poly->degree) {
        poly->degree = degree;
    }

    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_set_coeff_d(anewdsc_poly_t poly, size_t degree, double val) {
    if (!poly) return ANEWDSC_ERR_INVALID_ARG;
    
    anewdsc_status_t status = ensure_capacity(poly, degree);
    if (status != ANEWDSC_SUCCESS) return status;
    
    mpfi_set_d(poly->coeffs[degree], val);
    
    if (degree > poly->degree) {
        poly->degree = degree;
    }
    
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_eval_interval(const anewdsc_poly_t poly, const anewdsc_interval_t x, anewdsc_interval_t res) {
    if (!poly || !poly->coeffs || !x || !res) return ANEWDSC_ERR_INVALID_ARG;
    
    if (poly->degree == 0 && poly->alloc == 0) {
        mpfi_set_ui(res->bounds, 0);
        return ANEWDSC_SUCCESS;
    }
    
    // Horner's method: res = p_d + x * (p_{d-1} + x * (...))
    mpfi_set(res->bounds, poly->coeffs[poly->degree]);
    
    for (size_t i = poly->degree; i > 0; i--) {
        // res = res * x + p_{i-1}
        mpfi_mul(res->bounds, res->bounds, x->bounds);
        mpfi_add(res->bounds, res->bounds, poly->coeffs[i - 1]);
    }
    
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_copy(anewdsc_poly_t dest, const anewdsc_poly_t src) {
    if (!dest || !src) return ANEWDSC_ERR_INVALID_ARG;
    
    anewdsc_status_t status = ensure_capacity(dest, src->degree);
    if (status != ANEWDSC_SUCCESS) return status;
    
    dest->degree = src->degree;
    for (size_t i = 0; i <= src->degree; i++) {
        mpfi_set(dest->coeffs[i], src->coeffs[i]);
    }
    
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_taylor_shift(anewdsc_poly_t res, const anewdsc_poly_t in, const anewdsc_interval_t a) {
    if (!res || !in || !a) return ANEWDSC_ERR_INVALID_ARG;
    
    anewdsc_status_t status = anewdsc_poly_copy(res, in);
    if (status != ANEWDSC_SUCCESS) return status;
    
    if (res->degree == 0) return ANEWDSC_SUCCESS;
    
    mpfi_t tmp;
    mpfi_init(tmp);
    
    for (size_t i = 1; i <= res->degree; i++) {
        for (size_t k = res->degree; k >= i; k--) {
            size_t j = k - 1;
            mpfi_mul(tmp, a->bounds, res->coeffs[k]);
            mpfi_add(res->coeffs[j], res->coeffs[j], tmp);
        }
    }
    
    mpfi_clear(tmp);
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_scale_roots(anewdsc_poly_t res, const anewdsc_poly_t in, const anewdsc_interval_t scale) {
    if (!res || !in || !scale) return ANEWDSC_ERR_INVALID_ARG;
    anewdsc_status_t status = ensure_capacity(res, in->degree);
    if (status != ANEWDSC_SUCCESS) return status;
    
    res->degree = in->degree;
    
    mpfi_t current_scale;
    mpfi_init(current_scale);
    mpfi_set_ui(current_scale, 1);
    
    for (size_t i = 0; i <= in->degree; i++) {
        mpfi_mul(res->coeffs[i], in->coeffs[i], current_scale);
        if (i < in->degree) {
            mpfi_mul(current_scale, current_scale, scale->bounds);
        }
    }
    
    mpfi_clear(current_scale);
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_reverse(anewdsc_poly_t res, const anewdsc_poly_t in) {
    if (!res || !in) return ANEWDSC_ERR_INVALID_ARG;
    anewdsc_status_t status = ensure_capacity(res, in->degree);
    if (status != ANEWDSC_SUCCESS) return status;
    
    res->degree = in->degree;
    
    if (res == in) {
        size_t d = in->degree;
        mpfi_t tmp;
        mpfi_init(tmp);
        for (size_t i = 0; i < (d + 1) / 2; i++) {
            mpfi_set(tmp, res->coeffs[i]);
            mpfi_set(res->coeffs[i], res->coeffs[d - i]);
            mpfi_set(res->coeffs[d - i], tmp);
        }
        mpfi_clear(tmp);
    } else {
        for (size_t i = 0; i <= in->degree; i++) {
            mpfi_set(res->coeffs[i], in->coeffs[in->degree - i]);
        }
    }
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_poly_derivative(anewdsc_poly_t res, const anewdsc_poly_t in) {
    if (!res || !in) return ANEWDSC_ERR_INVALID_ARG;
    if (in->degree == 0) {
        anewdsc_status_t status = ensure_capacity(res, 0);
        if (status != ANEWDSC_SUCCESS) return status;
        res->degree = 0;
        mpfi_set_ui(res->coeffs[0], 0);
        return ANEWDSC_SUCCESS;
    }
    anewdsc_status_t status = ensure_capacity(res, in->degree - 1);
    if (status != ANEWDSC_SUCCESS) return status;
    
    res->degree = in->degree - 1;
    for (size_t i = 1; i <= in->degree; i++) {
        mpfi_mul_ui(res->coeffs[i-1], in->coeffs[i], i);
    }
    return ANEWDSC_SUCCESS;
}

