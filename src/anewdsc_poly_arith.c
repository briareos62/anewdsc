#include "anewdsc.h"
#include <stdlib.h>

// Forward declaration of internal function defined in anewdsc_poly.c
extern anewdsc_status_t ensure_capacity(anewdsc_poly_t poly, size_t required_degree);

anewdsc_status_t anewdsc_poly_add(anewdsc_poly_t res, const anewdsc_poly_t p1, const anewdsc_poly_t p2) {
    if (!res || !p1 || !p2) return ANEWDSC_ERR_INVALID_ARG;
    anewdsc_poly_t tmp;
    anewdsc_poly_init(tmp);
    size_t max_deg = p1->degree > p2->degree ? p1->degree : p2->degree;
    anewdsc_status_t status = ensure_capacity(tmp, max_deg);
    if (status != ANEWDSC_SUCCESS) { anewdsc_poly_clear(tmp); return status; }
    
    tmp->degree = max_deg;
    for (size_t i = 0; i <= max_deg; i++) {
        if (i <= p1->degree && i <= p2->degree) {
            mpfi_add(tmp->coeffs[i], p1->coeffs[i], p2->coeffs[i]);
        } else if (i <= p1->degree) {
            mpfi_set(tmp->coeffs[i], p1->coeffs[i]);
        } else {
            mpfi_set(tmp->coeffs[i], p2->coeffs[i]);
        }
    }
    
    status = anewdsc_poly_copy(res, tmp);
    anewdsc_poly_clear(tmp);
    return status;
}

anewdsc_status_t anewdsc_poly_sub(anewdsc_poly_t res, const anewdsc_poly_t p1, const anewdsc_poly_t p2) {
    if (!res || !p1 || !p2) return ANEWDSC_ERR_INVALID_ARG;
    anewdsc_poly_t tmp;
    anewdsc_poly_init(tmp);
    size_t max_deg = p1->degree > p2->degree ? p1->degree : p2->degree;
    anewdsc_status_t status = ensure_capacity(tmp, max_deg);
    if (status != ANEWDSC_SUCCESS) { anewdsc_poly_clear(tmp); return status; }
    
    tmp->degree = max_deg;
    for (size_t i = 0; i <= max_deg; i++) {
        if (i <= p1->degree && i <= p2->degree) {
            mpfi_sub(tmp->coeffs[i], p1->coeffs[i], p2->coeffs[i]);
        } else if (i <= p1->degree) {
            mpfi_set(tmp->coeffs[i], p1->coeffs[i]);
        } else {
            mpfi_neg(tmp->coeffs[i], p2->coeffs[i]);
        }
    }
    
    status = anewdsc_poly_copy(res, tmp);
    anewdsc_poly_clear(tmp);
    return status;
}

anewdsc_status_t anewdsc_poly_mul(anewdsc_poly_t res, const anewdsc_poly_t p1, const anewdsc_poly_t p2) {
    if (!res || !p1 || !p2) return ANEWDSC_ERR_INVALID_ARG;
    anewdsc_poly_t tmp;
    anewdsc_poly_init(tmp);
    
    size_t res_deg = p1->degree + p2->degree;
    anewdsc_status_t status = ensure_capacity(tmp, res_deg);
    if (status != ANEWDSC_SUCCESS) { anewdsc_poly_clear(tmp); return status; }
    
    tmp->degree = res_deg;
    for (size_t i = 0; i <= res_deg; i++) {
        mpfi_set_ui(tmp->coeffs[i], 0);
    }
    
    mpfi_t m_tmp;
    mpfi_init(m_tmp);
    
    for (size_t i = 0; i <= p1->degree; i++) {
        for (size_t j = 0; j <= p2->degree; j++) {
            mpfi_mul(m_tmp, p1->coeffs[i], p2->coeffs[j]);
            mpfi_add(tmp->coeffs[i+j], tmp->coeffs[i+j], m_tmp);
        }
    }
    mpfi_clear(m_tmp);
    
    status = anewdsc_poly_copy(res, tmp);
    anewdsc_poly_clear(tmp);
    return status;
}

anewdsc_status_t anewdsc_poly_divrem(
    anewdsc_poly_t quotient, 
    anewdsc_poly_t remainder, 
    const anewdsc_poly_t dividend, 
    const anewdsc_poly_t divisor
) {
    if (!quotient || !remainder || !dividend || !divisor) return ANEWDSC_ERR_INVALID_ARG;
    
    if (divisor->degree == 0 && mpfi_has_zero(divisor->coeffs[0])) {
        return ANEWDSC_ERR_INVALID_ARG;
    }
    if (mpfi_has_zero(divisor->coeffs[divisor->degree])) {
        return ANEWDSC_ERR_INVALID_ARG; // Leading coefficient cannot contain zero
    }
    
    anewdsc_poly_t q_tmp, r_tmp;
    anewdsc_poly_init(q_tmp);
    anewdsc_poly_init(r_tmp);
    
    anewdsc_poly_copy(r_tmp, dividend);
    
    if (dividend->degree < divisor->degree) {
        anewdsc_poly_set_coeff_d(q_tmp, 0, 0.0);
        anewdsc_poly_copy(quotient, q_tmp);
        anewdsc_poly_copy(remainder, r_tmp);
        anewdsc_poly_clear(q_tmp);
        anewdsc_poly_clear(r_tmp);
        return ANEWDSC_SUCCESS;
    }
    
    size_t q_deg = dividend->degree - divisor->degree;
    ensure_capacity(q_tmp, q_deg);
    q_tmp->degree = q_deg;
    for (size_t i = 0; i <= q_deg; i++) {
        mpfi_set_ui(q_tmp->coeffs[i], 0);
    }
    
    mpfi_t factor;
    mpfi_init(factor);
    
    while (r_tmp->degree >= divisor->degree) {
        size_t deg_diff = r_tmp->degree - divisor->degree;
        
        mpfi_div(factor, r_tmp->coeffs[r_tmp->degree], divisor->coeffs[divisor->degree]);
        mpfi_add(q_tmp->coeffs[deg_diff], q_tmp->coeffs[deg_diff], factor);
        
        mpfi_t term;
        mpfi_init(term);
        for (size_t i = 0; i <= divisor->degree; i++) {
            mpfi_mul(term, factor, divisor->coeffs[i]);
            mpfi_sub(r_tmp->coeffs[i + deg_diff], r_tmp->coeffs[i + deg_diff], term);
        }
        mpfi_clear(term);
        
        mpfi_set_ui(r_tmp->coeffs[r_tmp->degree], 0);
        
        if (r_tmp->degree == 0) break;
        r_tmp->degree--;
    }
    mpfi_clear(factor);
    
    anewdsc_poly_copy(quotient, q_tmp);
    anewdsc_poly_copy(remainder, r_tmp);
    anewdsc_poly_clear(q_tmp);
    anewdsc_poly_clear(r_tmp);
    
    return ANEWDSC_SUCCESS;
}
