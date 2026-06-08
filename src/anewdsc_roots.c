#include "anewdsc.h"
#include <stdlib.h>

anewdsc_status_t anewdsc_01_test(
    const anewdsc_poly_t poly, 
    const anewdsc_interval_t I, 
    anewdsc_var_result_t *result
) {
    if (!poly || !I || !result) return ANEWDSC_ERR_INVALID_ARG;
    
    mpfr_t a_fr, b_fr;
    mpfr_init(a_fr); mpfr_init(b_fr);
    mpfi_get_left(a_fr, I->bounds);
    mpfi_get_right(b_fr, I->bounds);
    
    anewdsc_interval_t a_iv, b_minus_a_iv, one_iv;
    anewdsc_interval_init(a_iv);
    anewdsc_interval_init(b_minus_a_iv);
    anewdsc_interval_init(one_iv);
    
    mpfi_set_fr(a_iv->bounds, a_fr);
    
    mpfi_t tmp_b;
    mpfi_init(tmp_b);
    mpfi_set_fr(tmp_b, b_fr);
    mpfi_sub(b_minus_a_iv->bounds, tmp_b, a_iv->bounds);
    mpfi_clear(tmp_b);
    
    mpfi_set_ui(one_iv->bounds, 1);
    
    anewdsc_poly_t P1, P2, P3, P_I;
    anewdsc_poly_init(P1); anewdsc_poly_init(P2); anewdsc_poly_init(P3); anewdsc_poly_init(P_I);
    
    anewdsc_poly_taylor_shift(P1, poly, a_iv);
    anewdsc_poly_scale_roots(P2, P1, b_minus_a_iv);
    anewdsc_poly_reverse(P3, P2);
    anewdsc_poly_taylor_shift(P_I, P3, one_iv);
    
    long min_var_pos = -1, max_var_pos = -1;
    long min_var_neg = -1, max_var_neg = -1;
    int has_first = 0;
    
    for (size_t i = 0; i <= P_I->degree; i++) {
        mpfi_get_right(b_fr, P_I->coeffs[i]);
        int can_be_pos_val = (mpfr_cmp_d(b_fr, 0.0) > 0);
        
        mpfi_get_left(a_fr, P_I->coeffs[i]);
        int can_be_neg_val = (mpfr_cmp_d(a_fr, 0.0) < 0);
        
        if (!can_be_pos_val && !can_be_neg_val) {
            continue;
        }
        
        if (!has_first) {
            if (can_be_pos_val) { min_var_pos = 0; max_var_pos = 0; }
            if (can_be_neg_val) { min_var_neg = 0; max_var_neg = 0; }
            has_first = 1;
        } else {
            long new_min_pos = -1, new_max_pos = -1;
            long new_min_neg = -1, new_max_neg = -1;
            
            if (can_be_pos_val) {
                long v1_min = (min_var_pos != -1) ? min_var_pos : -1;
                long v2_min = (min_var_neg != -1) ? min_var_neg + 1 : -1;
                new_min_pos = (v1_min == -1) ? v2_min : ((v2_min == -1) ? v1_min : (v1_min < v2_min ? v1_min : v2_min));
                
                long v1_max = (max_var_pos != -1) ? max_var_pos : -1;
                long v2_max = (max_var_neg != -1) ? max_var_neg + 1 : -1;
                new_max_pos = (v1_max == -1) ? v2_max : ((v2_max == -1) ? v1_max : (v1_max > v2_max ? v1_max : v2_max));
            }
            
            if (can_be_neg_val) {
                long v1_min = (min_var_neg != -1) ? min_var_neg : -1;
                long v2_min = (min_var_pos != -1) ? min_var_pos + 1 : -1;
                new_min_neg = (v1_min == -1) ? v2_min : ((v2_min == -1) ? v1_min : (v1_min < v2_min ? v1_min : v2_min));
                
                long v1_max = (max_var_neg != -1) ? max_var_neg : -1;
                long v2_max = (max_var_pos != -1) ? max_var_pos + 1 : -1;
                new_max_neg = (v1_max == -1) ? v2_max : ((v2_max == -1) ? v1_max : (v1_max > v2_max ? v1_max : v2_max));
            }
            
            min_var_pos = new_min_pos; max_var_pos = new_max_pos;
            min_var_neg = new_min_neg; max_var_neg = new_max_neg;
        }
    }
    
    long v_min = -1, v_max = -1;
    if (min_var_pos != -1 && min_var_neg != -1) v_min = (min_var_pos < min_var_neg) ? min_var_pos : min_var_neg;
    else if (min_var_pos != -1) v_min = min_var_pos;
    else if (min_var_neg != -1) v_min = min_var_neg;
    
    if (max_var_pos != -1 && max_var_neg != -1) v_max = (max_var_pos > max_var_neg) ? max_var_pos : max_var_neg;
    else if (max_var_pos != -1) v_max = max_var_pos;
    else if (max_var_neg != -1) v_max = max_var_neg;
    
    if (!has_first) {
        *result = ANEWDSC_VAR_ZERO;
    } else {
        if (v_max == 0) {
            *result = ANEWDSC_VAR_ZERO;
        } else if (v_min == 1 && v_max == 1) {
            *result = ANEWDSC_VAR_ONE;
        } else {
            *result = ANEWDSC_VAR_UNKNOWN;
        }
    }
    
    anewdsc_poly_clear(P1); anewdsc_poly_clear(P2); anewdsc_poly_clear(P3); anewdsc_poly_clear(P_I);
    anewdsc_interval_clear(a_iv); anewdsc_interval_clear(b_minus_a_iv); anewdsc_interval_clear(one_iv);
    mpfr_clear(a_fr); mpfr_clear(b_fr);
    
    return ANEWDSC_SUCCESS;
}

anewdsc_status_t anewdsc_find_admissible_point(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t I,
    anewdsc_interval_t m_star,
    int *out_found
) {
    if (!poly || !I || !m_star || !out_found) return ANEWDSC_ERR_INVALID_ARG;
    
    mpfr_t a, b, w, m;
    mpfr_init(a); mpfr_init(b); mpfr_init(w); mpfr_init(m);
    
    mpfi_get_left(a, I->bounds);
    mpfi_get_right(b, I->bounds);
    mpfr_sub(w, b, a, MPFR_RNDN); // w = b - a
    mpfr_add(m, a, b, MPFR_RNDN);
    mpfr_div_ui(m, m, 2, MPFR_RNDN); // m = (a + b) / 2
    
    anewdsc_interval_t eval_res;
    anewdsc_interval_init(eval_res);
    
    anewdsc_status_t status = ANEWDSC_SUCCESS;
    int found = 0;
    
    // Try up to 100 deterministic pseudo-random points (Weyl sequence)
    for (int iter = 0; iter < 100; iter++) {
        // Generate pseudo-random offset in [-0.25 * w, 0.25 * w] using Golden Ratio
        double frac = (iter * 0.6180339887498948482);
        frac -= (long)frac;
        double r = frac * 0.5 - 0.25;
        
        mpfr_t offset, mj;
        mpfr_init(offset); mpfr_init(mj);
        
        mpfr_mul_d(offset, w, r, MPFR_RNDN);
        mpfr_add(mj, m, offset, MPFR_RNDN);
        
        mpfi_set_fr(m_star->bounds, mj);
        
        anewdsc_poly_eval_interval(poly, m_star, eval_res);
        
        // Check if 0 is in eval_res
        if (!mpfi_has_zero(eval_res->bounds)) {
            found = 1;
            mpfr_clear(offset); mpfr_clear(mj);
            break;
        }
        
        mpfr_clear(offset); mpfr_clear(mj);
    }
    
    if (!found) {
        // Fallback to exact midpoint if not found
        mpfi_set_fr(m_star->bounds, m);
    }
    
    *out_found = found;
    
    anewdsc_interval_clear(eval_res);
    mpfr_clear(a); mpfr_clear(b); mpfr_clear(w); mpfr_clear(m);
    
    return status;
}

anewdsc_status_t anewdsc_newton_test(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t I,
    size_t N_I,
    anewdsc_interval_t I_prime,
    int *success
) {
    if (!poly || !I || !I_prime || !success) return ANEWDSC_ERR_INVALID_ARG;
    *success = 0;
    
    mpfr_t a, b, w;
    mpfr_init(a); mpfr_init(b); mpfr_init(w);
    mpfi_get_left(a, I->bounds);
    mpfi_get_right(b, I->bounds);
    mpfr_sub(w, b, a, MPFR_RNDN);
    
    anewdsc_poly_t deriv;
    anewdsc_poly_init(deriv);
    anewdsc_poly_derivative(deriv, poly);
    
    anewdsc_interval_t xi[3], xi_star[3], v[3];
    for (int i=0; i<3; i++) {
        anewdsc_interval_init(xi[i]);
        anewdsc_interval_init(xi_star[i]);
        anewdsc_interval_init(v[i]);
        
        mpfr_t tmp; mpfr_init(tmp);
        mpfr_mul_ui(tmp, w, i+1, MPFR_RNDN);
        mpfr_div_ui(tmp, tmp, 4, MPFR_RNDN);
        mpfr_add(tmp, a, tmp, MPFR_RNDN);
        mpfi_set_fr(xi[i]->bounds, tmp);
        mpfr_clear(tmp);
        
        int dummy_found;
        anewdsc_find_admissible_point(poly, xi[i], xi_star[i], &dummy_found);
        
        anewdsc_interval_t p_val, dp_val;
        anewdsc_interval_init(p_val); anewdsc_interval_init(dp_val);
        anewdsc_poly_eval_interval(poly, xi_star[i], p_val);
        anewdsc_poly_eval_interval(deriv, xi_star[i], dp_val);
        
        if (!mpfi_has_zero(dp_val->bounds)) {
            mpfi_div(v[i]->bounds, p_val->bounds, dp_val->bounds);
        } else {
            mpfi_set_ui(v[i]->bounds, 0);
        }
        
        anewdsc_interval_clear(p_val); anewdsc_interval_clear(dp_val);
    }
    
    for (int i = 0; i < 3 && !(*success); i++) {
        for (int j = i + 1; j < 3 && !(*success); j++) {
            anewdsc_interval_t num, den, k_tilde, lambda;
            anewdsc_interval_init(num); anewdsc_interval_init(den); 
            anewdsc_interval_init(k_tilde); anewdsc_interval_init(lambda);
            
            mpfi_sub(num->bounds, xi_star[j]->bounds, xi_star[i]->bounds);
            mpfi_sub(den->bounds, v[i]->bounds, v[j]->bounds);
            
            if (!mpfi_has_zero(den->bounds)) {
                mpfi_div(k_tilde->bounds, num->bounds, den->bounds);
                mpfi_mul(lambda->bounds, k_tilde->bounds, v[i]->bounds);
                mpfi_add(lambda->bounds, xi_star[i]->bounds, lambda->bounds);
                
                mpfr_t lam_left, lam_right;
                mpfr_init(lam_left); mpfr_init(lam_right);
                mpfi_get_left(lam_left, lambda->bounds);
                mpfi_get_right(lam_right, lambda->bounds);
                
                if (mpfr_cmp(lam_right, a) >= 0 && mpfr_cmp(lam_left, b) <= 0) {
                    mpfr_t w_prime; mpfr_init(w_prime);
                    mpfr_div_ui(w_prime, w, N_I, MPFR_RNDN);
                    
                    mpfr_t m_lam; mpfr_init(m_lam);
                    mpfr_add(m_lam, lam_left, lam_right, MPFR_RNDN);
                    mpfr_div_ui(m_lam, m_lam, 2, MPFR_RNDN);
                    
                    mpfr_t a_prime, b_prime;
                    mpfr_init(a_prime); mpfr_init(b_prime);
                    
                    mpfr_t half_w; mpfr_init(half_w);
                    mpfr_div_ui(half_w, w_prime, 2, MPFR_RNDN);
                    mpfr_sub(a_prime, m_lam, half_w, MPFR_RNDN);
                    mpfr_add(b_prime, m_lam, half_w, MPFR_RNDN);
                    
                    if (mpfr_cmp(a_prime, a) < 0) mpfr_set(a_prime, a, MPFR_RNDN);
                    if (mpfr_cmp(b_prime, b) > 0) mpfr_set(b_prime, b, MPFR_RNDN);
                    
                    anewdsc_interval_t I_prime_test;
                    anewdsc_interval_init(I_prime_test);
                    mpfi_interv_fr(I_prime_test->bounds, a_prime, b_prime);
                    
                    anewdsc_interval_t left_outer; anewdsc_interval_init(left_outer);
                    mpfi_interv_fr(left_outer->bounds, a, a_prime);
                    
                    anewdsc_interval_t right_outer; anewdsc_interval_init(right_outer);
                    mpfi_interv_fr(right_outer->bounds, b_prime, b);
                    
                    anewdsc_var_result_t left_res, right_res;
                    anewdsc_01_test(poly, left_outer, &left_res);
                    anewdsc_01_test(poly, right_outer, &right_res);
                    
                    if (left_res == ANEWDSC_VAR_ZERO && right_res == ANEWDSC_VAR_ZERO) {
                        *success = 1;
                        mpfi_set(I_prime->bounds, I_prime_test->bounds);
                    }
                    
                    anewdsc_interval_clear(left_outer); anewdsc_interval_clear(right_outer);
                    anewdsc_interval_clear(I_prime_test);
                    mpfr_clear(w_prime); mpfr_clear(m_lam); mpfr_clear(a_prime); mpfr_clear(b_prime); mpfr_clear(half_w);
                }
                mpfr_clear(lam_left); mpfr_clear(lam_right);
            }
            anewdsc_interval_clear(num); anewdsc_interval_clear(den);
            anewdsc_interval_clear(k_tilde); anewdsc_interval_clear(lambda);
        }
    }
    
    for (int i=0; i<3; i++) {
        anewdsc_interval_clear(xi[i]);
        anewdsc_interval_clear(xi_star[i]);
        anewdsc_interval_clear(v[i]);
    }
    
    anewdsc_poly_clear(deriv);
    mpfr_clear(a); mpfr_clear(b); mpfr_clear(w);
    
    return ANEWDSC_SUCCESS;
}

typedef struct {
    anewdsc_interval_t I;
    size_t N_I;
} isolate_task_t;

anewdsc_status_t anewdsc_isolate_roots(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t search_domain,
    anewdsc_interval_t **output_roots,
    size_t *num_roots
) {
    if (!poly || !search_domain || !output_roots || !num_roots) {
        return ANEWDSC_ERR_INVALID_ARG;
    }
    
    size_t capacity = 16;
    *output_roots = malloc(capacity * sizeof(anewdsc_interval_t));
    if (!*output_roots) return ANEWDSC_ERR_NOMEM;
    *num_roots = 0;
    
    size_t stack_cap = 32;
    isolate_task_t *stack = malloc(stack_cap * sizeof(isolate_task_t));
    if (!stack) {
        free(*output_roots);
        return ANEWDSC_ERR_NOMEM;
    }
    size_t stack_size = 0;
    
    anewdsc_interval_init(stack[stack_size].I);
    mpfi_set(stack[stack_size].I->bounds, search_domain->bounds);
    stack[stack_size].N_I = 4;
    stack_size++;
    
    anewdsc_status_t status = ANEWDSC_SUCCESS;
    
    while (stack_size > 0) {
        stack_size--;
        isolate_task_t current;
        anewdsc_interval_init(current.I);
        mpfi_set(current.I->bounds, stack[stack_size].I->bounds);
        current.N_I = stack[stack_size].N_I;
        anewdsc_interval_clear(stack[stack_size].I);
        
        mpfr_t a_chk, b_chk;
        mpfr_init(a_chk); mpfr_init(b_chk);
        mpfi_get_left(a_chk, current.I->bounds);
        mpfi_get_right(b_chk, current.I->bounds);
        if (mpfr_cmp(a_chk, b_chk) >= 0) {
            if (*num_roots >= capacity) {
                capacity *= 2;
                anewdsc_interval_t *new_out = realloc(*output_roots, capacity * sizeof(anewdsc_interval_t));
                if (!new_out) { status = ANEWDSC_ERR_NOMEM; mpfr_clear(a_chk); mpfr_clear(b_chk); break; }
                *output_roots = new_out;
            }
            anewdsc_interval_init((*output_roots)[*num_roots]);
            mpfi_set((*output_roots)[*num_roots]->bounds, current.I->bounds);
            (*num_roots)++;
            anewdsc_interval_clear(current.I);
            mpfr_clear(a_chk); mpfr_clear(b_chk);
            continue;
        }
        mpfr_clear(a_chk); mpfr_clear(b_chk);
        
        anewdsc_var_result_t v_res;
        status = anewdsc_01_test(poly, current.I, &v_res);
        if (status != ANEWDSC_SUCCESS) break;
        
        if (v_res == ANEWDSC_VAR_ZERO) {
            anewdsc_interval_clear(current.I);
            continue;
        } else if (v_res == ANEWDSC_VAR_ONE) {
            if (*num_roots >= capacity) {
                capacity *= 2;
                anewdsc_interval_t *new_out = realloc(*output_roots, capacity * sizeof(anewdsc_interval_t));
                if (!new_out) { status = ANEWDSC_ERR_NOMEM; break; }
                *output_roots = new_out;
            }
            anewdsc_interval_init((*output_roots)[*num_roots]);
            mpfi_set((*output_roots)[*num_roots]->bounds, current.I->bounds);
            (*num_roots)++;
            anewdsc_interval_clear(current.I);
            continue;
        } else {
            anewdsc_interval_t I_prime;
            anewdsc_interval_init(I_prime);
            int newton_success = 0;
            
            status = anewdsc_newton_test(poly, current.I, current.N_I, I_prime, &newton_success);
            if (status != ANEWDSC_SUCCESS) break;
            
            if (newton_success) {
                if (stack_size >= stack_cap) {
                    stack_cap *= 2;
                    isolate_task_t *new_stack = realloc(stack, stack_cap * sizeof(isolate_task_t));
                    if (!new_stack) { status = ANEWDSC_ERR_NOMEM; break; }
                    stack = new_stack;
                }
                anewdsc_interval_init(stack[stack_size].I);
                mpfi_set(stack[stack_size].I->bounds, I_prime->bounds);
                size_t next_N_I = current.N_I * current.N_I;
                if (next_N_I < current.N_I || next_N_I > 1000000000UL) {
                    next_N_I = 1000000000UL;
                }
                stack[stack_size].N_I = next_N_I;
                stack_size++;
            } else {
                anewdsc_interval_t m_star;
                anewdsc_interval_init(m_star);
                int found_admissible = 0;
                status = anewdsc_find_admissible_point(poly, current.I, m_star, &found_admissible);
                if (status != ANEWDSC_SUCCESS) break;
                
                mpfr_t a, b, m_val;
                mpfr_init(a); mpfr_init(b); mpfr_init(m_val);
                mpfi_get_left(a, current.I->bounds);
                mpfi_get_right(b, current.I->bounds);
                
                mpfi_get_left(m_val, m_star->bounds);
                
                if (!found_admissible || mpfr_cmp(m_val, a) <= 0 || mpfr_cmp(m_val, b) >= 0) {
                    if (*num_roots >= capacity) {
                        capacity *= 2;
                        anewdsc_interval_t *new_out = realloc(*output_roots, capacity * sizeof(anewdsc_interval_t));
                        if (!new_out) {
                            status = ANEWDSC_ERR_NOMEM;
                            mpfr_clear(a); mpfr_clear(b); mpfr_clear(m_val);
                            anewdsc_interval_clear(m_star);
                            break;
                        }
                        *output_roots = new_out;
                    }
                    anewdsc_interval_init((*output_roots)[*num_roots]);
                    mpfi_set((*output_roots)[*num_roots]->bounds, current.I->bounds);
                    (*num_roots)++;
                    
                    mpfr_clear(a); mpfr_clear(b); mpfr_clear(m_val);
                    anewdsc_interval_clear(m_star);
                    anewdsc_interval_clear(I_prime);
                    anewdsc_interval_clear(current.I);
                    continue;
                }
                
                if (stack_size + 2 > stack_cap) {
                    stack_cap *= 2;
                    isolate_task_t *new_stack = realloc(stack, stack_cap * sizeof(isolate_task_t));
                    if (!new_stack) { status = ANEWDSC_ERR_NOMEM; break; }
                    stack = new_stack;
                }
                
                anewdsc_interval_init(stack[stack_size].I);
                mpfi_interv_fr(stack[stack_size].I->bounds, m_val, b);
                stack[stack_size].N_I = current.N_I;
                stack_size++;
                
                anewdsc_interval_init(stack[stack_size].I);
                mpfi_interv_fr(stack[stack_size].I->bounds, a, m_val);
                stack[stack_size].N_I = current.N_I;
                stack_size++;
                
                mpfr_clear(a); mpfr_clear(b); mpfr_clear(m_val);
                anewdsc_interval_clear(m_star);
            }
            anewdsc_interval_clear(I_prime);
            anewdsc_interval_clear(current.I);
        }
    }
    
    for (size_t i = 0; i < stack_size; i++) {
        anewdsc_interval_clear(stack[i].I);
    }
    free(stack);
    
    if (status != ANEWDSC_SUCCESS) {
        anewdsc_free_roots_array(*output_roots, *num_roots);
        *output_roots = NULL;
        *num_roots = 0;
    }
    
    return status;
}

anewdsc_status_t anewdsc_first_root_in_range(
    const anewdsc_poly_t poly,
    const anewdsc_interval_t search_domain,
    anewdsc_interval_t first_root,
    int *found
) {
    if (!poly || !search_domain || !first_root || !found) {
        return ANEWDSC_ERR_INVALID_ARG;
    }
    
    *found = 0;
    
    size_t stack_cap = 128;
    isolate_task_t stack_buffer[128];
    isolate_task_t *stack = stack_buffer;
    
    size_t stack_size = 0;
    anewdsc_interval_init(stack[stack_size].I);
    mpfi_set(stack[stack_size].I->bounds, search_domain->bounds);
    stack[stack_size].N_I = 4;
    stack_size++;
    
    anewdsc_status_t status = ANEWDSC_SUCCESS;
    
    while (stack_size > 0) {
        stack_size--;
        isolate_task_t current;
        anewdsc_interval_init(current.I);
        mpfi_set(current.I->bounds, stack[stack_size].I->bounds);
        current.N_I = stack[stack_size].N_I;
        anewdsc_interval_clear(stack[stack_size].I);
        
        mpfr_t a_chk, b_chk;
        mpfr_init(a_chk); mpfr_init(b_chk);
        mpfi_get_left(a_chk, current.I->bounds);
        mpfi_get_right(b_chk, current.I->bounds);
        if (mpfr_cmp(a_chk, b_chk) >= 0) {
            mpfi_set(first_root->bounds, current.I->bounds);
            *found = 1;
            anewdsc_interval_clear(current.I);
            mpfr_clear(a_chk); mpfr_clear(b_chk);
            break;
        }
        mpfr_clear(a_chk); mpfr_clear(b_chk);
        
        anewdsc_var_result_t v_res;
        status = anewdsc_01_test(poly, current.I, &v_res);
        if (status != ANEWDSC_SUCCESS) {
            anewdsc_interval_clear(current.I);
            break;
        }
        
        if (v_res == ANEWDSC_VAR_ZERO) {
            anewdsc_interval_clear(current.I);
            continue;
        } else if (v_res == ANEWDSC_VAR_ONE) {
            mpfi_set(first_root->bounds, current.I->bounds);
            *found = 1;
            anewdsc_interval_clear(current.I);
            break; // Stop! We found the first root.
        } else {
            anewdsc_interval_t I_prime;
            anewdsc_interval_init(I_prime);
            int newton_success = 0;
            
            status = anewdsc_newton_test(poly, current.I, current.N_I, I_prime, &newton_success);
            if (status != ANEWDSC_SUCCESS) {
                anewdsc_interval_clear(I_prime);
                anewdsc_interval_clear(current.I);
                break;
            }
            
            if (newton_success) {
                if (stack_size >= stack_cap) {
                    status = ANEWDSC_ERR_NOMEM;
                    break;
                }
                anewdsc_interval_init(stack[stack_size].I);
                mpfi_set(stack[stack_size].I->bounds, I_prime->bounds);
                size_t next_N_I = current.N_I * current.N_I;
                if (next_N_I < current.N_I || next_N_I > 1000000000UL) {
                    next_N_I = 1000000000UL;
                }
                stack[stack_size].N_I = next_N_I;
                stack_size++;
            } else {
                anewdsc_interval_t m_star;
                anewdsc_interval_init(m_star);
                int found_admissible = 0;
                status = anewdsc_find_admissible_point(poly, current.I, m_star, &found_admissible);
                if (status != ANEWDSC_SUCCESS) {
                    anewdsc_interval_clear(m_star);
                    anewdsc_interval_clear(I_prime);
                    anewdsc_interval_clear(current.I);
                    break;
                }
                
                mpfr_t a, b, m_val;
                mpfr_init(a); mpfr_init(b); mpfr_init(m_val);
                mpfi_get_left(a, current.I->bounds);
                mpfi_get_right(b, current.I->bounds);
                mpfi_get_left(m_val, m_star->bounds);
                
                if (!found_admissible || mpfr_cmp(m_val, a) <= 0 || mpfr_cmp(m_val, b) >= 0) {
                    mpfi_set(first_root->bounds, current.I->bounds);
                    *found = 1;
                    mpfr_clear(a); mpfr_clear(b); mpfr_clear(m_val);
                    anewdsc_interval_clear(m_star);
                    anewdsc_interval_clear(I_prime);
                    anewdsc_interval_clear(current.I);
                    break;
                }
                
                if (stack_size + 2 > stack_cap) {
                    status = ANEWDSC_ERR_NOMEM;
                    break;
                }
                
                // Push right interval [m_val, b]
                anewdsc_interval_init(stack[stack_size].I);
                mpfi_interv_fr(stack[stack_size].I->bounds, m_val, b);
                stack[stack_size].N_I = current.N_I;
                stack_size++;
                
                // Push left interval [a, m_val] (Processed first!)
                anewdsc_interval_init(stack[stack_size].I);
                mpfi_interv_fr(stack[stack_size].I->bounds, a, m_val);
                stack[stack_size].N_I = current.N_I;
                stack_size++;
                
                mpfr_clear(a); mpfr_clear(b); mpfr_clear(m_val);
                anewdsc_interval_clear(m_star);
            }
            anewdsc_interval_clear(I_prime);
            anewdsc_interval_clear(current.I);
        }
    }
    
    for (size_t i = 0; i < stack_size; i++) {
        anewdsc_interval_clear(stack[i].I);
    }
    
    return status;
}

anewdsc_status_t anewdsc_free_roots_array(anewdsc_interval_t *output_roots, size_t num_roots) {
    if (output_roots) {
        for (size_t i = 0; i < num_roots; i++) {
            anewdsc_interval_clear(output_roots[i]);
        }
        free(output_roots);
    }
    return ANEWDSC_SUCCESS;
}
