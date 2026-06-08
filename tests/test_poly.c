#include <stdio.h>
#include <assert.h>
#include "anewdsc.h"

void test_eval() {
    anewdsc_poly_t poly;
    anewdsc_interval_t x, res;
    anewdsc_status_t status;
    
    status = anewdsc_poly_init(poly);
    assert(status == ANEWDSC_SUCCESS);
    status = anewdsc_interval_init(x);
    assert(status == ANEWDSC_SUCCESS);
    status = anewdsc_interval_init(res);
    assert(status == ANEWDSC_SUCCESS);
    
    // P(x) = 2 + 3x + x^2
    anewdsc_poly_set_coeff_d(poly, 0, 2.0);
    anewdsc_poly_set_coeff_d(poly, 1, 3.0);
    anewdsc_poly_set_coeff_d(poly, 2, 1.0);
    
    // Evaluate at x = 2
    mpfi_set_d(x->bounds, 2.0);
    anewdsc_poly_eval_interval(poly, x, res);
    
    // Expected: 2 + 6 + 4 = 12
    mpfr_t left, right;
    mpfr_init(left); mpfr_init(right);
    mpfi_get_left(left, res->bounds);
    mpfi_get_right(right, res->bounds);
    assert(mpfr_cmp_d(left, 12.0) == 0);
    assert(mpfr_cmp_d(right, 12.0) == 0);
    
    // Evaluate at interval [1, 2]
    mpfi_interv_d(x->bounds, 1.0, 2.0);
    anewdsc_poly_eval_interval(poly, x, res);
    // P(1) = 6, P(2) = 12. Since coeffs are positive, interval is [6, 12]
    mpfi_get_left(left, res->bounds);
    mpfi_get_right(right, res->bounds);
    assert(mpfr_cmp_d(left, 6.0) == 0);
    assert(mpfr_cmp_d(right, 12.0) == 0);
    
    mpfr_clear(left); mpfr_clear(right);
    anewdsc_interval_clear(x);
    anewdsc_interval_clear(res);
    anewdsc_poly_clear(poly);
}

void test_taylor_shift() {
    anewdsc_poly_t poly, shifted;
    anewdsc_interval_t a;
    anewdsc_status_t status;
    
    status = anewdsc_poly_init(poly);
    assert(status == ANEWDSC_SUCCESS);
    status = anewdsc_poly_init(shifted);
    assert(status == ANEWDSC_SUCCESS);
    status = anewdsc_interval_init(a);
    assert(status == ANEWDSC_SUCCESS);
    
    // P(x) = 1 + 2x + 3x^2
    anewdsc_poly_set_coeff_d(poly, 0, 1.0);
    anewdsc_poly_set_coeff_d(poly, 1, 2.0);
    anewdsc_poly_set_coeff_d(poly, 2, 3.0);
    
    // a = 1
    mpfi_set_d(a->bounds, 1.0);
    
    // Q(x) = P(x+1) = 3x^2 + 8x + 6
    status = anewdsc_poly_taylor_shift(shifted, poly, a);
    assert(status == ANEWDSC_SUCCESS);
    assert(shifted->degree == 2);
    
    mpfr_t left, right;
    mpfr_init(left); mpfr_init(right);
    
    mpfi_get_left(left, shifted->coeffs[0]);
    assert(mpfr_cmp_d(left, 6.0) == 0);
    
    mpfi_get_left(left, shifted->coeffs[1]);
    assert(mpfr_cmp_d(left, 8.0) == 0);
    
    mpfi_get_left(left, shifted->coeffs[2]);
    assert(mpfr_cmp_d(left, 3.0) == 0);
    
    mpfr_clear(left); mpfr_clear(right);
    anewdsc_interval_clear(a);
    anewdsc_poly_clear(shifted);
    anewdsc_poly_clear(poly);
}

void test_01_test() {
    anewdsc_poly_t poly;
    anewdsc_interval_t I;
    anewdsc_var_result_t result;
    anewdsc_status_t status;
    
    anewdsc_poly_init(poly);
    anewdsc_interval_init(I);
    
    // P(x) = x^2 - 2
    anewdsc_poly_set_coeff_d(poly, 0, -2.0);
    anewdsc_poly_set_coeff_d(poly, 1, 0.0);
    anewdsc_poly_set_coeff_d(poly, 2, 1.0);
    
    // Test [1, 2], expected ONE
    mpfi_interv_d(I->bounds, 1.0, 2.0);
    status = anewdsc_01_test(poly, I, &result);
    assert(status == ANEWDSC_SUCCESS);
    assert(result == ANEWDSC_VAR_ONE);
    
    // Test [2, 3], expected ZERO
    mpfi_interv_d(I->bounds, 2.0, 3.0);
    status = anewdsc_01_test(poly, I, &result);
    assert(status == ANEWDSC_SUCCESS);
    assert(result == ANEWDSC_VAR_ZERO);
    
    // Test [0, 3], expected ONE
    mpfi_interv_d(I->bounds, 0.0, 3.0);
    status = anewdsc_01_test(poly, I, &result);
    assert(status == ANEWDSC_SUCCESS);
    assert(result == ANEWDSC_VAR_ONE);
    
    // Test [-2, 2], expected UNKNOWN
    mpfi_interv_d(I->bounds, -2.0, 2.0);
    status = anewdsc_01_test(poly, I, &result);
    assert(status == ANEWDSC_SUCCESS);
    assert(result == ANEWDSC_VAR_UNKNOWN);
    
    anewdsc_interval_clear(I);
    anewdsc_poly_clear(poly);
}

void test_admissible_point() {
    anewdsc_poly_t poly;
    anewdsc_interval_t I, m_star;
    
    anewdsc_poly_init(poly);
    anewdsc_interval_init(I);
    anewdsc_interval_init(m_star);
    
    // P(x) = x^2 - 2
    anewdsc_poly_set_coeff_d(poly, 0, -2.0);
    anewdsc_poly_set_coeff_d(poly, 1, 0.0);
    anewdsc_poly_set_coeff_d(poly, 2, 1.0);
    
    mpfi_interv_d(I->bounds, 1.0, 2.0);
    
    int dummy_found;
    anewdsc_status_t status = anewdsc_find_admissible_point(poly, I, m_star, &dummy_found);
    assert(status == ANEWDSC_SUCCESS);
    
    anewdsc_interval_t eval_res;
    anewdsc_interval_init(eval_res);
    anewdsc_poly_eval_interval(poly, m_star, eval_res);
    
    // It should have found a point where P(m_star) does not contain 0
    assert(!mpfi_has_zero(eval_res->bounds));
    
    anewdsc_interval_clear(eval_res);
    anewdsc_interval_clear(m_star);
    anewdsc_interval_clear(I);
    anewdsc_poly_clear(poly);
}

void test_newton_test() {
    anewdsc_poly_t poly;
    anewdsc_interval_t I, I_prime;
    anewdsc_status_t status;
    int success = 0;
    
    status = anewdsc_poly_init(poly);
    assert(status == ANEWDSC_SUCCESS);
    status = anewdsc_interval_init(I);
    assert(status == ANEWDSC_SUCCESS);
    status = anewdsc_interval_init(I_prime);
    assert(status == ANEWDSC_SUCCESS);
    
    // P(x) = x^2 - 3x + 2.25 (root at 1.5, multiplicity 2)
    anewdsc_poly_set_coeff_d(poly, 0, 2.25);
    anewdsc_poly_set_coeff_d(poly, 1, -3.0);
    anewdsc_poly_set_coeff_d(poly, 2, 1.0);
    
    mpfi_interv_d(I->bounds, 1.0, 2.0);
    
    status = anewdsc_newton_test(poly, I, 4, I_prime, &success);
    assert(status == ANEWDSC_SUCCESS);
    
    if (success) {
        // Just verify that the returned interval is valid
        mpfr_t left, right;
        mpfr_init(left); mpfr_init(right);
        mpfi_get_left(left, I_prime->bounds);
        mpfi_get_right(right, I_prime->bounds);
        assert(mpfr_cmp(left, right) <= 0);
        mpfr_clear(left); mpfr_clear(right);
    }
    
    anewdsc_interval_clear(I_prime);
    anewdsc_interval_clear(I);
    anewdsc_poly_clear(poly);
}

void test_isolate_roots() {
    anewdsc_poly_t poly;
    anewdsc_interval_t search_domain;
    anewdsc_status_t status;
    anewdsc_interval_t *roots = NULL;
    size_t num_roots = 0;
    
    anewdsc_poly_init(poly);
    anewdsc_interval_init(search_domain);
    
    // P(x) = x^2 - 2 (roots at -1.414... and 1.414...)
    anewdsc_poly_set_coeff_d(poly, 0, -2.0);
    anewdsc_poly_set_coeff_d(poly, 1, 0.0);
    anewdsc_poly_set_coeff_d(poly, 2, 1.0);
    
    mpfi_interv_d(search_domain->bounds, -2.0, 2.0);
    
    status = anewdsc_isolate_roots(poly, search_domain, &roots, &num_roots);
    assert(status == ANEWDSC_SUCCESS);
    assert(num_roots == 2);
    
    anewdsc_free_roots_array(roots, num_roots);
    anewdsc_interval_clear(search_domain);
    anewdsc_poly_clear(poly);
}

void test_poly_arithmetic() {
    anewdsc_poly_t p1, p2, res, q, r;
    anewdsc_status_t status;
    
    anewdsc_poly_init(p1);
    anewdsc_poly_init(p2);
    anewdsc_poly_init(res);
    anewdsc_poly_init(q);
    anewdsc_poly_init(r);
    
    // p1 = x^2 + 2x + 1
    anewdsc_poly_set_coeff_d(p1, 0, 1.0);
    anewdsc_poly_set_coeff_d(p1, 1, 2.0);
    anewdsc_poly_set_coeff_d(p1, 2, 1.0);
    
    // p2 = x + 1
    anewdsc_poly_set_coeff_d(p2, 0, 1.0);
    anewdsc_poly_set_coeff_d(p2, 1, 1.0);
    
    // ADD: x^2 + 3x + 2
    status = anewdsc_poly_add(res, p1, p2);
    assert(status == ANEWDSC_SUCCESS);
    assert(res->degree == 2);
    
    // SUB: x^2 + x
    status = anewdsc_poly_sub(res, p1, p2);
    assert(status == ANEWDSC_SUCCESS);
    assert(res->degree == 2);
    
    // MUL: x^3 + 3x^2 + 3x + 1
    status = anewdsc_poly_mul(res, p1, p2);
    assert(status == ANEWDSC_SUCCESS);
    assert(res->degree == 3);
    
    // DIVREM: (x^2 + 2x + 1) / (x + 1) = x + 1, remainder 0
    status = anewdsc_poly_divrem(q, r, p1, p2);
    assert(status == ANEWDSC_SUCCESS);
    assert(q->degree == 1);
    assert(r->degree == 0);
    
    // Check if remainder is zero
    assert(mpfi_has_zero(r->coeffs[0]));
    
    anewdsc_poly_clear(p1);
    anewdsc_poly_clear(p2);
    anewdsc_poly_clear(res);
    anewdsc_poly_clear(q);
    anewdsc_poly_clear(r);
}

int main() {
    anewdsc_poly_t poly;
    anewdsc_status_t status;
    
    status = anewdsc_poly_init(poly);
    assert(status == ANEWDSC_SUCCESS);
    
    status = anewdsc_poly_set_coeff_d(poly, 0, -1.0);
    assert(status == ANEWDSC_SUCCESS);
    
    status = anewdsc_poly_set_coeff_str(poly, 2, "1.0", 10);
    assert(status == ANEWDSC_SUCCESS);
    
    assert(poly->degree == 2);
    
    status = anewdsc_poly_clear(poly);
    assert(status == ANEWDSC_SUCCESS);
    
    test_eval();
    test_taylor_shift();
    test_01_test();
    test_admissible_point();
    test_newton_test();
    test_isolate_roots();
    test_poly_arithmetic();
    
    printf("Tests passed successfully!\n");
    return 0;
}
