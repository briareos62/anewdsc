#include "anewdsc.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void poly_mul_linear(anewdsc_poly_t p, double r) {
    size_t old_deg = p->degree;
    anewdsc_poly_t old_p;
    anewdsc_poly_init(old_p);
    anewdsc_poly_copy(old_p, p);
    
    anewdsc_poly_set_coeff_d(p, old_deg + 1, 0.0);
    for(size_t i = 0; i <= old_deg + 1; i++) {
        mpfi_set_ui(p->coeffs[i], 0);
    }
    
    mpfi_t tmp; mpfi_init(tmp);
    for(size_t i = 0; i <= old_deg; i++) {
        mpfi_add(p->coeffs[i+1], p->coeffs[i+1], old_p->coeffs[i]);
        mpfi_mul_d(tmp, old_p->coeffs[i], -r);
        mpfi_add(p->coeffs[i], p->coeffs[i], tmp);
    }
    p->degree = old_deg + 1;
    
    mpfi_clear(tmp);
    anewdsc_poly_clear(old_p);
}

int compare_doubles(const void *a, const void *b) {
    double arg1 = *(const double *)a;
    double arg2 = *(const double *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

int compare_intervals(const void *a, const void *b) {
    anewdsc_interval_t *ia = (anewdsc_interval_t *)a;
    anewdsc_interval_t *ib = (anewdsc_interval_t *)b;
    mpfr_t la, lb;
    mpfr_init(la); mpfr_init(lb);
    mpfi_get_left(la, (*ia)->bounds);
    mpfi_get_left(lb, (*ib)->bounds);
    int res = mpfr_cmp(la, lb);
    mpfr_clear(la); mpfr_clear(lb);
    return res;
}

int main() {
    srand(time(NULL));
    
    for (int test = 1; test <= 10; test++) {
        int num_roots = 1 + rand() % 10; // 1 to 10 roots
        double roots[10];
        
        printf("--- Test %d ---\n", test);
        printf("Generating polynomial with %d roots...\n", num_roots);
        
        anewdsc_poly_t poly;
        anewdsc_poly_init(poly);
        anewdsc_poly_set_coeff_d(poly, 0, 1.0); // P(x) = 1
        
        for (int i = 0; i < num_roots; i++) {
            roots[i] = ((double)rand() / RAND_MAX) * 20.0 - 10.0;
            poly_mul_linear(poly, roots[i]);
        }
        
        qsort(roots, num_roots, sizeof(double), compare_doubles);
        
        printf("Original roots (sorted):\n");
        for (int i = 0; i < num_roots; i++) {
            printf("  Root %d: %f\n", i+1, roots[i]);
        }
        
        anewdsc_interval_t search_domain;
        anewdsc_interval_init(search_domain);
        mpfi_interv_d(search_domain->bounds, -20.0, 20.0);
        
        anewdsc_interval_t *found_roots = NULL;
        size_t num_found = 0;
        
        printf("Running anewdsc_isolate_roots...\n");
        anewdsc_status_t status = anewdsc_isolate_roots(poly, search_domain, &found_roots, &num_found);
        
        if (status != ANEWDSC_SUCCESS) {
            printf("Error in isolate_roots: %d\n", status);
        } else {
            printf("Found %zu root intervals.\n", num_found);
            qsort(found_roots, num_found, sizeof(anewdsc_interval_t), compare_intervals);
            
            for (size_t i = 0; i < num_found; i++) {
                mpfr_t left, right;
                mpfr_init(left); mpfr_init(right);
                mpfi_get_left(left, found_roots[i]->bounds);
                mpfi_get_right(right, found_roots[i]->bounds);
                
                double left_d = mpfr_get_d(left, MPFR_RNDD);
                double right_d = mpfr_get_d(right, MPFR_RNDU);
                
                printf("  Interval %zu: [%.6f, %.6f]\n", i+1, left_d, right_d);
                
                mpfr_clear(left); mpfr_clear(right);
            }
        }
        
        anewdsc_free_roots_array(found_roots, num_found);
        anewdsc_interval_clear(search_domain);
        anewdsc_poly_clear(poly);
        printf("\n");
    }
    
    return 0;
}
