// Define any basic numerical functions needed for the program
#include <cassert>
#include <functional>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <immintrin.h>

// A simple bisection method for root-finding.
/*double bisection(std::function<double(double)> func, double begin, double end, double tol, int num_iters) {
    int begin_sign = (func(begin) > 0) - (func(begin) < 0);
    int end_sign = (func(end) > 0) - (func(end) < 0);
    //printf("%d %d\n", begin_sign, end_sign);
    // Signs must be different to ensure existence of solution in the interval.
    if (begin_sign != end_sign)
        return -1;

    double mid = (end - begin) / 2;
    double eval = func(mid);
    double i = 0;
    while (eval > tol) {
        // Update endpoints of search interval based on the sign of the current midpoint
        if (((eval > 0) - (eval < 0)) == begin_sign) {
            begin = mid;
        } else {
            end = mid;
        }

        // Update midpoint
        mid = (end - begin) / 2; 
        eval = func(mid);

        // If no solution has been found after a set number of iterations, assume none exist.
        i++;
        if (i > num_iters) 
            return -1;
    }

    return mid;
}*/

// Calculate F based on register level operations. Final answers are returned in buffer_1
void calc_F_register(__m256d* buffer_1, __m256d* buffer_2, __m256d* buffer_3, __m256d* buffer_4,
                     __m256d* r_f_Gl, 
                     __m256d* r_f_Glambda,
                     __m256d* r_G_wl,
                     __m256d* r_G_wlambda,
                     __m256d* r_gamma,
                     __m256d* r_c_c,
                     __m256d* r_lead,
                     __m256d* r_m,
                     __m256d* r_one,
                     __m256d* r_c_l,
                     __m256d* r_A) {
    // F_2
    *buffer_1 = _mm256_sub_pd(*r_c_l, *r_gamma);
    *buffer_1 = _mm256_mul_pd(*r_G_wlambda, *buffer_1);
    *buffer_2 = _mm256_mul_pd(*r_f_Glambda, *r_A);
    *buffer_1 = _mm256_sub_pd(*buffer_1, *buffer_2);
    *buffer_2 = _mm256_sub_pd(*r_c_c, *r_c_l);
    *buffer_2 = _mm256_mul_pd(*r_G_wlambda, *buffer_2);
    *buffer_3 = _mm256_sub_pd(*r_f_Gl, *r_f_Glambda);
    *buffer_3 = _mm256_mul_pd(*r_A, *buffer_3);
    *buffer_2 = _mm256_sub_pd(*buffer_2, *buffer_3);
    *buffer_1 = _mm256_div_pd(*buffer_1, *buffer_2);

    // F_3
    *buffer_2 = _mm256_sub_pd(*r_c_c, *r_c_l);
    *buffer_2 = _mm256_mul_pd(*r_G_wlambda, *buffer_2);
    *buffer_3 = _mm256_mul_pd(*r_f_Glambda, *r_A);
    *buffer_3 = _mm256_sub_pd(*buffer_2, *buffer_3);
    *buffer_4 = _mm256_sub_pd(*r_f_Gl, *r_f_Glambda);
    *buffer_4 = _mm256_mul_pd(*buffer_4, *r_A);
    *buffer_2 = _mm256_add_pd(*buffer_2, *buffer_4);
    *buffer_2 = _mm256_div_pd(*buffer_3, *buffer_2);
    *buffer_2 = _mm256_mul_pd(*r_lead, *buffer_2);
    *buffer_2 = _mm256_add_pd(*r_one, *buffer_2);

    // F_1
    *buffer_3 = _mm256_div_pd(*r_G_wl, *r_G_wlambda);
    *buffer_3 = _mm256_mul_pd(*r_f_Glambda, *buffer_3);
    *buffer_3 = _mm256_sub_pd(*r_f_Gl, *buffer_3);
    *buffer_3 = _mm256_mul_pd(*buffer_3, *r_A);
    *buffer_4 = _mm256_add_pd(*r_c_c, *r_c_l);
    *buffer_4 = _mm256_mul_pd(*r_G_wl, *buffer_4);
    *buffer_3 = _mm256_sub_pd(*buffer_3, *buffer_4);
    *buffer_4 = _mm256_mul_pd(*r_m, *r_A);
    *buffer_3 = _mm256_div_pd(*buffer_3, *buffer_4);

    // F
    *buffer_1 = _mm256_mul_pd(*buffer_1, *buffer_2);
    *buffer_1 = _mm256_mul_pd(*buffer_1, *buffer_3);
    *buffer_1 = _mm256_sub_pd(*buffer_1, *r_one);
}

void bisection(double* roots, double* values, double* begin, double* end, double tol, int num_iters,
                 double f_gl, double f_glambda, double G_wl, double* G_wlambda, double gamma,
                 double c_c, double w_c, double* w_l, double d_w, double m, double* F_A, double* F_B, double* F_C, double* F_D) {
    // Takes a set of four values and simultaneously compute bisection.
    __m256d r_f_Gl = _mm256_set1_pd(f_gl);
    __m256d r_f_Glambda = _mm256_set1_pd(f_glambda);
    __m256d r_G_wl = _mm256_set1_pd(G_wl);
    __m256d r_G_wlambda = _mm256_loadu_pd(G_wlambda);
    __m256d r_gamma = _mm256_set1_pd(gamma);
    __m256d r_c_c = _mm256_set1_pd(c_c);

    __m256d r_lead = _mm256_set1_pd(w_c);
    r_lead = _mm256_sub_pd(r_lead, _mm256_loadu_pd(w_l));
    r_lead = _mm256_div_pd(r_lead, _mm256_set1_pd(d_w));

    __m256d r_m = _mm256_set1_pd(m);
    __m256d r_one = _mm256_set1_pd(1);


    __m256d buffer_1;
    __m256d buffer_2;
    __m256d buffer_3;
    __m256d buffer_4;

    // Here we calculate c_lk and A_k
    __m256d r_c_l;
    __m256d r_A;

    // Determine initial boundary signs
    r_c_l = _mm256_loadu_pd(begin);
    __m256d r_begin_signs;
    calc_F_register(&r_begin_signs, &buffer_2, &buffer_3, &buffer_4,
                    &r_f_Gl, 
                    &r_f_Glambda,
                    &r_G_wl,
                    &r_G_wlambda,
                    &r_gamma,
                    &r_c_c,
                    &r_lead,
                    &r_m,
                    &r_one,
                    &r_c_l,
                    &r_A);
    r_begin_signs = _mm256_cmp_pd(r_begin_signs, _mm256_set1_pd(0.0), _CMP_LT_OQ);

    buffer_1 = _mm256_loadu_pd(end);
    buffer_2 = _mm256_set1_pd(2);
    r_c_l = _mm256_div_pd(_mm256_add_pd(r_c_l, buffer_1), buffer_2);

    for (int i = 0; i < num_iters; i++) {
        buffer_1 = _mm256_loadu_pd(F_A);
        buffer_2 = _mm256_loadu_pd(F_B);
        buffer_3 = _mm256_loadu_pd(F_C);
        buffer_4 = _mm256_loadu_pd(F_D);

        buffer_1 = _mm256_mul_pd(buffer_1, r_c_l);
        buffer_3 = _mm256_mul_pd(buffer_3, r_c_l);
        buffer_1 = _mm256_add_pd(buffer_1, buffer_2);
        buffer_3 = _mm256_add_pd(buffer_3, buffer_4);
        r_A = _mm256_div_pd(buffer_1, buffer_2);

        // Now we can calculate the main functions.
        calc_F_register(&buffer_1, &buffer_2, &buffer_3, &buffer_4,
                        &r_f_Gl, 
                        &r_f_Glambda,
                        &r_G_wl,
                        &r_G_wlambda,
                        &r_gamma,
                        &r_c_c,
                        &r_lead,
                        &r_m,
                        &r_one,
                        &r_c_l,
                        &r_A);

        // Once calculated, we must now update c_l based one where the zero is expected to be.
        buffer_2 = _mm256_cmp_pd(buffer_1, _mm256_set1_pd(0.0), _CMP_GE_OQ); //bitmask for nen-negative values
        _mm256_storeu_pd(values, buffer_2);
        buffer_2 = _mm256_xor_pd(buffer_2, r_begin_signs); // Identify which numbers have the same sign as the begin values (or rather begin and estimate do not have opposite signs)
        
        // Load same sign end-points
        buffer_3 = _mm256_load_pd(begin);
        buffer_4 = _mm256_load_pd(end);

        // Save with masked data
        _mm256_maskstore_pd(begin, (__m256i) buffer_2, r_c_l);
        buffer_2 = _mm256_xor_pd(buffer_2, (__m256d) _mm256_set1_epi32(-1));
        _mm256_maskstore_pd(end, (__m256i) r_c_l, r_c_l);

        buffer_3 = _mm256_blendv_pd(buffer_3, buffer_4, buffer_2);

        // Calculate midpoints and save properly into memory
        buffer_3 = _mm256_add_pd(buffer_3, r_c_l);
        r_c_l = _mm256_div_pd(buffer_3, _mm256_set1_pd(2.0));
    }

    _mm256_storeu_pd(roots, r_c_l);

    for (int r = 0; r < 4; r++) {
        if (std::abs(values[r]) > tol) {
            roots[r] = -1;
        }
    }
}

