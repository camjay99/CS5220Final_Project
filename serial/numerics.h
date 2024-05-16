#ifndef __NUMERICS_H__
#define __NUMERICS_H__

#include <functional>

//double bisection(std::function<double(double)> func, double begin, double end, double tol, int num_iter);
double bisection(double* roots, double* values, double* begin, double* end, double tol, int num_iters,
                 double f_gl, double f_glambda, double G_wl, double* G_wlambda, double gamma,
                 double c_c, double w_c, double* w_l, double d_w, double m, double* F_A, double* F_B, double* F_C, double* F_D);
#endif