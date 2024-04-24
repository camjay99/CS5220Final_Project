#ifndef __NUMERICS_H__
#define __NUMERICS_H__

#include <functional>

__device__ double bisection(std::function<double(double)> func, double begin, double end, double tol, int num_iter);

#endif