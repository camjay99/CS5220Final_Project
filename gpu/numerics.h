#ifndef __NUMERICS_H__
#define __NUMERICS_H__

#include <functional>
template <typename Function>
__device__ double bisection(Function func, double begin, double end, double tol, int num_iter);

#endif