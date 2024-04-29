// Define any basic numerical functions needed for the program
#include <cassert>
#include <functional>
#include <cstdio>
#include <iostream>
#include <fstream>

// A simple bisection method for root-finding.
template <typename Function>
__device__ double bisection(std::function<double(double)> func, double begin, double end, double tol, int num_iters) {
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
}