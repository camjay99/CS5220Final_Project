// Define any basic numerical functions needed for the program
#include <cassert>
#include <functional>

// A simple bisection method for root-finding.
double bisection(std::function<double(double)> func, double begin, double end, double tol, int num_iters) {
    int begin_sign = (begin > 0) - (begin < 0);
    int end_sign = (end > 0) - (end < 0);

    // Signs must be different to ensure existence of solution in the interval.
    if (begin_sign != end_sign)
        return -1;

    double mid = (end - begin) / 2;
    double eval = func(mid);
    double i = 0;
    while (eval > tol) {
        // Update endpoints of search interval based on the sign of the current midpoint
        if (((mid > 0) - (mid < 0)) == begin_sign) {
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