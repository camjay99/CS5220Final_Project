#include <cmath>
#include "parameters.h"
#include "radiation_solver.h"
#include "../Eigen/Dense"

// This document contain functions for solving radiation transfer for PAR, NIR, and TIR

// Direct radiation profile solvers

// Solves direct radiation profile. Profile array should be of length num_cohorts+1.
void calculate_direct_profile(double* direct_profile, int num_cohorts, double incoming) {
    direct_profile[num_cohorts] = incoming;
    for (int k = num_cohorts-1; k >= 0; k--) {
        direct_profile[k] = direct_profile[k+1]*exp(-total_PAI / inv_op_depth_direct);
    }
}

// IN THE FUTURE: Calculate diffuse vertical profile (requires sparse linear algebra)
// A whole bunch of helper functions for setting up the diffuse vertical profile solver

double D_plus(double scat, double backscat) {
    return 0.5 * (1 + sqrt((1 - scat) / (1 - (1 - 2*backscat)*scat)));
}

double D_minus(double scat, double backscat) {
    return 0.5 * (1 - sqrt((1 - scat) / (1 - (1 - 2*backscat)*scat)));
}

double xi_sq(double scat, double backscat, double inv_op_depth) {
    return (1 - (1 - 2*backscat)*scat)*(1 - scat) / (inv_op_depth*inv_op_depth);
}

double kappa_plus(double scat, double backscat_diffuse, double inv_op_depth_dif, double backscat_direct, double inv_op_depth_dir, double interception) {
    return -((1 - (1 - 2*backscat_diffuse)*scat)/inv_op_depth_dif +
              (1 - 2*backscat_direct) / inv_op_depth_dir) *
            (scat*interception/inv_op_depth_dir);
}

double kappa_minus(double scat, double backscat_diffuse, double inv_op_depth_dif, double backscat_direct, double inv_op_depth_dir, double interception) {
    return -((1 - scat)*(1 - 2*backscat_diffuse)/inv_op_depth_dif +
              1 / inv_op_depth_dir) *
            (scat*interception/inv_op_depth_dir);
}

double p_plus(double scat, double backscat_diffuse, double inv_op_depth_dif, double backscat_direct, double inv_op_depth_dir, double interception) {
    return (kappa_plus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception) +
            kappa_minus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception)) * inv_op_depth_dir * inv_op_depth_dir /
            (2 * (1 - xi_sq(scat, backscat_diffuse, inv_op_depth_diffuse) * inv_op_depth_dir * inv_op_depth_dir));
}

double p_plus(double scat, double backscat_diffuse, double inv_op_depth_dif, double backscat_direct, double inv_op_depth_dir, double interception) {
    return (kappa_plus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception) -
            kappa_minus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception)) * inv_op_depth_dir * inv_op_depth_dir /
            (2 * (1 - xi_sq(scat, backscat_diffuse, inv_op_depth_diffuse) * inv_op_depth_dir * inv_op_depth_dir));
}


// Uses radiation profiles to determines cohort level absorbed radiance. Absorbed radiance array should be of length num_cohorts.
void calculate_absorbed_radiance(double* absorbed_radiance, double* direct_PAR_profile, double* direct_NIR_profile, double* direct_TIR_profile, int num_cohorts) {
    for (int k = 0; k < num_cohorts; k++) {
        absorbed_radiance[k] = direct_PAR_profile[k+1] - direct_PAR_profile[k] +
                               direct_NIR_profile[k+1] - direct_NIR_profile[k] +
                               direct_TIR_profile[k+1] - direct_TIR_profile[k];
    }
}