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

// Uses radiation profiles to determines cohort level absorbed radiance. Absorbed radiance array should be of length num_cohorts.
void calculate_absorbed_radiance(double* absorbed_radiance, double* direct_PAR_profile, double* direct_NIR_profile, double* direct_TIR_profile, int num_cohorts) {
    for (int k = 0; k < num_cohorts; k++) {
        absorbed_radiance[k] = direct_PAR_profile[k+1] - direct_PAR_profile[k] +
                               direct_NIR_profile[k+1] - direct_NIR_profile[k] +
                               direct_TIR_profile[k+1] - direct_TIR_profile[k];
    }
}

// NOTE: CURRENTLY NOT PLANNING ON USING THE FOLLOWING FUNCTIONS AS DIFFUSE RADIATION CAN BE SIMPLY MODELLED AS INCREASED EXTINCTION COEFFICIENT.
//       MAY DELETE THESE FUNCTION IN THE FUTURE.

void calculate_black_body_profile(double* black_body_profile, double* temp_profile, int num_cohorts) {
    for (int k = 0; k < num_cohorts; k++) {
        black_body_profile[k] = stef_boltz_constant*pow(temp_profile[k], 4);
    }
    black_body_profile[num_cohorts + 1] = 0;
}

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
    return -((1 - scat)*(1 - 2*backscat_direct)/inv_op_depth_dif +
              1 / inv_op_depth_dir) *
            (scat*interception/inv_op_depth_dir);
}

double p_plus(double scat, double backscat_diffuse, double inv_op_depth_dif, double backscat_direct, double inv_op_depth_dir, double interception) {
    return (kappa_plus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception) +
            kappa_minus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception)) * inv_op_depth_dir * inv_op_depth_dir /
            (2 * (1 - xi_sq(scat, backscat_diffuse, inv_op_depth_dif) * inv_op_depth_dir * inv_op_depth_dir));
}

double p_minus(double scat, double backscat_diffuse, double inv_op_depth_dif, double backscat_direct, double inv_op_depth_dir, double interception) {
    return (kappa_plus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception) -
            kappa_minus(scat, backscat_diffuse, inv_op_depth_dif, backscat_direct, inv_op_depth_dir, interception)) * inv_op_depth_dir * inv_op_depth_dir /
            (2 * (1 - xi_sq(scat, backscat_diffuse, inv_op_depth_dif) * inv_op_depth_dir * inv_op_depth_dir));
}

void calculate_diffuse_profile(double* direct_profile, double* up_diffuse_profile, double* down_diffuse_profile, double* black_body_profile, double black_body_ground,
                               double scat, double scat_ground, 
                               double backscat_dir, double backscat_dif, double inv_op_depth_dir, double inv_op_depth_dif, 
                               int num_cohorts, double incoming_dif) {
    Eigen::MatrixXd S(2*num_cohorts + 2, 2*num_cohorts + 2);
    Eigen::VectorXd y(2*num_cohorts + 2);
    Eigen::VectorXd x(2*num_cohorts + 2);
    
    // Fill in matrix S
    double D_plus_m_k;
    double D_minus_m_k;
    double D_plus_m_k1;
    double D_minus_m_k1;
    double pos_mod;
    double neg_mod;
    S(0, 0) = (D_minus(scat, backscat_dif) - scat_ground*D_plus(scat, backscat_dif))*
                exp(-sqrt(xi_sq(scat, backscat_dif, inv_op_depth_dif))*total_PAI);
    S(0, 1) = (D_plus(scat, backscat_dif) - scat_ground*D_minus(scat, backscat_dif))*
                exp(sqrt(xi_sq(scat, backscat_dif, inv_op_depth_dif))*total_PAI);
    for (int k = 1; k <= num_cohorts; k++) {
        D_plus_m_k = D_plus(scat, backscat_dif);
        D_minus_m_k = D_minus(scat, backscat_dif);
        D_plus_m_k1 = (k == num_cohorts ? 1 : D_plus(scat, backscat_dif));
        D_minus_m_k1 = (k == num_cohorts ? 0 : D_minus(scat, backscat_dif));
        neg_mod = exp(-sqrt((k == num_cohorts ? xi_sq(1, 0, 1) : xi_sq(scat, backscat_dif, inv_op_depth_dif)))*(k == num_cohorts ? 0 : total_PAI));
        pos_mod = exp(sqrt((k == num_cohorts ? xi_sq(1, 0, 1) : xi_sq(scat, backscat_dif, inv_op_depth_dif)))*(k == num_cohorts ? 0 : total_PAI));

        S(2*k - 1, 2*k - 2) = D_plus_m_k;
        S(2*k - 1, 2*k - 1) = D_minus_m_k;
        S(2*k - 1, 2*k)     = -D_plus_m_k1*neg_mod;
        S(2*k - 1, 2*k + 1) = -D_minus_m_k1*pos_mod;
        S(2*k, 2*k - 2)     = D_minus_m_k;
        S(2*k, 2*k - 1)     = D_plus_m_k;
        S(2*k, 2*k)         = -D_minus_m_k1*neg_mod;
        S(2*k + 1, 2*k + 1) = -D_plus_m_k1*pos_mod;
    }
    S(2*(num_cohorts+1) - 1, 2*(num_cohorts+1) - 2) = 1;
    S(2*(num_cohorts+1) - 1, 2*(num_cohorts+1) - 1) = 0;
    //S(2*(num_cohorts+1), 2*(num_cohorts+1) - 2)     = 0;
    //S(2*(num_cohorts+1), 2*(num_cohorts+1) - 1)     = 1;

    // Fill in vector y
    double p_m = p_minus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[1]);
    double p_p = p_plus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[1]);
    y(0) = scat_ground*direct_profile[0] + (1-scat_ground)*(black_body_ground - black_body_profile[0]) -
           (p_m - scat_ground*p_p)*exp(-total_PAI / inv_op_depth_dir);
    for (int k = 1; k <= num_cohorts; k++) {
        p_m = p_minus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[k]);
        p_p = p_plus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[k]);
        y(2*k - 1) = p_p*exp(-total_PAI / inv_op_depth_dir) -
                     p_p + black_body_profile[k] - black_body_profile[k-1];
        y(2*k)     = p_m*exp(-total_PAI / inv_op_depth_dir) -
                     p_m + black_body_profile[k] - black_body_profile[k-1];
    }
    p_p = p_plus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[num_cohorts]);
    y(2*num_cohorts + 1) = incoming_dif - p_p;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(S);
    x = dec.solve(y);

    // Use x coefficients to calculate upward and downward diffuse radiation.
    for (int k = 1; k <= num_cohorts + 1; k++) {
        double D_plus_m_k = D_plus(scat, backscat_dif);
        double D_minus_m_k = D_minus(scat, backscat_dif);
        double neg_mod = exp(-sqrt(xi_sq(scat, backscat_dif, inv_op_depth_dif))*total_PAI);
        double pos_mod = exp(sqrt(xi_sq(scat, backscat_dif, inv_op_depth_dif))*total_PAI);
        double p_m = p_minus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[k]);
        double p_p = p_plus(scat, backscat_dif, inv_op_depth_dif, backscat_dir, inv_op_depth_dir, direct_profile[k]);
        down_diffuse_profile[k-1] = x(2*k - 2)*D_plus_m_k*neg_mod + 
                                    x(2*k - 1)*D_minus_m_k*pos_mod +
                                    p_p*exp(-total_PAI/inv_op_depth_dir) + black_body_profile[k-1];
        up_diffuse_profile[k-1]   = x(2*k - 2)*D_minus_m_k*neg_mod + 
                                    x(2*k - 1)*D_plus_m_k*pos_mod +
                                    p_m*exp(-total_PAI/inv_op_depth_dir) + black_body_profile[k-1];

    }
}
