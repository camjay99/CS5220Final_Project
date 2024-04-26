#include <cmath>
#include "parameters.h"
#include "radiation_solver.h"

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
void calculate_absorbed_radiance(double* absorbed_radiance, double* direct_PAR_profile, double* direct_NIR_profile, int num_cohorts) {
    for (int k = 0; k < num_cohorts; k++) {
        absorbed_radiance[k] = direct_PAR_profile[k+1] - direct_PAR_profile[k] +
                               direct_NIR_profile[k+1] - direct_NIR_profile[k];
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