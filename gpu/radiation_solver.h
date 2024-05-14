#ifndef __RADIATION_SOLVER_H__
#define __RADIATION_SOLVER_H__

// Direct radiation profile solvers

// Solves PAR direct radiation profile. Profile array should be of length num_cohorts+1.
__device__ void calculate_direct_profile(double* direct_profile, int num_cohorts, double incoming);

// Calculate blackbody radiation based on temperature
__device__ void calculate_black_body_profile(double* black_body_profile, double* temp_profile, int num_cohorts);

// Uses radiation profiles to determines cohort level absorbed radiance. Absorbed radiance array should be of length num_cohorts.
__device__ void calculate_absorbed_radiance(double* absorbed_radiance, 
                                 double* direct_PAR_profile, double* direct_NIR_profile,
                                 int num_cohorts);
#endif