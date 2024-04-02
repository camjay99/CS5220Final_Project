#ifndef __RADIATION_SOLVER_H__
#define __RADIATION_SOLVER_H__

// Direct radiation profile solvers

// Solves PAR direct radiation profile. Profile array should be of length num_cohorts+1.
void calculate_direct_profile(double* direct_profile_PAR, int num_cohorts, double incoming);

// Calculate blackbody radiation based on temperature
void calculate_black_body_profile(double* black_body_profile, double* temp_profile, int num_cohorts);

// Indirect radiation profile solvers (to be implemented)
void calculate_diffuse_profile(double* direct_profile, double* up_diffuse_profile, double* down_diffuse_profile, double* black_body_profile, double black_body_ground,
                               double scat, double scat_ground, 
                               double backscat_dir, double backscat_dif, double inv_op_depth_dir, double inv_op_depth_dif, 
                               int num_cohorts, double incoming_dif);

// Uses radiation profiles to determines cohort level absorbed radiance. Absorbed radiance array should be of length num_cohorts.
void calculate_absorbed_radiance(double* absorbed_radiance, 
                                 double* direct_PAR_profile, double* direct_NIR_profile, double* direct_TIR_profile, 
                                 int num_cohorts);
#endif