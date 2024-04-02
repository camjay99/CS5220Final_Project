#include <cmath>
#include <cstdio>
#include "parameters.h"
#include "radiation_solver.h"

// As a first order operation, we are going to begin with a single cohort and implement a simple
// multi-layer two-stream model that considers one "big flat leaf" and the ground for three different bands. Future versions
// should incorporate additional cohorts that are each treated as seperate "big flat leaves."
// Considerations for evapotranspiration, growth, diurnal cycles, etc. will be implemented in the future.

int num_cohorts = 10;    // Number of cohorts to include in a grid cell

int main(int argc, char** argv) {
    double* direct_profile_PAR = new double[num_cohorts+1];
    double* direct_profile_NIR = new double[num_cohorts+1];
    double* direct_profile_TIR = new double[num_cohorts+1];
    double* up_diffuse_profile_PAR = new double[num_cohorts+1];
    double* up_diffuse_profile_NIR = new double[num_cohorts+1];
    double* up_diffuse_profile_TIR = new double[num_cohorts+1];
    double* down_diffuse_profile_PAR = new double[num_cohorts+1];
    double* down_diffuse_profile_NIR = new double[num_cohorts+1];
    double* down_diffuse_profile_TIR = new double[num_cohorts+1];
    double* black_body_profile_PAR = new double[num_cohorts+1];
    double* black_body_profile_NIR = new double[num_cohorts+1];
    double* black_body_profile_TIR = new double[num_cohorts+1];
    double* temp_profile = new double[num_cohorts+1];

    // Calculate direct radiation profile 
    calculate_direct_profile(direct_profile_PAR, num_cohorts, incoming_direct_PAR);
    calculate_direct_profile(direct_profile_NIR, num_cohorts, incoming_direct_NIR);
    calculate_direct_profile(direct_profile_TIR, num_cohorts, incoming_direct_TIR);

    // Set temperature
    for (int k = 0; k <= num_cohorts; k++) {
        temp_profile[k] = 0;//298.15;
    }

    // Calculate black body radiation for TIR
    calculate_black_body_profile(black_body_profile_TIR, temp_profile, num_cohorts);

    // Calculate diffuse radiation profiles
    calculate_diffuse_profile(direct_profile_NIR, up_diffuse_profile_NIR, down_diffuse_profile_NIR, black_body_profile_NIR,
                              0, scat_NIR, scat_ground_NIR, backscat_NIR_direct, backscat_NIR_diffuse,
                              inv_op_depth_direct, inv_op_depth_diffuse, num_cohorts, incoming_diffuse_NIR);

    calculate_diffuse_profile(direct_profile_PAR, up_diffuse_profile_PAR, down_diffuse_profile_PAR, black_body_profile_PAR,
                              0, scat_PAR, scat_ground_PAR, backscat_PAR_direct, backscat_PAR_diffuse,
                              inv_op_depth_direct, inv_op_depth_diffuse, num_cohorts, incoming_diffuse_PAR);

    calculate_diffuse_profile(direct_profile_TIR, up_diffuse_profile_TIR, down_diffuse_profile_TIR, black_body_profile_TIR,
                              0, scat_TIR, scat_ground_TIR, backscat_TIR_direct, backscat_TIR_diffuse,
                              inv_op_depth_direct, inv_op_depth_diffuse, num_cohorts, incoming_diffuse_TIR);

    double* absorbed_radiance = new double[num_cohorts];
    calculate_absorbed_radiance(absorbed_radiance, 
                                direct_profile_PAR, direct_profile_NIR, direct_profile_TIR,
                                num_cohorts);

    printf("Direct PAR");
    for (int i = 0; i <= num_cohorts; i++) {
        printf("%f\n", direct_profile_PAR[i]);
    }

    printf("Up Diffuse PAR");
    for (int i = 0; i <= num_cohorts; i++) {
        printf("%f\n", up_diffuse_profile_PAR[i]);
    }

    printf("Down Diffuse PAR");
    for (int i = 0; i <= num_cohorts; i++) {
        printf("%f\n", down_diffuse_profile_PAR[i]);
    }

    printf("Blackbody TIR");
    for (int i = 0; i <= num_cohorts; i++) {
        printf("%f\n", black_body_profile_PAR[i]);
    }

    printf("Total Absorbed PAR");
    for (int i = 0; i <= num_cohorts; i++) {
        printf("%f\n", absorbed_radiance[i]);
    }

    

    
    return 1;
}