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
    // Calculate direct radiation profile 
    double* direct_profile_PAR = new double[num_cohorts+1];
    double* direct_profile_NIR = new double[num_cohorts+1];
    double* direct_profile_TIR = new double[num_cohorts+1];

    calculate_direct_profile(direct_profile_PAR, num_cohorts, incoming_PAR);
    calculate_direct_profile(direct_profile_NIR, num_cohorts, incoming_NIR);
    calculate_direct_profile(direct_profile_TIR, num_cohorts, incoming_TIR);

    double* absorbed_radiance = new double[num_cohorts];
    calculate_absorbed_radiance(absorbed_radiance, direct_profile_PAR, direct_profile_NIR, direct_profile_TIR, num_cohorts);

    for (int i = 0; i < num_cohorts; i++) {
        printf("%f\n", absorbed_radiance[i]);
    }

    

    
    return 1;
}