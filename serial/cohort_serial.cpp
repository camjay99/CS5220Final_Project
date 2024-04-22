#include <cmath>
#include <cstdio>
#include "parameters.h"
#include "radiation_solver.h"
#include "photosynthesis.h"
#include "canopy_air_space.h"

// As a first order operation, we are going to begin with a single cohort and implement a simple
// multi-layer two-stream model that considers one "big flat leaf" and the ground for three different bands. Future versions
// should incorporate additional cohorts that are each treated as seperate "big flat leaves."
// Considerations for evapotranspiration, growth, diurnal cycles, etc. will be implemented in the future.

int num_patches = 10;    // Number of patches to include in a grid cell
int num_cohorts_per_patch = 500;
double air_temp = 300;
double wind_speed = 2;
double air_density = 1.225;
double G_Wl = 0.01;
double Vcmax15 = 6.25;
double Q10Vcmax = 2.40;
double f_cold = 0.40;
double f_hot = 0.40;
double T_cold = 283.15;
double T_hot = 318.15;
double f_clump = 0.80;
double f_R = 0.015;
double quantum_yield = 0.08;
double p_c = 1;
double stomata = 1;
double soil_temp = 298.15;
double q_leaf = 3218;
double leaf_size = 0.1;
double M = 9;
double c_c = 400;
double dw = 0.016;
double w_c = 0.017;
double dt = 0.5;

int main(int argc, char** argv) {
    // As a temporary starting point, we can assume that each patch has ~500 thin cohorts, if we want in the future we can also consider different plant functional types so that parameters can be variable.

    // Initalize patch structure in the future, randomize these
    double* leaf_area_profile = new double[num_patches*num_cohorts_per_patch];
    double* mass_profile = new double[num_patches*num_cohorts_per_patch];
    for (int k = 0; k < num_patches*num_cohorts_per_patch; k++) {
        leaf_area_profile[k] = 4;
        mass_profile[k] = 100;
    }
    double* temp_profile = new double[num_patches*num_cohorts_per_patch];
    
    // For each patch
    for (int p = 0; p < num_patches; p++) {
        double* direct_profile_PAR = new double[num_cohorts_per_patch+1];
        double* direct_profile_NIR = new double[num_cohorts_per_patch+1];
        double* black_body_profile_TIR = new double[num_cohorts_per_patch+1];

        // Calculate direct radiation profile 
        calculate_direct_profile(direct_profile_PAR, num_cohorts_per_patch, incoming_direct_PAR);
        calculate_direct_profile(direct_profile_NIR, num_cohorts_per_patch, incoming_direct_NIR);

        // Set temperature
        for (int k = 0; k <= num_cohorts_per_patch; k++) {
            temp_profile[p*num_cohorts_per_patch + k] = 298.15;
        }

        // Calculate black body radiation for TIR
        //calculate_black_body_profile(black_body_profile_TIR, temp_profile, num_cohorts);

        double* absorbed_radiance = new double[num_cohorts_per_patch];
        calculate_absorbed_radiance(absorbed_radiance, 
                                    direct_profile_PAR, direct_profile_NIR, black_body_profile_TIR,
                                    num_cohorts_per_patch);

        // Order could probably be changed to be more sensical, but this loops over all time steps.
        for (int i = 0; i < 1000; i++) {
            // For each cohort in this patch
            for (int k = 0; k < num_cohorts_per_patch; k++) {
                double G_Qlambda = calculate_G_Qlambda(air_temp, temp_profile[p*num_cohorts_per_patch + k], leaf_size, wind_speed);
                double G_Wlambda = calculate_G_Wlambda(G_Qlambda, air_density);
                double Q_cohort_canopy = calculate_Q_cohort_canopy(leaf_area_profile[p*num_cohorts_per_patch + k], air_density, G_Qlambda, air_temp, temp_profile[p*num_cohorts_per_patch + k], w_c);
                double G_Clambda = G_Wlambda/f_Glambda;
                // Calculate photosynthesis
                double A = co2_mixing_ratio_solver_C4(G_Wl, G_Wlambda, G_Clambda, M, c_c, w_c, dw, temp_profile[p*num_cohorts_per_patch + k], Vcmax15, Q10Vcmax,
                                           f_cold, f_hot, T_cold, T_hot,
                                           direct_profile_PAR[k+1] - direct_profile_PAR[k], f_clump, f_R, quantum_yield, p_c);

                double water_flux = calculate_water_flux(stomata, leaf_area_profile[p*num_cohorts_per_patch + k], G_Wlambda, G_Wl, w_c, temp_profile[p*num_cohorts_per_patch + k], p_c);
                double soil_cohort_enthalpy = calculate_soil_cohort_enthalpy(water_flux, soil_temp);
                double cohort_CAS_enthalpy = calculate_cohort_CAS_enthalpy(water_flux, temp_profile[p*num_cohorts_per_patch + k]);
                //printf("A: %f\n", A);
                double total_enthalpy_change = (absorbed_radiance[k] - Q_cohort_canopy) + (soil_cohort_enthalpy - cohort_CAS_enthalpy);
                //printf("%f\t%f\t%f\t%f\n", absorbed_radiance[k], -Q_cohort_canopy, soil_cohort_enthalpy, -cohort_CAS_enthalpy);
                // Update temperature.
                temp_profile[p*num_cohorts_per_patch + k] += ((total_enthalpy_change) / (mass_profile[p*num_cohorts_per_patch + k]*((0.7*q_l_water + q_leaf)/1.7)))*dt;
            }
        }
        //for (int k = 0; k < num_cohorts; k++) {
        //    printf("Cohort %d, Temp %f\n", k, temp_profile[k]);
        //}
    }
    return 1;
}