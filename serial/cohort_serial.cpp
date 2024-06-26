#include <cmath>
#include <cstdio>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include "parameters.h"
#include "radiation_solver.h"
#include "photosynthesis.h"
#include "canopy_air_space.h"

// As a first order operation, we are going to begin with a single cohort and implement a simple
// multi-layer two-stream model that considers one "big flat leaf" and the ground for three different bands. Future versions
// should incorporate additional cohorts that are each treated as seperate "big flat leaves."
// Considerations for evapotranspiration, growth, diurnal cycles, etc. will be implemented in the future.

int num_patches = 1;    // Number of patches to include in a grid cell
int num_cohorts_per_patch = 4;
int seed = 42;
double dt = 300;
bool print_output = true;

// PFT and related stats.
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

double calc_air_temp(double t) {
    return (max_air_temp - min_air_temp) / 2 * std::sin(2 * 3.1415 * t / 86400)  + (max_air_temp + min_air_temp) / 2;
}

double calc_incoming_PAR(double t) {
    return incoming_direct_PAR / 2 * std::sin(2 * 3.1415 * (t - 18000) / 86400 )  + incoming_direct_PAR / 2;
}

double calc_incoming_NIR(double t) {
    return incoming_direct_NIR / 2 * std::sin(2 * 3.1415 * (t - 18000) / 86400)  + incoming_direct_NIR/ 2;
}

// Calculate temperature increments four at a time.
void calc_temp_increment(double* d_temp, double* temp, double* incoming_radiation, double* incoming_PAR, double air_temp, double* leaf_area, double* mass) {
    double* G_Qlambda = new double[4];
    double* G_Wlambda = new double[4];
    double* Q_cohort_canopy = new double[4];
    double* G_Clambda = new double[4];
    for (int i = 0; i < 4; i++) {
        G_Qlambda[i] = calculate_G_Qlambda(air_temp, temp[i], leaf_size, wind_speed);
        G_Wlambda[i] = calculate_G_Wlambda(G_Qlambda[i], air_density);
        Q_cohort_canopy[i] = calculate_Q_cohort_canopy(leaf_area[i], air_density, G_Qlambda[i], air_temp, temp[i], w_c);
        G_Clambda[i] = G_Wlambda[i]/f_Glambda;
    }

    double* A = new double[4];
    co2_mixing_ratio_solver_C4(A, G_Wl, 
                               G_Wlambda, G_Clambda, M, 
                               c_c, w_c, dw,
                               temp, Vcmax15, Q10Vcmax,
                               f_cold, f_hot, T_cold, T_hot,
                               incoming_PAR, f_clump, f_R, quantum_yield, p_c);

    double water_flux;
    double soil_cohort_enthalpy;
    double cohort_CAS_enthalpy;
    double total_enthalpy_change;

    for (int i = 0; i < 4; i++) {
        water_flux = calculate_water_flux(stomata, leaf_area[i], G_Wlambda[i], G_Wl, w_c, temp[i], p_c);
        soil_cohort_enthalpy = calculate_soil_cohort_enthalpy(water_flux, soil_temp);
        cohort_CAS_enthalpy = calculate_cohort_CAS_enthalpy(water_flux, temp[i]);

        total_enthalpy_change = (incoming_radiation[i] - Q_cohort_canopy[i]) + (soil_cohort_enthalpy - cohort_CAS_enthalpy);
        d_temp[i] = ((total_enthalpy_change) / (mass[i]*((0.7*q_l_water + q_leaf)/1.7)));
    } 
}

int main(int argc, char** argv) {
    // As a temporary starting point, we can assume that each patch has ~500 thin cohorts, if we want in the future we can also consider different plant functional types so that parameters can be variable.

    // Initalize patch structure to random values
    std::random_device rd;
    std::mt19937 gen(seed ? seed : rd());
    std::uniform_real_distribution la_rd(3., 5.);
    std::uniform_real_distribution m_rd(90., 110.);

    double* leaf_area_profile = new double[num_patches*num_cohorts_per_patch];
    double* mass_profile = new double[num_patches*num_cohorts_per_patch];
    for (int k = 0; k < num_patches*num_cohorts_per_patch; k++) {
        leaf_area_profile[k] = la_rd(gen);
        mass_profile[k] = m_rd(gen);
    }

    // No need to randomize temp_profiles, as the will quickly converge to long-term behavior.
    double* temp_profile = new double[num_patches*num_cohorts_per_patch];
    
    // Set temperature
    for (int k = 0; k < num_patches*num_cohorts_per_patch; k++) {
        temp_profile[k] = 298.15;
    }

    // Starting simulation algorithm
    auto start_time = std::chrono::steady_clock::now();

    // For each patch
    for (int p = 0; p < num_patches; p++) {
        double* direct_profile_PAR = new double[num_cohorts_per_patch+1];
        double* direct_profile_NIR = new double[num_cohorts_per_patch+1];
        double* absorbed_radiance = new double[num_cohorts_per_patch];
        double* absorbed_PAR = new double[num_cohorts_per_patch];
        // Order could probably be changed to be more sensical, but this loops over all time steps.
        std::ofstream output;
        if ((p == 0) && print_output)
            output.open("C:/Users/camer/CS5220Final_Project/output.txt");
        for (int i = 0; i < 2016; i++) {
            // Calculate direct radiation profile 
            calculate_direct_profile(direct_profile_PAR, num_cohorts_per_patch, calc_incoming_PAR(i*dt));
            calculate_direct_profile(direct_profile_NIR, num_cohorts_per_patch, calc_incoming_NIR(i*dt));
            calculate_absorbed_radiance(absorbed_radiance, 
                                    direct_profile_PAR, direct_profile_NIR,
                                    num_cohorts_per_patch);
            for (int k = 0; k < num_cohorts_per_patch; k++) {
                absorbed_PAR[k] = direct_profile_PAR[k+1] - direct_profile_PAR[k];
            }
            // For each cohort in this patch
            double* k1 = new double[4];
            double* k2 = new double[4];
            double* k3 = new double[4];
            double* k4 = new double[4];
            double* intermediate_temp = new double[4];
            for (int k = 0; k < num_cohorts_per_patch; k += 4) {
                
                calc_temp_increment(k1, temp_profile + p*num_cohorts_per_patch + k, absorbed_radiance + k, absorbed_PAR + k, calc_air_temp(dt*i), 
                                    leaf_area_profile + p*num_cohorts_per_patch + k, mass_profile + p*num_cohorts_per_patch + k);

                for (int i = 0; i < 4; i++) {
                    intermediate_temp[i] = temp_profile[p*num_cohorts_per_patch + k + i] + k1[i]*dt/2;
                }
                calc_temp_increment(k2, intermediate_temp, absorbed_radiance + k, absorbed_PAR + k, calc_air_temp(dt*i + dt/2), 
                                    leaf_area_profile + p*num_cohorts_per_patch + k, mass_profile + p*num_cohorts_per_patch + k);

                for (int i = 0; i < 4; i++) {
                    intermediate_temp[i] = temp_profile[p*num_cohorts_per_patch + k + i] + k2[i]*dt/2;
                }
                calc_temp_increment(k3, intermediate_temp, absorbed_radiance + k, absorbed_PAR + k, calc_air_temp(dt*i + dt/2), 
                                    leaf_area_profile + p*num_cohorts_per_patch + k, mass_profile + p*num_cohorts_per_patch + k);

                for (int i = 0; i < 4; i++) {
                    intermediate_temp[i] = temp_profile[p*num_cohorts_per_patch + k + i] + k3[i]*dt;
                }
                calc_temp_increment(k4, intermediate_temp, absorbed_radiance + k, absorbed_PAR + k, calc_air_temp(dt*i + dt), 
                                    leaf_area_profile + p*num_cohorts_per_patch + k, mass_profile + p*num_cohorts_per_patch + k);

                for (int i = 0; i < 4; i++) {
                    temp_profile[p*num_cohorts_per_patch + k + i] += dt/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
                }

                if ((p == 0) && print_output)
                    for (int i = 0; i < 4; i++) {
                        output << temp_profile[k + i] << "\t";
                    }
            }
            if ((p == 0) && print_output)
                    output << "\n";
        }
        if ((p == 0) && print_output)
                output.close();

    }

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    std::cout << "Simulation Time = " << seconds << " seconds for " << num_patches
                  << " patches with " << num_cohorts_per_patch << " cohorts per patch.\n";

    return 0;
}
