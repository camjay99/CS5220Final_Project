#include <cmath>
#include <cstdio>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cuda.h>
#include "parameters.h"
#include "radiation_solver.h"
#include "photosynthesis.h"
#include "canopy_air_space.h"

#define NUM_THREADS 256
int blks;

// As a first order operation, we are going to begin with a single cohort and implement a simple
// multi-layer two-stream model that considers one "big flat leaf" and the ground for three different bands. Future versions
// should incorporate additional cohorts that are each treated as seperate "big flat leaves."
// Considerations for evapotranspiration, growth, diurnal cycles, etc. will be implemented in the future.

int num_patches;    // Number of patches to include in a grid cell
int num_cohorts_per_patch;
int seed = 42;
bool print_output = true;

__device__ double dt = 300;

__device__ double* direct_profile_PAR_dev;
__device__ double* direct_profile_NIR_dev;
__device__ double* absorbed_radiance_dev;

double* temp_profile_dev;
double* leaf_area_profile_dev;
double* mass_profile_dev;

// PFT and related stats.
__device__ double wind_speed = 2;
__device__ double air_density = 1.225;
__device__ double G_Wl = 0.01;
__device__ double Vcmax15 = 6.25;
__device__ double Q10Vcmax = 2.40;
__device__ double f_cold = 0.40;
__device__ double f_hot = 0.40;
__device__ double T_cold = 283.15;
__device__ double T_hot = 318.15;
__device__ double f_clump = 0.80;
__device__ double f_R = 0.015;
__device__ double quantum_yield = 0.08;
__device__ double p_c = 1;
__device__ double stomata = 1;
__device__ double soil_temp = 298.15;
__device__ double q_leaf = 3218;
__device__ double leaf_size = 0.1;
__device__ double M = 9;
__device__ double c_c = 400;
__device__ double dw = 0.016;
__device__ double w_c = 0.017;

__device__ double calc_air_temp(double t) {
    return (max_air_temp - min_air_temp) / 2 * std::sin(2 * 3.1415 * t / 86400)  + (max_air_temp + min_air_temp) / 2;
}

__device__ double calc_incoming_PAR(double t) {
    return incoming_direct_PAR / 2 * std::sin(2 * 3.1415 * (t - 7200) / 86400 )  + incoming_direct_PAR / 2;
}

__device__ double calc_incoming_NIR(double t) {
    return incoming_direct_NIR / 2 * std::sin(2 * 3.1415 * (t - 7200) / 86400)  + incoming_direct_NIR/ 2;
}

__device__ double calc_temp_increment(double temp, double incoming_radiation, double incoming_PAR, double air_temp, double leaf_area, double mass) {
    double G_Qlambda = calculate_G_Qlambda(air_temp, temp, leaf_size, wind_speed);
    double G_Wlambda = calculate_G_Wlambda(G_Qlambda, air_density);
    double Q_cohort_canopy = calculate_Q_cohort_canopy(leaf_area, air_density, G_Qlambda, air_temp, temp, w_c);
    double G_Clambda = G_Wlambda/f_Glambda;

    double A = co2_mixing_ratio_solver_C4(G_Wl, G_Wlambda, G_Clambda, M, c_c, w_c, dw, temp, Vcmax15, Q10Vcmax,
                                           f_cold, f_hot, T_cold, T_hot,
                                           incoming_PAR, f_clump, f_R, quantum_yield, p_c);

    double water_flux = calculate_water_flux(stomata, leaf_area, G_Wlambda, G_Wl, w_c, temp, p_c);
    double soil_cohort_enthalpy = calculate_soil_cohort_enthalpy(water_flux, soil_temp);
    double cohort_CAS_enthalpy = calculate_cohort_CAS_enthalpy(water_flux, temp);

    double total_enthalpy_change = (incoming_radiation - Q_cohort_canopy) + (soil_cohort_enthalpy - cohort_CAS_enthalpy);

    return ((total_enthalpy_change) / (mass*((0.7*q_l_water + q_leaf)/1.7)));
}


__global__ void simulate_one_step(int num_patches, int num_cohorts_per_patch, int i, double* temp_profile, double* leaf_area_profile, double* mass_profile) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= num_patches * num_cohorts_per_patch) {
        return;
    }
    //printf("%d ", tid);
    int p = tid / num_cohorts_per_patch;
    int k = tid % num_cohorts_per_patch;
    //printf("%d %d ", p, k);
    //printf("%f ", temp_profile[0]);
        
    //for (int p = 0; p < num_patches; p++) {
        double* direct_profile_PAR = direct_profile_PAR_dev;
        double* direct_profile_NIR = direct_profile_NIR_dev;
        double* absorbed_radiance = absorbed_radiance_dev;
        // Calculate direct radiation profile 
        calculate_direct_profile(direct_profile_PAR, num_cohorts_per_patch, calc_incoming_PAR(i*dt));
        calculate_direct_profile(direct_profile_NIR, num_cohorts_per_patch, calc_incoming_NIR(i*dt));
        calculate_absorbed_radiance(absorbed_radiance, 
                                direct_profile_PAR, direct_profile_NIR,
                                num_cohorts_per_patch);

        // For each cohort in this patch
        //for (int k = 0; k < num_cohorts_per_patch; k++) {
            double k1 = calc_temp_increment(temp_profile[p*num_cohorts_per_patch + k],           absorbed_radiance[k], direct_profile_PAR[k+1] - direct_profile_PAR[k], calc_air_temp(dt*i), 
                                leaf_area_profile[p*num_cohorts_per_patch + k], mass_profile[p*num_cohorts_per_patch + k]);

            double k2 = calc_temp_increment(temp_profile[p*num_cohorts_per_patch + k] + k1*dt/2, absorbed_radiance[k], direct_profile_PAR[k+1] - direct_profile_PAR[k], calc_air_temp(dt*i + dt/2), 
                                leaf_area_profile[p*num_cohorts_per_patch + k], mass_profile[p*num_cohorts_per_patch + k]);

            double k3 = calc_temp_increment(temp_profile[p*num_cohorts_per_patch + k] + k2*dt/2, absorbed_radiance[k], direct_profile_PAR[k+1] - direct_profile_PAR[k], calc_air_temp(dt*i + dt/2), 
                                leaf_area_profile[p*num_cohorts_per_patch + k], mass_profile[p*num_cohorts_per_patch + k]);

            double k4 = calc_temp_increment(temp_profile[p*num_cohorts_per_patch + k] + k3*dt,   absorbed_radiance[k], direct_profile_PAR[k+1] - direct_profile_PAR[k], calc_air_temp(dt*i + dt), 
                                leaf_area_profile[p*num_cohorts_per_patch + k], mass_profile[p*num_cohorts_per_patch + k]);

            temp_profile[p*num_cohorts_per_patch + k] += dt/6 * (k1 + 2*k2 + 2*k3 + k4);

        //}

    //}

}


int main(int argc, char** argv) {
    // As a temporary starting point, we can assume that each patch has ~500 thin cohorts, if we want in the future we can also consider different plant functional types so that parameters can be variable.
    if (argc < 3) {
        printf("Error: please argue num_patches and num_cohorts_per_patch (in that order)\n");
        return -1;
    }
    num_patches = atoi((argv[1]));
    num_cohorts_per_patch = atoi((argv[2]));

    blks = (num_patches * num_cohorts_per_patch + NUM_THREADS - 1) / NUM_THREADS;

    // Initalize patch structure to random values
    std::random_device rd;
    std::mt19937 gen(seed ? seed : rd());
    std::uniform_real_distribution<> la_rd(3., 5.);
    std::uniform_real_distribution<> m_rd(90., 110.);

    double* leaf_area_profile = new double[num_patches*num_cohorts_per_patch];
    double* mass_profile = new double[num_patches*num_cohorts_per_patch];
    for (int k = 0; k < num_patches*num_cohorts_per_patch; k++) {
        leaf_area_profile[k] = la_rd(gen);
        mass_profile[k] = m_rd(gen);
        //printf("%d\t%f\t%f\n", k, leaf_area_profile[k], mass_profile[k]);
    }

    // No need to randomize temp_profiles, as the will quickly converge to long-term behavior.
    double* temp_profile = new double[num_patches*num_cohorts_per_patch];
    
    // Set temperature
    for (int k = 0; k < num_patches*num_cohorts_per_patch; k++) {
        temp_profile[k] = 298.15 + k;
    }


    // alloc GPU space
    cudaMalloc((void **)&direct_profile_PAR_dev,       (num_cohorts_per_patch + 1) * sizeof(double) );
    cudaMalloc((void **)&direct_profile_NIR_dev,       (num_cohorts_per_patch + 1) * sizeof(double) );
    cudaMalloc((void **)&absorbed_radiance_dev ,       (num_cohorts_per_patch)     * sizeof(double) );

    cudaMalloc((void **)&temp_profile_dev ,       (num_cohorts_per_patch)     * sizeof(double) );
    cudaMemcpy(temp_profile_dev, temp_profile, num_patches * num_cohorts_per_patch * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&leaf_area_profile_dev ,       (num_cohorts_per_patch)     * sizeof(double) );
    cudaMemcpy(leaf_area_profile_dev, leaf_area_profile, num_patches * num_cohorts_per_patch * sizeof(double), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&mass_profile_dev ,       (num_cohorts_per_patch)     * sizeof(double) );
    cudaMemcpy(mass_profile_dev, mass_profile, num_patches * num_cohorts_per_patch * sizeof(double), cudaMemcpyHostToDevice);
	
    // Starting simulation algorithm
    auto start_time = std::chrono::stead_clock::now()

    // for each time step
    for (int i = 0; i < 2016; i++) {

        //printf("%d", i);

        simulate_one_step<<<blks, NUM_THREADS>>>(num_patches, num_cohorts_per_patch, i, temp_profile_dev, leaf_area_profile_dev, mass_profile_dev);

    }

    //printf("\n%d\n", blks);
    //printf("%d\n", NUM_THREADS);

    cudaMemcpy(temp_profile, temp_profile_dev, num_patches * num_cohorts_per_patch * sizeof(double), cudaMemcpyDeviceToHost); // copy data back from gpu

    printf("\n%f %f %f %f\n", temp_profile[0], temp_profile[1], temp_profile[2], temp_profile[3]);

    auto end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double> diff = end_time - start_time;
    double seconds = diff.count();

    std::cout << "Simulation Time = " << seconds << " seconds for " << num_patches
                  << " patches with " << num_cohorts_per_patch << " cohorts per patch.\n";
    
    return 0;
}
