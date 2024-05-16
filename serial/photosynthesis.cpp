#include <cmath>
#include "parameters.h"
#include "numerics.h"
#include <cstdio>
#include <iostream>
#include <fstream>

// Functions used for calculating photosynthesis of cohorts

double calculate_metabolic_rate(double temp, double x15, double Q10x) {
    return x15*std::pow(Q10x, (temp - 288.15)/10);
}

double calculate_michaelis_constants(double temp) {
    double K_C = calculate_metabolic_rate(temp, C15, Q10C);
    double K_O = calculate_metabolic_rate(temp, O15, Q10O);
    return K_C * (1 + o2_mixing_ratio / K_O);
}


double calculate_Vcmax(double temp, double x15, double Q10x, double f_cold, double f_hot, double T_cold, double T_hot) {
    return calculate_metabolic_rate(temp, x15, Q10x) / 
            ((1 + exp(-f_cold*(temp - T_cold)))*(1 + exp(-f_hot*(temp - T_hot))));
}

double calculate_photon_flux(double absorbed_PAR, double f_clump) {
    return (1 / ein) * (f_clump / total_PAI) * absorbed_PAR;
}

double calculate_leaf_respiration(double f_R, double Vcmax) {
    return f_R*Vcmax;
}

// Note: this function assumes a C3 plant
double calculate_gamma(double temp) {
    double CO_ratio = calculate_metabolic_rate(temp, COratio15, Q10COratio);
    return o2_mixing_ratio / (2*CO_ratio);
}

double calculate_p_sat(double temp) {
    double Y_1 = 54.842763 - 6763.22 / temp - 4.210*std::log(temp) + 0.000367*temp;
    double Y_2 = 53.878 - 1331.22 / temp - 9.44523*std::log(temp) + 0.014925*temp;
    double p_vl = std::exp(Y_1 + Y_2*std::tanh(0.0415*(temp - 218.8)));
    double p_i = std::exp(9.550426 - 5723.265 / temp + 3.53068*std::log(temp) - 0.00728332*temp);

    return (p_i < p_vl ? p_i : p_vl);
}

double calculate_w_sat(double temp, double p) {
    double p_sat = calculate_p_sat(temp);
    return M_w*p_sat / (M_d*(p - p_sat) + M_w*p_sat);
}

// The following functions are used for the CO2 mixing ratio solver, which is important for extracting A' and 
// ultimately calculating photosynthesis with respect to different limiting resources

double uptake(double F_A, double F_B, double F_C, double F_D, double co2_mixing_ratio) {
    return (F_A*co2_mixing_ratio + F_B) / (F_C*co2_mixing_ratio + F_D);
}
/*****************************************************************
double F_1(double co2_mixing_ratio, double G_Wl, double G_Wlambda, double M, double c_c,
           double F_A, double F_B, double F_C, double F_D) {
    double A = uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio);
    return ((f_Gl - f_Glambda*G_Wl/G_Wlambda)*A - G_Wl*(c_c - co2_mixing_ratio)) / (M*A);
}

double F_2(double co2_mixing_ratio, double G_Wl, double G_Wlambda, double M, double gamma, double c_c,
           double F_A, double F_B, double F_C, double F_D) {
    double A = uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio);
    return (G_Wlambda*(c_c - gamma) - f_Glambda*A) / (G_Wlambda*(c_c - co2_mixing_ratio) + (f_Gl - f_Glambda)*A);
}

double F_3(double co2_mixing_ratio, double G_Wl, 
           double G_Wlambda, double M, double c_c, double w_c, double w_l, double dw, 
           double F_A, double F_B, double F_C, double F_D) {
    double A = uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio);
    return 1 + (w_c - w_l) / (dw) * (G_Wlambda*(c_c - co2_mixing_ratio) - f_Glambda*A) / (G_Wlambda*(c_c - co2_mixing_ratio) + (f_Gl - f_Glambda)*A);
}

double F(double co2_mixing_ratio, double G_Wl, 
           double G_Wlambda, double M, double gamma, double c_c, double w_c, double w_l, double dw, 
           double F_A, double F_B, double F_C, double F_D) {
    return F_1(co2_mixing_ratio, G_Wl, G_Wlambda, M, c_c,
           F_A, F_B, F_C, F_D) * 
           F_2(co2_mixing_ratio, G_Wl, G_Wlambda, M, gamma, c_c,
               F_A, F_B, F_C, F_D) * 
           F_3(co2_mixing_ratio, G_Wl, 
               G_Wlambda, M, c_c, w_c, w_l, dw, 
               F_A, F_B, F_C, F_D) - 1;
}

double calculate_A(double F_A, double F_B, double F_C, double F_D,
                   double G_Wl, 
                   double G_Wlambda, double G_Clambda, double M, double gamma, 
                   double c_c, double w_c, double w_l, double dw, double leaf_respiration) {
    // Create lambda function to fill in parameters of F
    auto F_spef = [&G_Wl, &G_Wlambda, &M, &gamma, 
                   &c_c, &w_c, &w_l, &dw, &F_A, &F_B, &F_C, &F_D](double co2_mixing_ratio) 
                    -> double {return F(co2_mixing_ratio, G_Wl, G_Wlambda, M, gamma,
                                        c_c, w_c, w_l, dw, F_A, F_B, F_C, F_D);};

    double co2_mixing_ratio;

    if (F_A == 0 && F_C == 0) {
        co2_mixing_ratio = bisection(F_spef, 0, 5000, 1e-3, 1000);
    } else {
        double c_lmin = -(F_D*M - F_B) / (F_C*M - F_A);

        double c_lmax;
        if (F_C == 0) {
            // Eq S198 reduces to a linear function of c_lmax
            c_lmax = (F_B - F_D*(G_Clambda*c_c + M))/(-F_B - F_D*G_Clambda);
        } else {
            double b = (G_Clambda*F_D + F_B - F_C*(G_Clambda*c_c + M)) / (G_Clambda*F_C);
            double c = (F_B - F_D*(G_Clambda*c_c + M)) / (G_Clambda*F_C);
            c_lmax = (-b - std::sqrt(b*b - 4*c)) / 2;
        }
        c_lmax = (c_lmax < c_c? c_lmax : c_c); // Ensure singularity is within realm of plausible values

        co2_mixing_ratio = bisection(F_spef, c_lmin, c_lmax, 1e-6, 100);
    }

    // If no roots in this range, assume stomata are closed and no photosynthesis is occurring.
    double A = (co2_mixing_ratio > 0 ? uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio) : 0);
    return A - leaf_respiration;
}

// Note... c_c is canopy airspace carbon, and will be assumed a driver.
double co2_mixing_ratio_solver_C4(double G_Wl, 
                                  double G_Wlambda, double G_Clambda, double M, 
                                  double c_c, double w_c, double dw,
                                  double temp, double Vcmax15, double Q10Vcmax,
                                  double f_cold, double f_hot, double T_cold, double T_hot,
                                  double absorbed_PAR, double f_clump, double f_R, double quantum_yield, double p_c) {
    // Calculate relevant parameters for this analysis
    double Vcmax = calculate_Vcmax(temp, Vcmax15, Q10Vcmax, f_cold, f_hot, T_cold, T_hot);
    double photon_flux = calculate_photon_flux(absorbed_PAR, f_clump);
    double leaf_respiration = calculate_leaf_respiration(f_R, Vcmax);
    double w_l = calculate_w_sat(temp, p_c); // In the future, may attempt to model changes in p_c based on atmospheric forcings
    
    // Find the root to S193 for each of the potential limiting cases (RuBisCo, light, and CO2) for C4 plants.
    /// RuBisCo
    double A_rubisco = calculate_A(0, Vcmax, 0, 1,
                                   G_Wl, 
                                   G_Wlambda, G_Clambda, M, 0, 
                                   c_c, w_c, w_l, dw, leaf_respiration);

    /// CO2
    double A_CO2     = calculate_A(k_PEP*Vcmax, 0, 0, 1,
                                   G_Wl, 
                                   G_Wlambda, G_Clambda, M, 0, 
                                   c_c, w_c, w_l, dw, leaf_respiration);

    /// Light
    double A_light   = calculate_A(0, quantum_yield*photon_flux, 0, 1,
                                   G_Wl, 
                                   G_Wlambda, G_Clambda, M, 0, 
                                   c_c, w_c, w_l, dw, leaf_respiration);

    //printf("%f %f %f\n", A_rubisco, A_CO2, A_light);
    // Return the minimum estimate of A to reflect limiting conditions (Farquhar model)
    //printf("Rubisco %f CO2 %f Light %f\n", A_rubisco, A_CO2, A_light);
    if ((A_rubisco < A_CO2) && (A_rubisco < A_light))
        return A_rubisco;
    else if (A_CO2 < A_light)
        return A_CO2;
    else 
        return A_light;
}
*************************************************************************/

void calculate_A(double* A, double* F_A, double* F_B, double* F_C, double* F_D,
                   double G_Wl, 
                   double* G_Wlambda, double* G_Clambda, double M, double gamma, 
                   double c_c, double w_c, double* w_l, double dw, double* leaf_respiration) {

    double co2_mixing_ratio;
    double* values = new double[4];
    double* begin = new double[4];
    double * end = new double[4];
    double* roots = new double[4];
    double b;
    double c;
    for (int i = 0; i < 4; i++) {
        if (F_A[i] == 0 && F_C[i] == 0) {
            begin[i] = 0;
            end[i] = 5000;
        } else {
            begin[i] = -(F_D[i]*M - F_B[i]) / (F_C[i]*M - F_A[i]);
        
            if (F_C[i] == 0) {
            // Eq S198 reduces to a linear function of c_lmax
                end[i] = (F_B[i] - F_D[i]*(G_Clambda[i]*c_c + M))/(-F_B[i] - F_D[i]*G_Clambda[i]);
            } else {
                b = (G_Clambda[i]*F_D[i] + F_B[i] - F_C[i]*(G_Clambda[i]*c_c + M)) / (G_Clambda[i]*F_C[i]);
                c = (F_B[i] - F_D[i]*(G_Clambda[i]*c_c + M)) / (G_Clambda[i]*F_C[i]);
                end[i] = (-b - std::sqrt(b*b - 4*c)) / 2;
            }
        }

        end[i] = (end[i] < c_c ? end[i] : c_c); // Ensure singularity is within realm of plausible values
    }

    co2_mixing_ratio = bisection(roots, values, begin, end, 1e-6, 100,
                 f_Gl, f_Glambda, G_Wl, G_Wlambda, gamma,
                 c_c, w_c, w_l, dw, M, F_A, F_B, F_C, F_D);

    // If no roots in this range, assume stomata are closed and no photosynthesis is occurring.
    for (int i = 0; i < 4; i++) {
        A[i] = (co2_mixing_ratio > 0 ? uptake(F_A[i], F_B[i], F_C[i], F_D[i], co2_mixing_ratio) : 0) - leaf_respiration[i];
    }
}

void co2_mixing_ratio_solver_C4(double* A, double G_Wl, 
                                  double* G_Wlambda, double* G_Clambda, double M, 
                                  double c_c, double w_c, double dw,
                                  double* temp, double Vcmax15, double Q10Vcmax,
                                  double f_cold, double f_hot, double T_cold, double T_hot,
                                  double* absorbed_PAR, double f_clump, double f_R, double quantum_yield, double p_c) {
    // Calculate relevant parameters for this analysis
    double* Vcmax = new double[4];
    double* photon_flux = new double[4];
    double* leaf_respiration = new double[4];
    double* w_l = new double[4];
    for (int i = 0; i < 4; i++) { 
        Vcmax[i] = calculate_Vcmax(temp[i], Vcmax15, Q10Vcmax, f_cold, f_hot, T_cold, T_hot);
        photon_flux[i] = calculate_photon_flux(absorbed_PAR[i], f_clump);
        leaf_respiration[i] = calculate_leaf_respiration(f_R, Vcmax[i]);
        w_l[i] = calculate_w_sat(temp[i], p_c); // In the future, may attempt to model changes in p_c based on atmospheric forcings
    }
    // Find the root to S193 for each of the potential limiting cases (RuBisCo, light, and CO2) for C4 plants.
    /// RuBisCo
    double* A_rubisco = new double[4];
    double* F_A = new double[4];
    double* F_B = new double[4];
    double* F_C = new double[4];
    double* F_D = new double[4];

    for (int i = 0; i < 4; i++) {
        F_A[i] = 0;
        F_B[i] = Vcmax[i];
        F_C[i] = 0;
        F_D[i] = 1;
    }

    calculate_A(A_rubisco, F_A, F_B, F_C, F_D,
                G_Wl, 
                G_Wlambda, G_Clambda, M, 0, 
                c_c, w_c, w_l, dw, leaf_respiration);

    /// CO2
    double* A_CO2 = new double[4];

    for (int i = 0; i < 4; i++) {
        F_A[i] = k_PEP*Vcmax[i];
        F_B[i] = 0;
        F_C[i] = 0;
        F_D[i] = 1;
    }

    calculate_A(A_CO2, F_A, F_B, F_C, F_D,
                G_Wl, 
                G_Wlambda, G_Clambda, M, 0, 
                c_c, w_c, w_l, dw, leaf_respiration);

    /// Light
    double* A_light = new double[4];

    for (int i = 0; i < 4; i++) {
        F_A[i] = 0;
        F_B[i] = quantum_yield*photon_flux[i];
        F_C[i] = 0;
        F_D[i] = 1;
    }

    calculate_A(A_light, F_A, F_B, F_C, F_D,
                G_Wl, 
                G_Wlambda, G_Clambda, M, 0, 
                c_c, w_c, w_l, dw, leaf_respiration);

    //printf("%f %f %f\n", A_rubisco, A_CO2, A_light);
    // Return the minimum estimate of A to reflect limiting conditions (Farquhar model)
    //printf("Rubisco %f CO2 %f Light %f\n", A_rubisco, A_CO2, A_light);
    for (int i = 0; i < 4; i++) {
        if ((A_rubisco < A_CO2) && (A_rubisco < A_light))
            A[i] = A_rubisco[i];
        else if (A_CO2 < A_light)
            A[i] = A_CO2[i];
        else 
            A[i] = A_light[i];
    }
}

double calculate_water_flux(int stomata, double LAI, double G_wlambda, double G_wl, double w_c, double temp, double p_c) {
    double w_l = calculate_w_sat(temp, p_c);
    double E = (G_wlambda*G_wl) / (G_wlambda + G_wl) * (w_c - w_l);
    
    return stomata*LAI*M_w*E;
}

double calculate_soil_cohort_enthalpy(double water_flux, double soil_temp) {
    return water_flux*q_l_water*(soil_temp - T_l0);
}

double calculate_cohort_CAS_enthalpy(double water_flux, double temp) {
    return water_flux*q_pv*(temp - T_v0);
}