#include <cmath>
#include "parameters.h"
#include "numerics.h"

// Functions used for calculating photosynthesis of cohorts

double calculate_metabolic_rate(double temp, double x15, double Q10x) {
    return x15*std::pow(Q10x, (temp - 288.15)/10);
}

double calculate_michaelis_constants(double temp, double o2_mixing_ratio, 
                                     double C15, double Q10C, double O15, double Q10O) {
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

// The following functions are used for the CO2 mixing ratio solver, which is important for extracting A' and 
// ultimately calculating photosynthesis with respect to different limiting resources

double uptake(double F_A, double F_B, double F_C, double F_D, double co2_mixing_ratio) {
    return (F_A*co2_mixing_ratio - F_B) / (F_C*co2_mixing_ratio - F_D);
}

double F_1(double co2_mixing_ratio, double f_Gl, double f_Glambda, double G_Wl, double G_Wlambda, double M, double c_c,
           double F_A, double F_B, double F_C, double F_D) {
    double A = uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio);
    return ((f_Gl - f_Glambda*G_Wl/G_Wlambda)*A - G_Wl*(c_c - co2_mixing_ratio)) / (M*A);
}

double F_2(double co2_mixing_ratio, double f_Gl, double f_Glambda, double G_Wl, double G_Wlambda, double M, double gamma, double c_c,
           double F_A, double F_B, double F_C, double F_D) {
    double A = uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio);
    return (G_Wlambda*(c_c - gamma) - f_Glambda*A) / (G_Wlambda*(c_c - co2_mixing_ratio) + (f_Gl - f_Glambda)*A);
}

double F_3(double co2_mixing_ratio, double f_Gl, double f_Glambda, double G_Wl, 
           double G_Wlambda, double M, double gamma, double c_c, double w_c, double w_l, double dw, 
           double F_A, double F_B, double F_C, double F_D) {
    return 1 + (w_c - w_l) / (dw) * F_2(co2_mixing_ratio, f_Gl, f_Glambda, G_Wl, G_Wlambda, M, gamma, c_c, F_A, F_B, F_C, F_D);
}

double F(double co2_mixing_ratio, double f_Gl, double f_Glambda, double G_Wl, 
           double G_Wlambda, double M, double gamma, double c_c, double w_c, double w_l, double dw, 
           double F_A, double F_B, double F_C, double F_D) {
    return F_1(co2_mixing_ratio, f_Gl, f_Glambda, G_Wl, G_Wlambda, M, c_c,
           F_A, F_B, F_C, F_D) * 
           F_2(co2_mixing_ratio, f_Gl, f_Glambda, G_Wl, G_Wlambda, M, gamma, c_c,
               F_A, F_B, F_C, F_D) * 
           F_3(co2_mixing_ratio, f_Gl, f_Glambda, G_Wl, 
               G_Wlambda, M, gamma, c_c, w_c, w_l, dw, 
               F_A, F_B, F_C, F_D) - 1;
}

double calculate_A(double F_A, double F_B, double F_C, double F_D,
                   double f_Gl, double f_Glambda, double G_Wl, 
                   double G_Wlambda, double G_Clambda, double M, double gamma, 
                   double c_c, double w_c, double w_l, double dw, double leaf_respiration) {
    double c_lmin = -(F_D*M - F_B) / (F_C*M - F_A);

    double b = (G_Clambda*F_D + F_B - F_C*(G_Clambda*c_c + M)) / (G_Clambda*F_C);
    double c = (F_B - F_D*(G_Clambda*c_c + M)) / (G_Clambda*F_C);
    double c_lmax = (-b - std::sqrt(b*b - 4*c)) / 2;
    c_lmax = (c_lmax < c_c ? c_lmax : c_c); // Ensure singularity is within realm of plausible values

    // Create lambda function to fill in parameters of F
    auto F_spef = [&f_Gl, &f_Glambda, &G_Wl, &G_Wlambda, &M, &gamma, 
                   &c_c, &w_c, &w_l, &dw, &F_A, &F_B, &F_C, &F_D](double co2_mixing_ratio) 
                    -> double {return F(co2_mixing_ratio, f_Gl, f_Glambda, G_Wl, G_Wlambda, M, gamma,
                                        c_c, w_c, w_l, dw, F_A, F_B, F_C, F_D);};
    double co2_mixing_ratio = bisection(F_spef, c_lmin, c_lmax, 1e-6, 100);

    // If no roots in this range, assume stomata are closed and no photosynthesis is occurring.
    double A = (co2_mixing_ratio > 0 ? uptake(F_A, F_B, F_C, F_D, co2_mixing_ratio) : 0);
    return A - leaf_respiration;
}

// In reality, only noteworthy C3 pfts are grasses, so could probably just remove C3 functions

double co2_mixing_ratio_solver_C3(double f_Gl, double f_Glambda, double G_Wl, 
                                  double G_Wlambda, double G_Clambda, double M, double gamma, 
                                  double c_c, double w_c, double w_l, double dw,
                                  double temp, double Vcmax15, double Q10Vcmax, double k_PEP,
                                  double f_cold, double f_hot, double T_cold, double T_hot,
                                  double absorbed_PAR, double f_clump, double f_R, double quantum_yield,
                                  double o2_mixing_ratio, double C15, double Q10C,
                                  double O15, double Q10O) {
    // Calculate relevant parameters for this analysis
    double Vcmax = calculate_Vcmax(temp, Vcmax15, Q10Vcmax, f_cold, f_hot, T_cold, T_hot);
    double photon_flux = calculate_photon_flux(absorbed_PAR, f_clump);
    double leaf_respiration = calculate_leaf_respiration(f_R, Vcmax);
    double K_ME = calculate_michaelis_constants(temp, o2_mixing_ratio, 
                                                C15, Q10C, O15, Q10O);

    // Find the root to S193 for each of the potential limiting cases (RuBisCo, light, and CO2) for C3 plants.
    /// RuBisCo
    double A_rubisco = calculate_A(Vcmax, -Vcmax*gamma, 1, K_ME,
                                   f_Gl, f_Glambda, G_Wl, 
                                   G_Wlambda, G_Clambda, M, gamma, 
                                   c_c, w_c, w_l, dw, leaf_respiration);

    /// CO2, same as RuBisCo in this case

    /// Light
    double A_light   = calculate_A(quantum_yield*photon_flux, -quantum_yield*photon_flux*gamma, 1, 2*gamma,
                                   f_Gl, f_Glambda, G_Wl, 
                                   G_Wlambda, G_Clambda, M, gamma, 
                                   c_c, w_c, w_l, dw, leaf_respiration);


    // Return the minimum estimate of A to reflect limiting conditions (Farquhar model)
    if (A_rubisco < A_light)
        return A_rubisco;
    else 
        return A_light;
}

double co2_mixing_ratio_solver_C4(double f_Gl, double f_Glambda, double G_Wl, 
                                  double G_Wlambda, double G_Clambda, double M, double gamma, 
                                  double c_c, double w_c, double w_l, double dw,
                                  double temp, double Vcmax15, double Q10Vcmax, double k_PEP,
                                  double f_cold, double f_hot, double T_cold, double T_hot,
                                  double absorbed_PAR, double f_clump, double f_R, double quantum_yield) {
    // Calculate relevant parameters for this analysis
    double Vcmax = calculate_Vcmax(temp, Vcmax15, Q10Vcmax, f_cold, f_hot, T_cold, T_hot);
    double photon_flux = calculate_photon_flux(absorbed_PAR, f_clump);
    double leaf_respiration = calculate_leaf_respiration(f_R, Vcmax);
    
    // Find the root to S193 for each of the potential limiting cases (RuBisCo, light, and CO2) for C4 plants.
    /// RuBisCo
    double A_rubisco = calculate_A(0, Vcmax, 0, 1,
                                   f_Gl, f_Glambda, G_Wl, 
                                   G_Wlambda, G_Clambda, M, gamma, 
                                   c_c, w_c, w_l, dw, leaf_respiration);

    /// CO2
    double A_CO2     = calculate_A(k_PEP*Vcmax, 0, 0, 1,
                                   f_Gl, f_Glambda, G_Wl, 
                                   G_Wlambda, G_Clambda, M, gamma, 
                                   c_c, w_c, w_l, dw, leaf_respiration);

    /// Light
    double A_light   = calculate_A(0, quantum_yield*photon_flux, 0, 1,
                                   f_Gl, f_Glambda, G_Wl, 
                                   G_Wlambda, G_Clambda, M, gamma, 
                                   c_c, w_c, w_l, dw, leaf_respiration);


    // Return the minimum estimate of A to reflect limiting conditions (Farquhar model)
    if ((A_rubisco < A_CO2) && (A_rubisco < A_light))
        return A_rubisco;
    else if (A_CO2 < A_light)
        return A_CO2;
    else 
        return A_light;
}