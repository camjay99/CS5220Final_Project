#ifndef __PHOTOSYNTHESIS_H__
#define __PHOTOSYNTHESIS_H__

double co2_mixing_ratio_solver_C3(double G_Wl, 
                                  double G_Wlambda, double G_Clambda, double M, 
                                  double c_c, double w_c, double dw,
                                  double temp, double Vcmax15, double Q10Vcmax,
                                  double f_cold, double f_hot, double T_cold, double T_hot,
                                  double absorbed_PAR, double f_clump, double f_R, double quantum_yield, double p_c);

double co2_mixing_ratio_solver_C4(double G_Wl, 
                                  double G_Wlambda, double G_Clambda, double M, 
                                  double c_c, double w_c, double dw,
                                  double temp, double Vcmax15, double Q10Vcmax,
                                  double f_cold, double f_hot, double T_cold, double T_hot,
                                  double absorbed_PAR, double f_clump, double f_R, double quantum_yield, double p_c);

double calculate_water_flux(int stomata, double LAI, double A);
double calculate_soil_cohort_enthalpy(double water_flux, double soil_temp);
double calculate_cohort_CAS_enthalpy(double water_flux, double temp);
#endif