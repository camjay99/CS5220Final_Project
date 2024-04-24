#include "parameters.h"
#include <cmath>

// Calculate relevant values for canopy air space

double calculate_grashof_number(double air_temp, double leaf_temp, double leaf_size) {
    
    double v_c = (1.33e-5)*(1 + 0007*(air_temp - 273.15));
    return (1 / air_temp)*g*leaf_size*leaf_size*leaf_size / (v_c*v_c) * std::abs(leaf_temp - air_temp);
}

double calculate_free_nusselt_number(double air_temp, double leaf_temp, double leaf_size) {
    double Gr = calculate_grashof_number(air_temp, leaf_temp, leaf_size);
    double laminar = 0.50*std::pow(Gr, 1/2);
    double turbulent = 0.13*std::pow(Gr, 1/3);
    return (laminar > turbulent ? laminar : turbulent);
}

double calculate_reynolds_number(double wind_speed, double leaf_size, double thermal_dif) {
    return wind_speed*leaf_size / thermal_dif;
}

double calculate_forced_nusselt_number(double wind_speed, double leaf_size, double thermal_dif) {
    double Re = calculate_reynolds_number(wind_speed, leaf_size, thermal_dif);
    double laminar = 0.60*std::pow(Re, 0.5);
    double turbulent = 0.032*std::pow(Re, 0.8);

    return (laminar > turbulent ? laminar : turbulent);
}

// Note wind speed is assumed to be a driver in this model, and not impacted by drag.
double calculate_G_Qlambda(double air_temp, double leaf_temp, double leaf_size, double wind_speed) {
    double thermal_dif = (1.89e-5)*(1 + 0007*(air_temp - 273.15));
    double free_Nu = calculate_free_nusselt_number(air_temp, leaf_temp, leaf_size);
    double forced_Nu = calculate_forced_nusselt_number(wind_speed, leaf_size, thermal_dif);

    double heat_conductance = thermal_dif * (free_Nu + forced_Nu) / leaf_size;
    return heat_conductance;
}

double calculate_G_Wlambda(double G_Qlambda, double air_density) {
    double water_conductance = 1.075*G_Qlambda;
    return air_density*water_conductance/M_d;
}

double calculate_Q_cohort_canopy(double leaf_area, double air_density, double G_Qlambda, double air_temp, double leaf_temp, double w_c) {
    double q_p_c = (1 - w_c)*q_pd + w_c*q_pv;
    double q_cohort_canopy = G_Qlambda*air_density*q_p_c*(leaf_temp - air_temp);
    return 2*leaf_area*q_cohort_canopy;
}