#ifndef __CANOPY_AIR_SPACE_H__
#define __CANOPY_AIR_SPACE_H__

__device__ double calculate_G_Qlambda(double air_temp, double leaf_temp, double leaf_size, double wind_speed);

__device__ double calculate_G_Wlambda(double G_Qlambda, double air_density);

__device__ double calculate_Q_cohort_canopy(double leaf_area, double air_density, double G_Qlambda, double air_temp, double leaf_temp, double w_c);

#endif