#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

// Define program constants

// For now, we will assume all individuals are the same PFT (plant functional type), so there are only one set of parameters
// per spectral band to consider. Future versions should be altered to effectively account for variability in coefficients.

// PAR parameters
#define scat_PAR                0.14   // Scattering coefficient for PAR (1 - scat_PAR is the absorbtion of PAR).
#define scat_ground_PAR         0.2   // Scattering coefficient for PAR on the ground (1 - scat_PAR is the absorbtion of PAR).
#define backscat_PAR_direct     0.09   // Proporation of scattered direct radiation that is reflected.
#define backscat_PAR_diffuse    0.09   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_direct_PAR     20000.0 // Amount of incoming direct PAR (in the future, will be driven by climatic details)
#define incoming_diffuse_PAR    0 // Amount of incoming diffuse PAR (in the future, will be driven by climatic details)

// NIR parameters
#define scat_NIR                0.825   // Scattering coefficient for NIR (1 - scat_NIR is the absorbtion of NIR).
#define scat_ground_NIR         0.2   // Scattering coefficient for NIR on the ground (1 - scat_NIR is the absorbtion of NIR).
#define backscat_NIR_direct     0.577   // Proporation of scattered direct radiation that is reflected.
#define backscat_NIR_diffuse    0.577   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_direct_NIR     500.0  // Amount of incoming direct NIR (in the future, will be driven by climatic details)
#define incoming_diffuse_NIR    0  // Amount of incoming diffuse NIR (in the future, will be driven by climatic details)

// TIR parameters
#define scat_TIR                0.03   // Scattering coefficient for TIR (1 - scat_TIR is the absorbtion of PAR).
#define scat_ground_TIR         0.2   // Scattering coefficient for TIR on the ground (1 - scat_TIR is the absorbtion of PAR).
#define backscat_TIR_direct     0.03   // Proporation of scattered direct radiation that is reflected.
#define backscat_TIR_diffuse    0.03   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_direct_TIR     0.0   // Amount of incoming direct TIR (in the future, will be driven by climatic details)
#define incoming_diffuse_TIR    0.0   // Amount of incoming diffuse TIR (in the future, will be driven by climatic details)
#define stef_boltz_constant     0.00000005670374419

// General parameters
#define inv_op_depth_direct     10.0  // Inverse of the optical depth of direct radiation per unit of effective plant area index (measure of how far light travels through canopy)
#define inv_op_depth_diffuse    10.0  // Inverse of the optical depth of diffuse radiation per unit of effective plant area index (measure of how far light travels through canopy)
#define total_PAI               5.0   // Total plant area index for the entire layer.


// Photosynthesis parameters
#define ein                     0.0000217 // Average photon specific energy in the PAR band
#define o2_mixing_ratio         0.209     // Regerence oxygen mixing ratio
#define f_Gl                    1.6       // Water:co2 diffusivity ratio
#define f_Glambda               1.4       // Water:co2 leaf-boundary-layer conductance ratio
#define k_PEP                   17949     // Initial slope for PEP carboxylase in C4 photosynthesis
#define Q10C                    2.1       // Temperature factor for Michaelis constant of carboxylation
#define Q10O                    1.2       // Temperature factor for Michaelis constant of oxygenation
#define Q10COratio              0.57      // Temperature factor for carboxylase:oxygenase ratio
#define C15                     214.2     // Michaelis consant for carboxylation at 15 degrees C
#define O15                     0.2725    // MIchaelis consant for oxygenation at 15 degrees C
#define COratio15               4561      // Carboxylse:oxygenase ratio at 15 degrees C
#define M_d                     0.02897   // Molar mass of dry air
#define M_w                     0.01802   // Molar mass of water
#define M_C                     0.01201   // Molar mass of carbon
#define f_Glambda               1.4       // Water:CO2 leaf-boundary-layer conductance ratio

// Enthalpy transfer equations
#define T_l0                    56.79
#define T_v0                    -1558.86
#define q_pv                    1859      // Specific heat of water vapor at constant pressure
#define q_pd                    1005      // Specific heat of dry air at constant pressure
#define q_l_water               4186      // Specific heat of liquid water at constant pressure
#define g                       9.807     // Gravity
#endif