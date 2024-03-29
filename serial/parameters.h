#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

// Define program constants

// For now, we will assume all individuals are the same PFT (plant functional type), so there are only one set of parameters
// per spectral band to consider. Future versions should be altered to effectively account for variability in coefficients.

// PAR parameters
#define scat_PAR                0.2   // Scattering coefficient for PAR (1 - scat_PAR is the absorbtion of PAR).
#define backscat_PAR_direct     0.1   // Proporation of scattered direct radiation that is reflected.
#define backscat_PAR_diffuse    0.1   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_PAR            100.0 // Amount of incoming PAR (in the future, will be driven by climatic details)

// NIR parameters
#define scat_NIR                0.2   // Scattering coefficient for PAR (1 - scat_PAR is the absorbtion of PAR).
#define backscat_NIR_direct     0.1   // Proporation of scattered direct radiation that is reflected.
#define backscat_NIR_diffuse    0.1   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_NIR            50.0  // Amount of incoming PAR (in the future, will be driven by climatic details)

// TIR parameters
#define scat_TIR                0.2   // Scattering coefficient for PAR (1 - scat_PAR is the absorbtion of PAR).
#define backscat_TIR_direct     0.1   // Proporation of scattered direct radiation that is reflected.
#define backscat_TIR_diffuse    0.1   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_TIR            0.0   // Amount of incoming PAR (in the future, will be driven by climatic details)


// General parameters
#define inv_op_depth_direct     10.0  // Inverse of the optical depth of direct radiation per unit of effective plant area index (measure of how far light travels through canopy)
#define inv_op_depth_diffuse    10.0  // Inverse of the optical depth of diffuse radiation per unit of effective plant area index (measure of how far light travels through canopy)
#define total_PAI               5.0   // Total plant area index for the entire layer.
#endif