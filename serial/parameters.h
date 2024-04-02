#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

// Define program constants

// For now, we will assume all individuals are the same PFT (plant functional type), so there are only one set of parameters
// per spectral band to consider. Future versions should be altered to effectively account for variability in coefficients.

// PAR parameters
#define scat_PAR                0.2   // Scattering coefficient for PAR (1 - scat_PAR is the absorbtion of PAR).
#define scat_ground_PAR         0.2   // Scattering coefficient for PAR on the ground (1 - scat_PAR is the absorbtion of PAR).
#define backscat_PAR_direct     0.1   // Proporation of scattered direct radiation that is reflected.
#define backscat_PAR_diffuse    0.1   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_direct_PAR     100.0 // Amount of incoming direct PAR (in the future, will be driven by climatic details)
#define incoming_diffuse_PAR    0 // Amount of incoming diffuse PAR (in the future, will be driven by climatic details)

// NIR parameters
#define scat_NIR                0.2   // Scattering coefficient for NIR (1 - scat_NIR is the absorbtion of NIR).
#define scat_ground_NIR         0.2   // Scattering coefficient for NIR on the ground (1 - scat_NIR is the absorbtion of NIR).
#define backscat_NIR_direct     0.1   // Proporation of scattered direct radiation that is reflected.
#define backscat_NIR_diffuse    0.1   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_direct_NIR     50.0  // Amount of incoming direct NIR (in the future, will be driven by climatic details)
#define incoming_diffuse_NIR    0  // Amount of incoming diffuse NIR (in the future, will be driven by climatic details)

// TIR parameters
#define scat_TIR                0.2   // Scattering coefficient for TIR (1 - scat_TIR is the absorbtion of PAR).
#define scat_ground_TIR         0.2   // Scattering coefficient for TIR on the ground (1 - scat_TIR is the absorbtion of PAR).
#define backscat_TIR_direct     0.1   // Proporation of scattered direct radiation that is reflected.
#define backscat_TIR_diffuse    0.1   // Proporation of scattered diffuse radiation that is reflected.
#define incoming_direct_TIR     0.0   // Amount of incoming direct TIR (in the future, will be driven by climatic details)
#define incoming_diffuse_TIR    0.0   // Amount of incoming diffuse TIR (in the future, will be driven by climatic details)
#define stef_boltz_constant     0.00000005670374419

// General parameters
#define inv_op_depth_direct     10.0  // Inverse of the optical depth of direct radiation per unit of effective plant area index (measure of how far light travels through canopy)
#define inv_op_depth_diffuse    10.0  // Inverse of the optical depth of diffuse radiation per unit of effective plant area index (measure of how far light travels through canopy)
#define total_PAI               5.0   // Total plant area index for the entire layer.
#define ein                     1
#endif