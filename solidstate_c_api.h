#pragma once
// ============================================================================
//  solidstate_c_api.h  —  C API wrapper for SolidStateSystem
// ============================================================================

#ifdef __cplusplus
extern "C" {
#endif

// Lattice
double ss_lattice_spacing(int h, int k, int l, double a);
double ss_bragg_angle(double d, double wavelength, int n);
double ss_monatomic_dispersion(double kval, double a, double springK, double mass);

// Thermal
double ss_debye_specific_heat(double T, double thetaD);
double ss_einstein_specific_heat(double T, double thetaE);
double ss_thermal_conductivity(double Cv, double v, double mfp);

// Electronic
double ss_free_electron_energy(double k, double mStar);
double ss_fermi_energy(double n_density, double mStar);
double ss_fermi_dirac(double E, double EF, double T);
double ss_tight_binding_1d(double E0, double t, double k, double a);

// Semiconductor
double ss_intrinsic_carrier_conc(double Eg, double T, double mStar_e, double mStar_h);
double ss_built_in_voltage(double Na, double Nd, double ni, double T);
double ss_diode_current(double Is, double V, double T, double n_ideal);

// Magnetism
double ss_brillouin_function(double x, double J);
double ss_curie_weiss_temperature(double n, double J, double g, double lambda_mf);
double ss_ising_2d_tc(double J);

// Superconductivity
double ss_bcs_tc(double omegaD, double N_EF, double V_eff);
double ss_bcs_gap(double T, double Delta0, double Tc);
double ss_flux_quantum();
double ss_josephson_ac_freq(double V);

// Transport
double ss_landauer_conductance(double T_trans, int numChannels);
double ss_conductance_quantum();
double ss_hall_coefficient(double n_carrier, double charge);

#ifdef __cplusplus
}
#endif
