// ============================================================================
//  solidstate_c_api.cpp  —  C API wrapper implementation
// ============================================================================

#include "solidstate_c_api.h"
#include "SolidStateSystem.h"

extern "C" {

double ss_lattice_spacing(int h, int k, int l, double a) {
    return SolidStateSystem::latticeSpacing(h, k, l, a);
}
double ss_bragg_angle(double d, double wavelength, int n) {
    return SolidStateSystem::braggAngle(d, wavelength, n);
}
double ss_monatomic_dispersion(double kval, double a, double springK, double mass) {
    return SolidStateSystem::monatomicDispersion(kval, a, springK, mass);
}
double ss_debye_specific_heat(double T, double thetaD) {
    return SolidStateSystem::debyeSpecificHeat(T, thetaD);
}
double ss_einstein_specific_heat(double T, double thetaE) {
    return SolidStateSystem::einsteinSpecificHeat(T, thetaE);
}
double ss_thermal_conductivity(double Cv, double v, double mfp) {
    return SolidStateSystem::kineticTheoryThermalConductivity(Cv, v, mfp);
}
double ss_free_electron_energy(double k, double mStar) {
    return SolidStateSystem::freeElectronEnergy(k, mStar);
}
double ss_fermi_energy(double n_density, double mStar) {
    return SolidStateSystem::fermiEnergy(n_density, mStar);
}
double ss_fermi_dirac(double E, double EF, double T) {
    return SolidStateSystem::fermiDiracDistribution(E, EF, T);
}
double ss_tight_binding_1d(double E0, double t, double k, double a) {
    return SolidStateSystem::tightBindingEnergy1D(E0, t, k, a);
}
double ss_intrinsic_carrier_conc(double Eg, double T, double mStar_e, double mStar_h) {
    return SolidStateSystem::intrinsicCarrierConcentration(Eg, T, mStar_e, mStar_h);
}
double ss_built_in_voltage(double Na, double Nd, double ni, double T) {
    return SolidStateSystem::builtInVoltage(Na, Nd, ni, T);
}
double ss_diode_current(double Is, double V, double T, double n_ideal) {
    return SolidStateSystem::diodeCurrent(Is, V, T, n_ideal);
}
double ss_brillouin_function(double x, double J) {
    return SolidStateSystem::brillouinFunction(x, J);
}
double ss_curie_weiss_temperature(double n, double J, double g, double lambda_mf) {
    return SolidStateSystem::curieWeissTemperature(n, J, g, lambda_mf);
}
double ss_ising_2d_tc(double J) {
    return SolidStateSystem::ising2DCriticalTemperature(J);
}
double ss_bcs_tc(double omegaD, double N_EF, double V_eff) {
    return SolidStateSystem::bcsCriticalTemperature(omegaD, N_EF, V_eff);
}
double ss_bcs_gap(double T, double Delta0, double Tc) {
    return SolidStateSystem::bcsGap(T, Delta0, Tc);
}
double ss_flux_quantum() {
    return SolidStateSystem::fluxQuantum();
}
double ss_josephson_ac_freq(double V) {
    return SolidStateSystem::josephsonACFrequency(V);
}
double ss_landauer_conductance(double T_trans, int numChannels) {
    return SolidStateSystem::landauerConductance(T_trans, numChannels);
}
double ss_conductance_quantum() {
    return SolidStateSystem::conductanceQuantum();
}
double ss_hall_coefficient(double n_carrier, double charge) {
    return SolidStateSystem::hallCoefficient(n_carrier, charge);
}

} // extern "C"
