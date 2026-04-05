#pragma once
// ============================================================================
//  SolidStateConstants.h  —  Physical constants for condensed-matter physics
// ============================================================================

#include <cmath>

namespace SSConst {

// ── Fundamental ─────────────────────────────────────────────────────────────
constexpr double hbar       = 1.054571817e-34;   // J·s
constexpr double h_planck   = 6.62607015e-34;    // J·s
constexpr double e_charge   = 1.602176634e-19;   // C
constexpr double me         = 9.1093837015e-31;  // kg  (free electron mass)
constexpr double kB         = 1.380649e-23;      // J/K
constexpr double c_light    = 2.99792458e8;      // m/s
constexpr double epsilon0   = 8.8541878128e-12;  // F/m
constexpr double mu0        = 1.25663706212e-6;  // N/A^2
constexpr double NA         = 6.02214076e23;     // 1/mol

// ── Derived / Useful ────────────────────────────────────────────────────────
constexpr double eV         = 1.602176634e-19;   // J per eV
constexpr double Angstrom   = 1.0e-10;           // m
constexpr double a0_bohr    = 5.29177210903e-11;  // Bohr radius  (m)
constexpr double Ry         = 13.605693122994;    // Rydberg energy (eV)
constexpr double muB        = 9.2740100783e-24;   // Bohr magneton (J/T)
constexpr double muN        = 5.0507837461e-27;   // Nuclear magneton (J/T)
constexpr double alpha_fs   = 7.2973525693e-3;    // Fine-structure constant

// ── Lattice reference values ────────────────────────────────────────────────
constexpr double a_Si       = 5.431e-10;         // Si lattice constant (m)
constexpr double a_GaAs     = 5.653e-10;         // GaAs lattice constant (m)
constexpr double a_Cu       = 3.615e-10;         // Cu lattice constant (m)

// ── Semiconductor effective masses (in units of me) ─────────────────────────
constexpr double mStar_Si_e = 0.26;              // Si electron eff mass
constexpr double mStar_Si_h = 0.36;              // Si hole eff mass
constexpr double mStar_GaAs = 0.067;             // GaAs electron eff mass

// ── Math ────────────────────────────────────────────────────────────────────
constexpr double PI         = 3.14159265358979323846;
constexpr double TWO_PI     = 2.0 * PI;

} // namespace SSConst
