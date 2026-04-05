// ============================================================================
//  SolidStateSystem.cpp  —  Core physics implementations
// ============================================================================

#include "SolidStateSystem.h"
#include "SolidStateConstants.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <cassert>

using namespace std;
using namespace SSConst;

// ============================================================================
//  Constructor & Utility
// ============================================================================

SolidStateSystem::SolidStateSystem(const string& materialName,
                                   double latticeConstant,
                                   double effectiveMass)
    : m_materialName(materialName), m_a(latticeConstant), m_mStar(effectiveMass) {}

void SolidStateSystem::setMaterial(const string& name, double a, double mStar) {
    m_materialName = name;
    m_a = a;
    m_mStar = mStar;
}

// ============================================================================
//  1. Bravais Lattices
// ============================================================================

BravaisLattice SolidStateSystem::getBravaisLattice(int type, double a, double b, double c,
                                                    double alpha_a, double beta_a, double gamma_a) {
    BravaisLattice lat{};
    if (b == 0) b = a;
    if (c == 0) c = a;
    switch (type) {
        case 0: // Simple cubic
            lat.a1 = {a, 0, 0}; lat.a2 = {0, a, 0}; lat.a3 = {0, 0, a};
            lat.name = "Simple Cubic";
            break;
        case 1: // BCC
            lat.a1 = {a/2, a/2, -a/2}; lat.a2 = {-a/2, a/2, a/2}; lat.a3 = {a/2, -a/2, a/2};
            lat.name = "BCC";
            break;
        case 2: // FCC
            lat.a1 = {0, a/2, a/2}; lat.a2 = {a/2, 0, a/2}; lat.a3 = {a/2, a/2, 0};
            lat.name = "FCC";
            break;
        case 3: // Hexagonal
            lat.a1 = {a, 0, 0}; lat.a2 = {a/2, a*sqrt(3.0)/2.0, 0}; lat.a3 = {0, 0, c};
            lat.name = "Hexagonal";
            break;
        default:
            lat.a1 = {a, 0, 0}; lat.a2 = {0, b, 0}; lat.a3 = {0, 0, c};
            lat.name = "Orthorhombic";
    }
    return lat;
}

LatticeVector SolidStateSystem::computeReciprocalVector(const LatticeVector& a1,
                                                         const LatticeVector& a2,
                                                         const LatticeVector& a3, int which) {
    // b_i = 2*pi * (a_j x a_k) / (a_i . (a_j x a_k))
    auto cross = [](const LatticeVector& u, const LatticeVector& v) -> LatticeVector {
        return {u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x};
    };
    auto dot = [](const LatticeVector& u, const LatticeVector& v) -> double {
        return u.x*v.x + u.y*v.y + u.z*v.z;
    };
    LatticeVector c23 = cross(a2, a3);
    double vol = dot(a1, c23);
    LatticeVector num;
    if (which == 0) num = cross(a2, a3);
    else if (which == 1) num = cross(a3, a1);
    else num = cross(a1, a2);
    double f = TWO_PI / vol;
    return {num.x * f, num.y * f, num.z * f};
}

double SolidStateSystem::latticeSpacing(int h, int k, int l, double a) {
    return a / sqrt(static_cast<double>(h*h + k*k + l*l));
}

void SolidStateSystem::exportBravaisLatticeCSV(const string& filename, int type, double a) {
    auto lat = getBravaisLattice(type, a);
    ofstream out(filename);
    out << "vector,x,y,z\n";
    out << "a1," << lat.a1.x << "," << lat.a1.y << "," << lat.a1.z << "\n";
    out << "a2," << lat.a2.x << "," << lat.a2.y << "," << lat.a2.z << "\n";
    out << "a3," << lat.a3.x << "," << lat.a3.y << "," << lat.a3.z << "\n";
}

// ============================================================================
//  2. Reciprocal Lattice
// ============================================================================

LatticeVector SolidStateSystem::reciprocalLatticePoint(const LatticeVector& b1,
                                                        const LatticeVector& b2,
                                                        const LatticeVector& b3,
                                                        int n1, int n2, int n3) {
    return {n1*b1.x + n2*b2.x + n3*b3.x,
            n1*b1.y + n2*b2.y + n3*b3.y,
            n1*b1.z + n2*b2.z + n3*b3.z};
}

double SolidStateSystem::brillouinZoneBoundary1D(double a) {
    return PI / a;
}

void SolidStateSystem::exportReciprocalLatticeCSV(const string& filename, double a, int maxIndex) {
    ofstream out(filename);
    out << "n1,n2,n3,Gx,Gy,Gz\n";
    double b = TWO_PI / a;
    for (int i = -maxIndex; i <= maxIndex; ++i)
        for (int j = -maxIndex; j <= maxIndex; ++j)
            for (int k = -maxIndex; k <= maxIndex; ++k)
                out << i << "," << j << "," << k << ","
                    << i*b << "," << j*b << "," << k*b << "\n";
}

// ============================================================================
//  3. X-Ray Diffraction
// ============================================================================

double SolidStateSystem::braggAngle(double d, double wavelength, int n) {
    double sinTheta = n * wavelength / (2.0 * d);
    if (sinTheta > 1.0 || sinTheta < -1.0) return -1.0;
    return asin(sinTheta);
}

double SolidStateSystem::structureFactor(int h, int k, int l, const string& latticeType) {
    if (latticeType == "BCC") {
        return ((h + k + l) % 2 == 0) ? 2.0 : 0.0;
    } else if (latticeType == "FCC") {
        bool allEven = (h%2==0 && k%2==0 && l%2==0);
        bool allOdd  = (h%2!=0 && k%2!=0 && l%2!=0);
        return (allEven || allOdd) ? 4.0 : 0.0;
    }
    return 1.0;  // Simple cubic
}

vector<DiffractionPeak> SolidStateSystem::computeDiffractionPattern(double a, double wavelength,
                                                                      const string& latticeType, int maxHKL) {
    vector<DiffractionPeak> peaks;
    for (int h = 0; h <= maxHKL; ++h)
        for (int k = 0; k <= maxHKL; ++k)
            for (int l = 0; l <= maxHKL; ++l) {
                if (h == 0 && k == 0 && l == 0) continue;
                double sf = structureFactor(h, k, l, latticeType);
                if (sf < 1e-12) continue;
                double d = latticeSpacing(h, k, l, a);
                double theta = braggAngle(d, wavelength);
                if (theta < 0) continue;
                peaks.push_back({h, k, l, d, 2.0*theta, sf*sf});
            }
    return peaks;
}

double SolidStateSystem::scatteringVector(double theta, double wavelength) {
    return 4.0 * PI * sin(theta) / wavelength;
}

double SolidStateSystem::debyeWallerFactor(double B, double theta, double wavelength) {
    double q = scatteringVector(theta, wavelength);
    return exp(-B * q * q / (16.0 * PI * PI));
}

void SolidStateSystem::exportDiffractionCSV(const string& filename, double a, double wavelength,
                                              const string& latticeType, int maxHKL) {
    auto peaks = computeDiffractionPattern(a, wavelength, latticeType, maxHKL);
    ofstream out(filename);
    out << "h,k,l,d_spacing,twoTheta_rad,intensity\n";
    for (auto& p : peaks)
        out << p.h << "," << p.k << "," << p.l << ","
            << p.d_spacing << "," << p.twoTheta << "," << p.intensity << "\n";
}

// ============================================================================
//  4. Lattice Vibrations — 1D Monatomic
// ============================================================================

double SolidStateSystem::monatomicDispersion(double k, double a, double springK, double mass) {
    return 2.0 * sqrt(springK / mass) * abs(sin(k * a / 2.0));
}

double SolidStateSystem::monatomicGroupVelocity(double k, double a, double springK, double mass) {
    double omegaMax = 2.0 * sqrt(springK / mass);
    return (omegaMax * a / 2.0) * abs(cos(k * a / 2.0));
}

double SolidStateSystem::monatomicDOS1D(double omega, double omegaMax) {
    if (omega <= 0 || omega >= omegaMax) return 0.0;
    return 2.0 / (PI * sqrt(omegaMax * omegaMax - omega * omega));
}

double SolidStateSystem::debyeCutoffFrequency(double springK, double mass) {
    return 2.0 * sqrt(springK / mass);
}

void SolidStateSystem::exportMonatomicDispersionCSV(const string& filename, double a,
                                                      double springK, double mass, int numPoints) {
    ofstream out(filename);
    out << "k,omega,v_group\n";
    double kMax = PI / a;
    for (int i = 0; i <= numPoints; ++i) {
        double k = -kMax + 2.0 * kMax * i / numPoints;
        out << k << "," << monatomicDispersion(k, a, springK, mass) << ","
            << monatomicGroupVelocity(k, a, springK, mass) << "\n";
    }
}

// ============================================================================
//  5. Lattice Vibrations — 1D Diatomic
// ============================================================================

pair<double, double> SolidStateSystem::diatomicDispersion(double k, double a,
                                                            double springK, double m1, double m2) {
    double mu_inv = 1.0/m1 + 1.0/m2;
    double det = mu_inv * mu_inv - 4.0 * sin(k*a/2.0) * sin(k*a/2.0) / (m1*m2);
    if (det < 0) det = 0;
    double avg = springK * mu_inv;
    double diff = springK * sqrt(det);
    double acoustic = sqrt(avg - diff);
    double optical  = sqrt(avg + diff);
    return {acoustic, optical};
}

double SolidStateSystem::diatomicBandGap(double springK, double m1, double m2) {
    double omega_ac_max = sqrt(2.0 * springK / max(m1, m2));
    double omega_op_min = sqrt(2.0 * springK / min(m1, m2));
    return omega_op_min - omega_ac_max;
}

void SolidStateSystem::exportDiatomicDispersionCSV(const string& filename, double a,
                                                     double springK, double m1, double m2, int numPoints) {
    ofstream out(filename);
    out << "k,omega_acoustic,omega_optical\n";
    double kMax = PI / a;
    for (int i = 0; i <= numPoints; ++i) {
        double k = -kMax + 2.0 * kMax * i / numPoints;
        auto [ac, op] = diatomicDispersion(k, a, springK, m1, m2);
        out << k << "," << ac << "," << op << "\n";
    }
}

// ============================================================================
//  6. Debye Model
// ============================================================================

double SolidStateSystem::debyeFunction(double x, int n) {
    // D_n(x) = (n/x^n) * integral_0^x t^n / (e^t - 1) dt
    if (x < 1e-12) return 1.0;
    int numSteps = 1000;
    double dt = x / numSteps;
    double sum = 0.0;
    for (int i = 1; i <= numSteps; ++i) {
        double t = i * dt;
        double et = exp(t);
        if (et > 1.0 + 1e-15)
            sum += pow(t, n) / (et - 1.0);
    }
    sum *= dt;
    return (static_cast<double>(n) / pow(x, n)) * sum;
}

double SolidStateSystem::debyeSpecificHeat(double T, double thetaD) {
    if (T < 1e-12) return 0.0;
    double x = thetaD / T;
    return 3.0 * kB * NA * debyeFunction(x, 3) * 3.0; // 3NkB * D3
}

double SolidStateSystem::debyeEnergy(double T, double thetaD, int N) {
    if (T < 1e-12) return 0.0;
    double x = thetaD / T;
    return 3.0 * N * kB * T * debyeFunction(x, 3);
}

double SolidStateSystem::debyeEntropy(double T, double thetaD, int N) {
    if (T < 1e-12) return 0.0;
    double x = thetaD / T;
    // S = Nk[4D3(x) - 3 ln(1 - e^-x)]
    return N * kB * (4.0 * debyeFunction(x, 3) - 3.0 * log(1.0 - exp(-x)));
}

void SolidStateSystem::exportDebyeModelCSV(const string& filename, double thetaD, int N,
                                             double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,Cv,E,S\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        out << T << "," << debyeSpecificHeat(T, thetaD) << ","
            << debyeEnergy(T, thetaD, N) << "," << debyeEntropy(T, thetaD, N) << "\n";
    }
}

// ============================================================================
//  7. Einstein Model
// ============================================================================

double SolidStateSystem::einsteinSpecificHeat(double T, double thetaE) {
    if (T < 1e-12) return 0.0;
    double x = thetaE / T;
    double ex = exp(x);
    return 3.0 * kB * NA * (x * x * ex) / ((ex - 1.0) * (ex - 1.0));
}

double SolidStateSystem::einsteinEnergy(double T, double thetaE, int N) {
    if (T < 1e-12) return 0.0;
    double x = thetaE / T;
    return 3.0 * N * kB * thetaE / (exp(x) - 1.0);
}

double SolidStateSystem::einsteinEntropy(double T, double thetaE, int N) {
    if (T < 1e-12) return 0.0;
    double x = thetaE / T;
    double ex = exp(x);
    return 3.0 * N * kB * (x / (ex - 1.0) - log(1.0 - 1.0/ex));
}

void SolidStateSystem::exportEinsteinModelCSV(const string& filename, double thetaE, int N,
                                                double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,Cv,E,S\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        out << T << "," << einsteinSpecificHeat(T, thetaE) << ","
            << einsteinEnergy(T, thetaE, N) << "," << einsteinEntropy(T, thetaE, N) << "\n";
    }
}

// ============================================================================
//  8. Phonon Density of States
// ============================================================================

double SolidStateSystem::debyeDOS3D(double omega, double omegaD) {
    if (omega < 0 || omega > omegaD) return 0.0;
    return 9.0 * omega * omega / (omegaD * omegaD * omegaD);
}

double SolidStateSystem::vanHoveSingularity1D(double omega, double omegaMax) {
    return monatomicDOS1D(omega, omegaMax);
}

vector<double> SolidStateSystem::computePhononDOS(const vector<double>& frequencies,
                                                    double broadening, int numBins) {
    if (frequencies.empty()) return {};
    double fMin = *min_element(frequencies.begin(), frequencies.end());
    double fMax = *max_element(frequencies.begin(), frequencies.end());
    double df = (fMax - fMin) / numBins;
    vector<double> dos(numBins, 0.0);
    for (double f : frequencies) {
        for (int i = 0; i < numBins; ++i) {
            double fc = fMin + (i + 0.5) * df;
            double x = (f - fc) / broadening;
            dos[i] += exp(-0.5 * x * x) / (broadening * sqrt(TWO_PI));
        }
    }
    return dos;
}

void SolidStateSystem::exportPhononDOSCSV(const string& filename, double omegaD, int numPoints) {
    ofstream out(filename);
    out << "omega,DOS_Debye\n";
    for (int i = 0; i <= numPoints; ++i) {
        double omega = omegaD * i / numPoints;
        out << omega << "," << debyeDOS3D(omega, omegaD) << "\n";
    }
}

// ============================================================================
//  9. Thermal Expansion
// ============================================================================

double SolidStateSystem::grueneisenParameter(double bulkModulus, double Cv, double volume, double alpha_t) {
    if (Cv < 1e-30) return 0.0;
    return alpha_t * bulkModulus * volume / Cv;
}

double SolidStateSystem::thermalExpansionCoeff(double gamma, double Cv, double bulkModulus, double volume) {
    if (bulkModulus < 1e-30 || volume < 1e-30) return 0.0;
    return gamma * Cv / (bulkModulus * volume);
}

double SolidStateSystem::anharmonicPotential(double x, double k2, double k3, double k4) {
    return 0.5 * k2 * x * x + k3 * x * x * x / 6.0 + k4 * x * x * x * x / 24.0;
}

void SolidStateSystem::exportThermalExpansionCSV(const string& filename, double thetaD, double gamma,
                                                   double V0, double B, double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,alpha,Cv\n";
    int N = static_cast<int>(NA);
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        double Cv = debyeSpecificHeat(T, thetaD);
        double alpha_t = thermalExpansionCoeff(gamma, Cv, B, V0);
        out << T << "," << alpha_t << "," << Cv << "\n";
    }
}

// ============================================================================
//  10. Lattice Thermal Conductivity
// ============================================================================

double SolidStateSystem::kineticTheoryThermalConductivity(double Cv, double v, double mfp) {
    return Cv * v * mfp / 3.0;
}

double SolidStateSystem::umklappScatteringRate(double T, double thetaD, double omega) {
    return (omega * omega * T) / (thetaD * thetaD) * exp(-thetaD / (3.0 * T));
}

double SolidStateSystem::phononMeanFreePath(double T, double thetaD) {
    if (T < 1e-12) return 1e-3;
    return 1e-8 * (thetaD / T) * exp(thetaD / (3.0 * T));
}

void SolidStateSystem::exportThermalConductivityCSV(const string& filename, double thetaD,
                                                      double v_sound, double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,kappa,mfp\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        double Cv = debyeSpecificHeat(T, thetaD);
        double mfp = phononMeanFreePath(T, thetaD);
        double kappa = kineticTheoryThermalConductivity(Cv, v_sound, mfp);
        out << T << "," << kappa << "," << mfp << "\n";
    }
}

// ============================================================================
//  11. Free Electron Model
// ============================================================================

double SolidStateSystem::freeElectronEnergy(double k, double mStar) {
    return hbar * hbar * k * k / (2.0 * mStar);
}

double SolidStateSystem::fermiEnergy(double n_density, double mStar) {
    return (hbar * hbar / (2.0 * mStar)) * pow(3.0 * PI * PI * n_density, 2.0/3.0);
}

double SolidStateSystem::fermiWavevector(double EF, double mStar) {
    return sqrt(2.0 * mStar * EF) / hbar;
}

double SolidStateSystem::fermiVelocity(double EF, double mStar) {
    return sqrt(2.0 * EF / mStar);
}

double SolidStateSystem::fermiTemperature(double EF) {
    return EF / kB;
}

double SolidStateSystem::freeElectronDOS3D(double E, double mStar) {
    if (E < 0) return 0.0;
    return (1.0 / (2.0 * PI * PI)) * pow(2.0 * mStar / (hbar * hbar), 1.5) * sqrt(E);
}

double SolidStateSystem::freeElectronDOS2D(double mStar) {
    return mStar / (PI * hbar * hbar);
}

double SolidStateSystem::freeElectronDOS1D(double E, double mStar) {
    if (E <= 0) return 0.0;
    return (1.0 / PI) * sqrt(mStar / (2.0 * hbar * hbar * E));
}

double SolidStateSystem::fermiDiracDistribution(double E, double EF, double T) {
    if (T < 1e-12) return (E < EF) ? 1.0 : 0.0;
    double x = (E - EF) / (kB * T);
    if (x > 500) return 0.0;
    if (x < -500) return 1.0;
    return 1.0 / (exp(x) + 1.0);
}

double SolidStateSystem::electronicSpecificHeat(double T, double EF, double n_density) {
    if (T < 1e-12 || EF < 1e-30) return 0.0;
    return (PI * PI / 3.0) * n_density * kB * kB * T / EF;
}

void SolidStateSystem::exportFreeElectronCSV(const string& filename, double n_density, double mStar,
                                               double Emax, int numPoints) {
    ofstream out(filename);
    double EF = fermiEnergy(n_density, mStar);
    out << "E,DOS,f_FD_300K,f_FD_1000K\n";
    for (int i = 0; i <= numPoints; ++i) {
        double E = Emax * i / numPoints;
        out << E << "," << freeElectronDOS3D(E, mStar) << ","
            << fermiDiracDistribution(E, EF, 300.0) << ","
            << fermiDiracDistribution(E, EF, 1000.0) << "\n";
    }
}

// ============================================================================
//  12. Nearly Free Electron Model
// ============================================================================

pair<double, double> SolidStateSystem::NFEBandEdges(double k, double G, double V_G, double mStar) {
    double E0 = hbar * hbar * k * k / (2.0 * mStar);
    double E1 = hbar * hbar * (k - G) * (k - G) / (2.0 * mStar);
    double avg = (E0 + E1) / 2.0;
    double diff = sqrt((E0 - E1) * (E0 - E1) / 4.0 + V_G * V_G);
    return {avg - diff, avg + diff};
}

double SolidStateSystem::NFEBandGap(double V_G) {
    return 2.0 * abs(V_G);
}

BandResult SolidStateSystem::computeNFEBands(double a, double V_G, double mStar, int numBands, int numPoints) {
    BandResult result;
    double kMax = numBands * PI / a;
    result.k_points.resize(numPoints);
    result.bands.resize(2);
    result.bands[0].resize(numPoints);
    result.bands[1].resize(numPoints);
    double G = TWO_PI / a;
    for (int i = 0; i < numPoints; ++i) {
        double k = -kMax + 2.0 * kMax * i / (numPoints - 1);
        result.k_points[i] = k;
        auto [lower, upper] = NFEBandEdges(k, G, V_G * eV, mStar);
        result.bands[0][i] = lower / eV;
        result.bands[1][i] = upper / eV;
    }
    return result;
}

void SolidStateSystem::exportNFEBandsCSV(const string& filename, double a, double V_G,
                                           double mStar, int numBands, int numPoints) {
    auto result = computeNFEBands(a, V_G, mStar, numBands, numPoints);
    ofstream out(filename);
    out << "k,E_lower_eV,E_upper_eV\n";
    for (int i = 0; i < numPoints; ++i)
        out << result.k_points[i] << "," << result.bands[0][i] << "," << result.bands[1][i] << "\n";
}

// ============================================================================
//  13. Tight-Binding (Extended)
// ============================================================================

double SolidStateSystem::tightBindingEnergy1D(double E0, double t, double k, double a) {
    return E0 - 2.0 * t * cos(k * a);
}

double SolidStateSystem::tightBindingEnergy2D(double E0, double t, double kx, double ky, double a) {
    return E0 - 2.0 * t * (cos(kx * a) + cos(ky * a));
}

double SolidStateSystem::tightBindingEffectiveMass(double t, double a, double mStar) {
    (void)mStar;
    return hbar * hbar / (2.0 * t * a * a);
}

double SolidStateSystem::tightBindingDOS1D(double E, double E0, double t) {
    double x = (E - E0) / (2.0 * t);
    if (abs(x) >= 1.0) return 0.0;
    return 1.0 / (PI * 2.0 * t * sqrt(1.0 - x * x));
}

double SolidStateSystem::tightBindingBandwidth(double t) {
    return 4.0 * abs(t);
}

BandResult SolidStateSystem::computeTBBands1D(double E0, double t, double a, int numPoints) {
    BandResult result;
    result.k_points.resize(numPoints);
    result.bands.resize(1);
    result.bands[0].resize(numPoints);
    double kMax = PI / a;
    for (int i = 0; i < numPoints; ++i) {
        double k = -kMax + 2.0 * kMax * i / (numPoints - 1);
        result.k_points[i] = k;
        result.bands[0][i] = tightBindingEnergy1D(E0, t, k, a);
    }
    return result;
}

BandResult SolidStateSystem::computeTBBands2D(double E0, double t, double a, int numKPoints) {
    // High-symmetry path: Gamma -> X -> M -> Gamma
    BandResult result;
    int seg = numKPoints / 3;
    result.bands.resize(1);
    for (int i = 0; i < seg; ++i) { // Gamma to X
        double kx = PI / a * i / seg;
        result.k_points.push_back(static_cast<double>(i));
        result.bands[0].push_back(tightBindingEnergy2D(E0, t, kx, 0, a));
    }
    for (int i = 0; i < seg; ++i) { // X to M
        double ky = PI / a * i / seg;
        result.k_points.push_back(static_cast<double>(seg + i));
        result.bands[0].push_back(tightBindingEnergy2D(E0, t, PI/a, ky, a));
    }
    for (int i = 0; i <= seg; ++i) { // M to Gamma
        double k = PI / a * (1.0 - static_cast<double>(i) / seg);
        result.k_points.push_back(static_cast<double>(2*seg + i));
        result.bands[0].push_back(tightBindingEnergy2D(E0, t, k, k, a));
    }
    return result;
}

void SolidStateSystem::exportTBBandsCSV(const string& filename, double E0, double t, double a, int numPoints) {
    auto result = computeTBBands1D(E0, t, a, numPoints);
    ofstream out(filename);
    out << "k,E\n";
    for (int i = 0; i < numPoints; ++i)
        out << result.k_points[i] << "," << result.bands[0][i] << "\n";
}

// ============================================================================
//  14. Kronig-Penney (Extended)
// ============================================================================

double SolidStateSystem::kronigPenneyDispersion(double E, double V0, double a, double b, double mStar) {
    double d = a + b;
    if (E < V0) {
        double alpha = sqrt(2.0 * mStar * E) / hbar;
        double beta  = sqrt(2.0 * mStar * (V0 - E)) / hbar;
        if (alpha < 1e-30 || beta < 1e-30) return 2.0;
        return cos(alpha * a) * cosh(beta * b)
             - ((beta*beta - alpha*alpha) / (2.0*alpha*beta)) * sin(alpha*a) * sinh(beta*b);
    } else {
        double alpha = sqrt(2.0 * mStar * E) / hbar;
        double beta  = sqrt(2.0 * mStar * (E - V0)) / hbar;
        return cos(alpha * a) * cos(beta * b)
             - ((beta*beta + alpha*alpha) / (2.0*alpha*beta)) * sin(alpha*a) * sin(beta*b);
    }
}

double SolidStateSystem::kronigPenneyDeltaDispersion(double E, double P, double a) {
    double k = sqrt(2.0 * me * E) / hbar;
    if (k < 1e-30) return 2.0;
    return cos(k * a) + P * sin(k * a) / (k * a);
}

BandResult SolidStateSystem::computeKPBands(double V0, double a, double b, double mStar,
                                              int numBands, int numPoints) {
    // Scan energy to find allowed bands
    BandResult result;
    result.bands.resize(numBands);
    double Emax = V0 * (numBands + 1);
    int totalScan = numPoints * 100;
    vector<double> allowedE;
    for (int i = 0; i < totalScan; ++i) {
        double E = Emax * (i + 0.5) / totalScan;
        double val = kronigPenneyDispersion(E, V0, a, b, mStar);
        if (abs(val) <= 1.0) allowedE.push_back(E);
    }
    // Simple band assignment (placeholder for full implementation)
    for (int b_idx = 0; b_idx < numBands && !allowedE.empty(); ++b_idx) {
        result.bands[b_idx] = allowedE;
    }
    return result;
}

void SolidStateSystem::exportKPBandsCSV(const string& filename, double V0, double a, double b,
                                          double mStar, int numBands, int numPoints) {
    ofstream out(filename);
    out << "E,cos_ka\n";
    double Emax = V0 * (numBands + 1);
    for (int i = 0; i <= numPoints; ++i) {
        double E = Emax * i / numPoints;
        double val = kronigPenneyDispersion(E, V0, a, b, mStar);
        out << E << "," << val << "\n";
    }
}

// ============================================================================
//  15. APW / OPW Methods
// ============================================================================

double SolidStateSystem::muffinTinPotential(double r, double rMT, double V0) {
    return (r < rMT) ? -V0 : 0.0;
}

double SolidStateSystem::apwPhaseShift(int l, double E, double rMT, double V0) {
    // Simplified: phase shift for muffin-tin
    double k_in = sqrt(2.0 * me * (E + V0 * eV)) / hbar;
    double k_out = sqrt(2.0 * me * E) / hbar;
    if (k_in < 1e-30 || k_out < 1e-30) return 0.0;
    return atan2(k_out * cos(k_in * rMT), k_in * cos(k_out * rMT));
}

BandResult SolidStateSystem::computeAPWBands(double a, double rMT, double V0, int lMax, int numPoints) {
    (void)lMax;
    BandResult result;
    result.k_points.resize(numPoints);
    result.bands.resize(1);
    result.bands[0].resize(numPoints);
    double kMax = PI / a;
    for (int i = 0; i < numPoints; ++i) {
        double k = kMax * i / (numPoints - 1);
        result.k_points[i] = k;
        result.bands[0][i] = hbar*hbar*k*k / (2.0*me*eV) - V0; // simplified
    }
    return result;
}

void SolidStateSystem::exportAPWBandsCSV(const string& filename, double a, double rMT,
                                           double V0, int lMax, int numPoints) {
    auto result = computeAPWBands(a, rMT, V0, lMax, numPoints);
    ofstream out(filename);
    out << "k,E_eV\n";
    for (int i = 0; i < numPoints; ++i)
        out << result.k_points[i] << "," << result.bands[0][i] << "\n";
}

// ============================================================================
//  16. DOS & Fermi Surface
// ============================================================================

vector<double> SolidStateSystem::computeDOSFromBands(const BandResult& bands, double broadening, int numBins) {
    if (bands.bands.empty()) return {};
    double Emin = 1e30, Emax = -1e30;
    for (auto& b : bands.bands)
        for (double e : b) { Emin = min(Emin, e); Emax = max(Emax, e); }
    double dE = (Emax - Emin) / numBins;
    vector<double> dos(numBins, 0.0);
    for (auto& b : bands.bands)
        for (double e : b)
            for (int i = 0; i < numBins; ++i) {
                double Ec = Emin + (i + 0.5) * dE;
                double x = (e - Ec) / broadening;
                dos[i] += exp(-0.5*x*x) / (broadening*sqrt(TWO_PI));
            }
    return dos;
}

vector<FermiSurfacePoint> SolidStateSystem::computeFermiSurface2D(
    function<double(double,double)> Ek, double EF, int numPoints) {
    vector<FermiSurfacePoint> pts;
    double kMax = 2e10;
    for (int i = 0; i < numPoints; ++i) {
        double kx = -kMax + 2.0*kMax*i/(numPoints-1);
        for (int j = 0; j < numPoints; ++j) {
            double ky = -kMax + 2.0*kMax*j/(numPoints-1);
            double E = Ek(kx, ky);
            if (abs(E - EF) < 0.01 * abs(EF))
                pts.push_back({kx, ky, 0.0, E});
        }
    }
    return pts;
}

void SolidStateSystem::exportDOSFromBandsCSV(const string& filename, const BandResult& bands,
                                               double broadening, int numBins) {
    auto dos = computeDOSFromBands(bands, broadening, numBins);
    ofstream out(filename);
    out << "bin,DOS\n";
    for (int i = 0; i < (int)dos.size(); ++i)
        out << i << "," << dos[i] << "\n";
}

void SolidStateSystem::exportFermiSurface2DCSV(const string& filename, double E0, double t,
                                                  double a, double EF, int numPoints) {
    auto Ek = [&](double kx, double ky) { return tightBindingEnergy2D(E0, t, kx, ky, a); };
    auto pts = computeFermiSurface2D(Ek, EF, numPoints);
    ofstream out(filename);
    out << "kx,ky,E\n";
    for (auto& p : pts) out << p.kx << "," << p.ky << "," << p.energy << "\n";
}

// ============================================================================
//  17. Effective Mass & k·p Theory
// ============================================================================

double SolidStateSystem::effectiveMassFromBand(double d2Edk2) {
    if (abs(d2Edk2) < 1e-50) return 1e30;
    return hbar * hbar / d2Edk2;
}

double SolidStateSystem::kpBandEnergy(double k, double Eg, double P, double mStar) {
    (void)mStar;
    return Eg / 2.0 + sqrt(Eg*Eg/4.0 + P*P*k*k);
}

pair<double, double> SolidStateSystem::kpTwoBand(double k, double Eg, double P) {
    double disc = sqrt(Eg*Eg/4.0 + P*P*hbar*hbar*k*k);
    return {Eg/2.0 - disc, Eg/2.0 + disc};
}

double SolidStateSystem::cyclotronMass(double EF, double dAdE) {
    return hbar * hbar * dAdE / (TWO_PI);
}

void SolidStateSystem::exportKPTheoryCSV(const string& filename, double Eg, double P, int numPoints) {
    ofstream out(filename);
    out << "k,E_conduction,E_valence\n";
    double kMax = 1e10;
    for (int i = 0; i <= numPoints; ++i) {
        double k = -kMax + 2.0*kMax*i/numPoints;
        auto [Ev, Ec] = kpTwoBand(k, Eg, P);
        out << k << "," << Ec << "," << Ev << "\n";
    }
}

// ============================================================================
//  18. Pseudopotential Method
// ============================================================================

double SolidStateSystem::pseudopotentialFormFactor(double q, double V_s, double r_c) {
    if (q < 1e-30) return V_s;
    return V_s * cos(q * r_c);
}

BandResult SolidStateSystem::computePseudopotentialBands(double a, const vector<double>& V_G,
                                                           int numBands, int numPoints) {
    (void)V_G; (void)numBands;
    BandResult result;
    result.k_points.resize(numPoints);
    result.bands.resize(1);
    result.bands[0].resize(numPoints);
    double kMax = PI / a;
    for (int i = 0; i < numPoints; ++i) {
        double k = kMax * i / (numPoints - 1);
        result.k_points[i] = k;
        result.bands[0][i] = hbar*hbar*k*k / (2.0*me*eV);
    }
    return result;
}

void SolidStateSystem::exportPseudopotentialCSV(const string& filename, double a,
                                                  const vector<double>& V_G, int numBands, int numPoints) {
    auto result = computePseudopotentialBands(a, V_G, numBands, numPoints);
    ofstream out(filename);
    out << "k,E_eV\n";
    for (int i = 0; i < numPoints; ++i)
        out << result.k_points[i] << "," << result.bands[0][i] << "\n";
}

// ============================================================================
//  19. Wannier Functions
// ============================================================================

double SolidStateSystem::wannierFunction1D(double x, double a, int nBands,
                                             const vector<double>& blochCoeffs) {
    (void)blochCoeffs;
    double sum = 0.0;
    int N = max(nBands, 1) * 50;
    for (int n = -N; n <= N; ++n) {
        double k = TWO_PI * n / (a * (2*N+1));
        sum += cos(k * x);
    }
    return sum / (2*N+1);
}

double SolidStateSystem::wannierSpread(double a, int nBands) {
    (void)nBands;
    return a / (2.0 * PI);
}

void SolidStateSystem::exportWannierCSV(const string& filename, double a, int nBands, int numPoints) {
    ofstream out(filename);
    out << "x,W\n";
    vector<double> dummy;
    for (int i = 0; i <= numPoints; ++i) {
        double x = -3.0*a + 6.0*a*i/numPoints;
        out << x << "," << wannierFunction1D(x, a, nBands, dummy) << "\n";
    }
}

// ============================================================================
//  20. Band Structure Visualization
// ============================================================================

BandResult SolidStateSystem::computeModelBandStructure(const string& model, double a,
                                                         double param1, double param2,
                                                         int numBands, int numPoints) {
    if (model == "TB")
        return computeTBBands1D(0.0, param1, a, numPoints);
    else if (model == "NFE")
        return computeNFEBands(a, param1, me, numBands, numPoints);
    else {
        BandResult r;
        r.k_points.resize(numPoints);
        r.bands.resize(1);
        r.bands[0].resize(numPoints);
        double kMax = PI / a;
        for (int i = 0; i < numPoints; ++i) {
            double k = -kMax + 2.0*kMax*i/(numPoints-1);
            r.k_points[i] = k;
            r.bands[0][i] = hbar*hbar*k*k/(2.0*me*eV);
        }
        return r;
    }
}

void SolidStateSystem::exportBandStructureCSV(const string& filename, const string& model,
                                                double a, double param1, double param2,
                                                int numBands, int numPoints) {
    auto r = computeModelBandStructure(model, a, param1, param2, numBands, numPoints);
    ofstream out(filename);
    out << "k";
    for (int b = 0; b < (int)r.bands.size(); ++b) out << ",band_" << b;
    out << "\n";
    for (int i = 0; i < (int)r.k_points.size(); ++i) {
        out << r.k_points[i];
        for (int b = 0; b < (int)r.bands.size(); ++b) out << "," << r.bands[b][i];
        out << "\n";
    }
}

// ============================================================================
//  21. Intrinsic Semiconductors
// ============================================================================

double SolidStateSystem::intrinsicCarrierConcentration(double Eg, double T, double mStar_e, double mStar_h) {
    if (T < 1e-12) return 0.0;
    double Nc = 2.0 * pow(mStar_e * kB * T / (TWO_PI * hbar * hbar), 1.5);
    double Nv = 2.0 * pow(mStar_h * kB * T / (TWO_PI * hbar * hbar), 1.5);
    return sqrt(Nc * Nv) * exp(-Eg / (2.0 * kB * T));
}

double SolidStateSystem::intrinsicFermiLevel(double Eg, double T, double mStar_e, double mStar_h) {
    return Eg / 2.0 + 0.75 * kB * T * log(mStar_h / mStar_e);
}

double SolidStateSystem::conductivityIntrinsic(double ni, double mu_e, double mu_h) {
    return ni * e_charge * (mu_e + mu_h);
}

void SolidStateSystem::exportIntrinsicSemiCSV(const string& filename, double Eg, double mStar_e,
                                                double mStar_h, double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,ni,EF\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 10.0) T = 10.0;
        out << T << "," << intrinsicCarrierConcentration(Eg, T, mStar_e, mStar_h) << ","
            << intrinsicFermiLevel(Eg, T, mStar_e, mStar_h) << "\n";
    }
}

// ============================================================================
//  22. Extrinsic Semiconductors
// ============================================================================

double SolidStateSystem::donorIonizationFraction(double Nd, double Ed, double T, double mStar) {
    if (T < 1e-12) return 0.0;
    double Nc = 2.0 * pow(mStar * kB * T / (TWO_PI * hbar * hbar), 1.5);
    double ratio = Nd / Nc * exp(Ed / (kB * T));
    return 1.0 / (1.0 + ratio);
}

double SolidStateSystem::electronConcentration(double Nd, double Na, double Eg, double T,
                                                 double mStar_e, double mStar_h) {
    double ni = intrinsicCarrierConcentration(Eg, T, mStar_e, mStar_h);
    double net = Nd - Na;
    if (net > 0) return 0.5 * (net + sqrt(net * net + 4.0 * ni * ni));
    return ni * ni / (0.5 * (-net + sqrt(net * net + 4.0 * ni * ni)));
}

double SolidStateSystem::holeConcentration(double n, double ni) {
    if (n < 1e-30) return 0.0;
    return ni * ni / n;
}

double SolidStateSystem::extrinsicFermiLevel(double Eg, double T, double Nd, double Na,
                                               double mStar_e, double mStar_h) {
    double ni = intrinsicCarrierConcentration(Eg, T, mStar_e, mStar_h);
    double n = electronConcentration(Nd, Na, Eg, T, mStar_e, mStar_h);
    double Nc = 2.0 * pow(mStar_e * kB * T / (TWO_PI * hbar * hbar), 1.5);
    if (n < 1e-30 || Nc < 1e-30) return Eg / 2.0;
    return kB * T * log(n / Nc) + Eg;  // measured from valence band top
}

void SolidStateSystem::exportExtrinsicSemiCSV(const string& filename, double Eg, double Nd, double Na,
                                                double mStar_e, double mStar_h,
                                                double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,n,p,EF\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 10.0) T = 10.0;
        double ni = intrinsicCarrierConcentration(Eg, T, mStar_e, mStar_h);
        double n = electronConcentration(Nd, Na, Eg, T, mStar_e, mStar_h);
        double p = holeConcentration(n, ni);
        double EF = extrinsicFermiLevel(Eg, T, Nd, Na, mStar_e, mStar_h);
        out << T << "," << n << "," << p << "," << EF << "\n";
    }
}

// ============================================================================
//  23. p-n Junction
// ============================================================================

double SolidStateSystem::builtInVoltage(double Na, double Nd, double ni, double T) {
    return kB * T / e_charge * log(Na * Nd / (ni * ni));
}

double SolidStateSystem::depletionWidth(double Na, double Nd, double epsilon_r, double V_bi, double V_app) {
    double eps = epsilon_r * epsilon0;
    return sqrt(2.0 * eps * (V_bi - V_app) * (Na + Nd) / (e_charge * Na * Nd));
}

PNJunctionResult SolidStateSystem::computePNJunction(double Na, double Nd, double ni, double epsilon_r,
                                                       double T, double V_app, int numPoints) {
    PNJunctionResult res;
    res.V_bi = builtInVoltage(Na, Nd, ni, T);
    res.W = depletionWidth(Na, Nd, epsilon_r, res.V_bi, V_app);
    res.xp = -res.W * Nd / (Na + Nd);
    res.xn = res.W * Na / (Na + Nd);
    double xMin = 3.0 * res.xp;
    double xMax = 3.0 * res.xn;
    res.x.resize(numPoints);
    res.potential.resize(numPoints);
    res.efield.resize(numPoints);
    res.charge.resize(numPoints);
    double eps = epsilon_r * epsilon0;
    for (int i = 0; i < numPoints; ++i) {
        double x = xMin + (xMax - xMin) * i / (numPoints - 1);
        res.x[i] = x;
        if (x < res.xp) {
            res.charge[i] = 0;
            res.efield[i] = 0;
            res.potential[i] = 0;
        } else if (x < 0) {
            res.charge[i] = -e_charge * Na;
            res.efield[i] = -e_charge * Na * (x - res.xp) / eps;
            res.potential[i] = -e_charge * Na * (x - res.xp) * (x - res.xp) / (2.0 * eps);
        } else if (x < res.xn) {
            res.charge[i] = e_charge * Nd;
            res.efield[i] = e_charge * Nd * (res.xn - x) / eps;
            res.potential[i] = res.V_bi - e_charge * Nd * (res.xn - x) * (res.xn - x) / (2.0 * eps);
        } else {
            res.charge[i] = 0;
            res.efield[i] = 0;
            res.potential[i] = res.V_bi;
        }
    }
    return res;
}

double SolidStateSystem::diodeCurrent(double Is, double V, double T, double n_ideal) {
    return Is * (exp(e_charge * V / (n_ideal * kB * T)) - 1.0);
}

void SolidStateSystem::exportPNJunctionCSV(const string& filename, double Na, double Nd, double ni,
                                             double epsilon_r, double T, double V_app, int numPoints) {
    auto res = computePNJunction(Na, Nd, ni, epsilon_r, T, V_app, numPoints);
    ofstream out(filename);
    out << "x,potential,efield,charge\n";
    for (int i = 0; i < numPoints; ++i)
        out << res.x[i] << "," << res.potential[i] << "," << res.efield[i] << "," << res.charge[i] << "\n";
}

// ============================================================================
//  24. Semiconductor Optical Properties
// ============================================================================

double SolidStateSystem::directAbsorptionCoeff(double E_photon, double Eg, double mStar_r) {
    if (E_photon <= Eg) return 0.0;
    return sqrt(2.0 * mStar_r) / (hbar * hbar) * sqrt(E_photon - Eg);
}

double SolidStateSystem::indirectAbsorptionCoeff(double E_photon, double Eg, double E_phonon, double T) {
    double Eabs = E_photon - Eg + E_phonon;
    double Eem = E_photon - Eg - E_phonon;
    double nBE = 1.0 / (exp(E_phonon / (kB * T)) - 1.0);
    double alpha = 0.0;
    if (Eabs > 0) alpha += (nBE + 1.0) * Eabs * Eabs;
    if (Eem > 0) alpha += nBE * Eem * Eem;
    return alpha;
}

double SolidStateSystem::luminescenceSpectrum(double E_photon, double Eg, double T) {
    if (E_photon < Eg) return 0.0;
    return sqrt(E_photon - Eg) * exp(-(E_photon - Eg) / (kB * T));
}

double SolidStateSystem::jointDOS(double E, double Eg, double mStar_r) {
    if (E <= Eg) return 0.0;
    return (1.0 / (2.0 * PI * PI)) * pow(2.0 * mStar_r / (hbar * hbar), 1.5) * sqrt(E - Eg);
}

void SolidStateSystem::exportOpticalPropertiesCSV(const string& filename, double Eg, double mStar_r,
                                                    double T, double Emin, double Emax, int numPoints) {
    ofstream out(filename);
    out << "E_eV,absorption,luminescence,JDOS\n";
    for (int i = 0; i <= numPoints; ++i) {
        double E = Emin + (Emax - Emin) * i / numPoints;
        double E_J = E * eV;
        double Eg_J = Eg * eV;
        out << E << "," << directAbsorptionCoeff(E_J, Eg_J, mStar_r) << ","
            << luminescenceSpectrum(E_J, Eg_J, T) << ","
            << jointDOS(E_J, Eg_J, mStar_r) << "\n";
    }
}

// ============================================================================
//  25. Quantum Wells & Heterostructures
// ============================================================================

double SolidStateSystem::quantumWellEnergy(int n, double Lz, double mStar_w) {
    return n * n * PI * PI * hbar * hbar / (2.0 * mStar_w * Lz * Lz);
}

double SolidStateSystem::quantumWellWavefunction(int n, double z, double Lz) {
    return sqrt(2.0 / Lz) * sin(n * PI * z / Lz);
}

pair<double, double> SolidStateSystem::finiteQWBoundState(double V0, double Lz, double mStar_w,
                                                            double mStar_b, int n) {
    // Simplified: iterate to find bound state energy
    double E = quantumWellEnergy(n, Lz, mStar_w);
    if (E > V0) E = V0 * 0.9;
    // Simple bisection
    double Elow = 0, Ehigh = V0;
    for (int iter = 0; iter < 100; ++iter) {
        E = (Elow + Ehigh) / 2.0;
        double k = sqrt(2.0 * mStar_w * E) / hbar;
        double kappa = sqrt(2.0 * mStar_b * (V0 - E)) / hbar;
        double lhs = k * tan(k * Lz / 2.0);
        double rhs = kappa * mStar_w / mStar_b;
        if (lhs < rhs) Elow = E; else Ehigh = E;
    }
    return {E, E / eV};
}

double SolidStateSystem::superlatticeMiniBand(double t, double k, double d) {
    return -2.0 * t * cos(k * d);
}

void SolidStateSystem::exportQuantumWellCSV(const string& filename, double Lz, double mStar_w,
                                              int maxN, int numPoints) {
    ofstream out(filename);
    out << "z";
    for (int n = 1; n <= maxN; ++n) out << ",psi_" << n;
    out << "\n";
    for (int i = 0; i <= numPoints; ++i) {
        double z = Lz * i / numPoints;
        out << z;
        for (int n = 1; n <= maxN; ++n)
            out << "," << quantumWellWavefunction(n, z, Lz);
        out << "\n";
    }
}

// ============================================================================
//  26. Hall Effect & Magnetoresistance
// ============================================================================

double SolidStateSystem::hallCoefficient(double n_carrier, double charge) {
    return -1.0 / (n_carrier * charge);
}

double SolidStateSystem::hallVoltage(double I, double B, double n_carrier, double t) {
    return I * B / (n_carrier * e_charge * t);
}

double SolidStateSystem::hallMobility(double RH, double sigma) {
    return abs(RH) * sigma;
}

double SolidStateSystem::magnetoresistanceRatio(double mu, double B) {
    return (mu * B) * (mu * B);
}

double SolidStateSystem::cyclotronFrequency(double B, double mStar) {
    return e_charge * B / mStar;
}

void SolidStateSystem::exportHallEffectCSV(const string& filename, double n_carrier, double mu,
                                             double Bmin, double Bmax, int numPoints) {
    ofstream out(filename);
    out << "B,RH,VH,MR_ratio\n";
    double RH = hallCoefficient(n_carrier, e_charge);
    for (int i = 0; i <= numPoints; ++i) {
        double B = Bmin + (Bmax - Bmin) * i / numPoints;
        out << B << "," << RH << "," << hallVoltage(1.0, B, n_carrier, 1e-3) << ","
            << magnetoresistanceRatio(mu, B) << "\n";
    }
}

// ============================================================================
//  27. Boltzmann Transport
// ============================================================================

double SolidStateSystem::drudeCondutivity(double n, double tau, double mStar) {
    return n * e_charge * e_charge * tau / mStar;
}

double SolidStateSystem::drudeThermalConductivity(double n, double tau, double T, double mStar) {
    return PI * PI * n * kB * kB * T * tau / (3.0 * mStar);
}

double SolidStateSystem::wiedemannFranzRatio(double T) {
    return lorenzNumber() * T;
}

double SolidStateSystem::lorenzNumber() {
    return PI * PI / 3.0 * (kB / e_charge) * (kB / e_charge);
}

double SolidStateSystem::relaxationTimeApprox(double mu, double mStar) {
    return mu * mStar / e_charge;
}

void SolidStateSystem::exportBoltzmannTransportCSV(const string& filename, double n, double mStar,
                                                      double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,sigma,kappa_e,L_ratio\n";
    double tau = 1e-14; // typical relaxation time
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        double sigma = drudeCondutivity(n, tau, mStar);
        double kappa_e = drudeThermalConductivity(n, tau, T, mStar);
        double L = (sigma > 0) ? kappa_e / (sigma * T) : 0.0;
        out << T << "," << sigma << "," << kappa_e << "," << L << "\n";
    }
}

// ============================================================================
//  28. Thermoelectric Effects
// ============================================================================

double SolidStateSystem::seebeckCoefficient(double T, double EF, double mStar) {
    (void)mStar;
    if (EF < 1e-30) return 0.0;
    return -PI * PI * kB * kB * T / (3.0 * e_charge * EF);
}

double SolidStateSystem::peltierCoefficient(double S, double T) {
    return S * T;
}

double SolidStateSystem::thomsonCoefficient(double S, double T, double dSdT) {
    (void)S;
    return T * dSdT;
}

double SolidStateSystem::thermoelectricZT(double S, double sigma, double kappa, double T) {
    if (kappa < 1e-30) return 0.0;
    return S * S * sigma * T / kappa;
}

double SolidStateSystem::maxCoolingDeltaT(double ZT, double Tc) {
    return 0.5 * ZT * Tc;
}

void SolidStateSystem::exportThermoelectricCSV(const string& filename, double EF, double mStar,
                                                  double sigma, double kappa,
                                                  double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,S,Pi,ZT\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        double S = seebeckCoefficient(T, EF, mStar);
        double Pi = peltierCoefficient(S, T);
        double ZT = thermoelectricZT(S, sigma, kappa, T);
        out << T << "," << S << "," << Pi << "," << ZT << "\n";
    }
}

// ============================================================================
//  29. Diamagnetism & Paramagnetism
// ============================================================================

double SolidStateSystem::langevinFunction(double x) {
    if (abs(x) < 1e-8) return x / 3.0;
    return 1.0 / tanh(x) - 1.0 / x;
}

double SolidStateSystem::brillouinFunction(double x, double J) {
    double a = (2.0*J+1.0)/(2.0*J);
    double b = 1.0/(2.0*J);
    if (abs(x) < 1e-8) return (J+1.0)*x/(3.0*J);
    return a / tanh(a * x) - b / tanh(b * x);
}

double SolidStateSystem::curieParaMagnetization(double B, double T, double J, double g, int N) {
    if (T < 1e-12) return N * g * muB * J;
    double x = g * muB * J * B / (kB * T);
    return N * g * muB * J * brillouinFunction(x, J);
}

double SolidStateSystem::curieParaSusceptibility(double T, double C) {
    if (T < 1e-12) return 1e30;
    return C / T;
}

double SolidStateSystem::landauDiamagneticSusceptibility(double n, double mStar) {
    return -n * e_charge * e_charge / (12.0 * PI * PI * mStar);
}

double SolidStateSystem::pauliParamagneticSusceptibility(double dosEF) {
    return mu0 * muB * muB * dosEF;
}

void SolidStateSystem::exportParamagnetismCSV(const string& filename, double J, double g, int N,
                                                double Bmax, double T, int numPoints) {
    ofstream out(filename);
    out << "B,M,chi\n";
    for (int i = 0; i <= numPoints; ++i) {
        double B = Bmax * i / numPoints;
        double M = curieParaMagnetization(B, T, J, g, N);
        double chi = (B > 1e-12) ? mu0 * M / B : 0.0;
        out << B << "," << M << "," << chi << "\n";
    }
}

// ============================================================================
//  30. Ferromagnetism
// ============================================================================

double SolidStateSystem::weissFieldMagnetization(double M, double n, double J, double g,
                                                    double T, double lambda_mf) {
    double Beff = mu0 * lambda_mf * M;
    return curieParaMagnetization(Beff, T, J, g, static_cast<int>(n));
}

double SolidStateSystem::curieWeissTemperature(double n, double J, double g, double lambda_mf) {
    double C = n * mu0 * g * g * muB * muB * J * (J + 1.0) / (3.0 * kB);
    return C * lambda_mf;
}

MagnetizationResult SolidStateSystem::solveMeanFieldFerromagnet(double n, double J, double g,
                                                                  double lambda_mf, double Tmin,
                                                                  double Tmax, int numPoints) {
    MagnetizationResult res;
    double Msat = n * g * muB * J;
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        res.temperature.push_back(T);
        // Self-consistent iteration
        double M = Msat * 0.99;
        for (int iter = 0; iter < 200; ++iter) {
            double Beff = mu0 * lambda_mf * M;
            double x = g * muB * J * Beff / (kB * T);
            double Mnew = Msat * brillouinFunction(x, J);
            if (abs(Mnew - M) < 1e-6 * Msat) { M = Mnew; break; }
            M = 0.5 * (M + Mnew);
        }
        res.magnetization.push_back(M);
        // Susceptibility
        double dM = 0.01 * Msat;
        double Bp = mu0 * lambda_mf * (M + dM);
        double xp = g * muB * J * Bp / (kB * T);
        double Mp = Msat * brillouinFunction(xp, J);
        res.susceptibility.push_back((dM > 0) ? (Mp - M) / (mu0 * lambda_mf * dM) : 0.0);
    }
    return res;
}

double SolidStateSystem::spinWaveDispersion(double k, double a, double J, double S) {
    return 4.0 * J * S * (1.0 - cos(k * a));
}

void SolidStateSystem::exportFerromagnetismCSV(const string& filename, double n, double J, double g,
                                                  double lambda_mf, double Tmin, double Tmax, int numPoints) {
    auto res = solveMeanFieldFerromagnet(n, J, g, lambda_mf, Tmin, Tmax, numPoints);
    ofstream out(filename);
    out << "T,M,chi\n";
    for (int i = 0; i < (int)res.temperature.size(); ++i)
        out << res.temperature[i] << "," << res.magnetization[i] << "," << res.susceptibility[i] << "\n";
}

// ============================================================================
//  31. Heisenberg Model
// ============================================================================

double SolidStateSystem::heisenbergExchangeEnergy(double J, double S1_dot_S2) {
    return -2.0 * J * S1_dot_S2;
}

double SolidStateSystem::heisenbergGroundStateEnergy1D(int N, double J, double S) {
    return -2.0 * J * S * S * N;
}

double SolidStateSystem::spinWaveEnergy(double k, double a, double J, double S) {
    return spinWaveDispersion(k, a, J, S);
}

double SolidStateSystem::blochT32Law(double T, double Tc, double M0) {
    if (T >= Tc) return 0.0;
    return M0 * (1.0 - pow(T / Tc, 1.5));
}

void SolidStateSystem::exportHeisenbergCSV(const string& filename, double J, double S, double a,
                                              int numPoints) {
    ofstream out(filename);
    out << "k,E_spinwave\n";
    double kMax = PI / a;
    for (int i = 0; i <= numPoints; ++i) {
        double k = -kMax + 2.0 * kMax * i / numPoints;
        out << k << "," << spinWaveEnergy(k, a, J, S) << "\n";
    }
}

// ============================================================================
//  32. Ising Model (1D/2D)
// ============================================================================

double SolidStateSystem::ising1DPartitionFunction(int N, double J, double B, double T) {
    if (T < 1e-12) return 1.0;
    double K = J / (kB * T);
    double h = muB * B / (kB * T);
    double l1 = exp(K) * cosh(h) + sqrt(exp(2.0*K) * cosh(h) * cosh(h) - 2.0 * sinh(2.0*K));
    double l2 = exp(K) * cosh(h) - sqrt(exp(2.0*K) * cosh(h) * cosh(h) - 2.0 * sinh(2.0*K));
    return pow(l1, N) + pow(l2, N);
}

double SolidStateSystem::ising1DMagnetization(double J, double B, double T) {
    if (T < 1e-12) return 1.0;
    double K = J / (kB * T);
    double h = muB * B / (kB * T);
    return sinh(h) / sqrt(sinh(h)*sinh(h) + exp(-4.0*K));
}

double SolidStateSystem::ising1DCorrelationLength(double J, double T) {
    if (T < 1e-12) return 1e30;
    return -1.0 / log(tanh(J / (kB * T)));
}

double SolidStateSystem::ising2DCriticalTemperature(double J) {
    return 2.0 * J / (kB * log(1.0 + sqrt(2.0)));
}

double SolidStateSystem::ising2DExactMagnetization(double T, double Tc) {
    if (T >= Tc) return 0.0;
    double x = sinh(2.0 * Tc * kB * log(1.0 + sqrt(2.0)) / (2.0 * kB * T));
    double m4 = 1.0 - 1.0 / (x * x * x * x);
    if (m4 < 0) return 0.0;
    return pow(m4, 0.125);
}

IsingResult SolidStateSystem::ising2DMonteCarloSimulation(int L, double J, double Tmin, double Tmax,
                                                            int numTemps, int mcSteps) {
    IsingResult res;
    mt19937 rng(42);
    uniform_int_distribution<int> posDist(0, L-1);
    uniform_real_distribution<double> prob(0.0, 1.0);

    for (int ti = 0; ti < numTemps; ++ti) {
        double T = Tmin + (Tmax - Tmin) * ti / (numTemps - 1);
        if (T < 0.1) T = 0.1;
        res.temperature.push_back(T);
        double beta = 1.0 / (kB * T);

        // Initialize lattice
        vector<vector<int>> spin(L, vector<int>(L, 1));
        double Esum = 0, Msum = 0, E2sum = 0, M2sum = 0;
        int nSamples = 0;

        // Compute initial energy
        double Einit = 0;
        for (int i = 0; i < L; ++i)
            for (int j = 0; j < L; ++j) {
                int sR = spin[(i+1)%L][j];
                int sD = spin[i][(j+1)%L];
                Einit -= J * spin[i][j] * (sR + sD);
            }

        double E = Einit;
        int M_total = L * L;

        for (int step = 0; step < mcSteps; ++step) {
            for (int flip = 0; flip < L*L; ++flip) {
                int i = posDist(rng), j = posDist(rng);
                int nb = spin[(i+1)%L][j] + spin[(i-1+L)%L][j]
                       + spin[i][(j+1)%L] + spin[i][(j-1+L)%L];
                double dE = 2.0 * J * spin[i][j] * nb;
                if (dE <= 0 || prob(rng) < exp(-beta * dE)) {
                    spin[i][j] = -spin[i][j];
                    E += dE;
                    M_total += 2 * spin[i][j];
                }
            }
            if (step > mcSteps / 2) {
                double Eper = E / (L*L);
                double Mper = abs(M_total) / (double)(L*L);
                Esum += Eper; E2sum += Eper*Eper;
                Msum += Mper; M2sum += Mper*Mper;
                nSamples++;
            }
        }
        double Eavg = Esum / nSamples;
        double Mavg = Msum / nSamples;
        res.energy.push_back(Eavg);
        res.magnetization.push_back(Mavg);
        res.specificHeat.push_back((E2sum/nSamples - Eavg*Eavg) * beta * beta);
        res.susceptibility.push_back((M2sum/nSamples - Mavg*Mavg) * beta);
    }
    return res;
}

void SolidStateSystem::exportIsing1DCSV(const string& filename, double J, double Bmax,
                                          double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,M_B0,corr_length\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 0.1) T = 0.1;
        out << T << "," << ising1DMagnetization(J, 0.0, T) << ","
            << ising1DCorrelationLength(J, T) << "\n";
    }
}

void SolidStateSystem::exportIsing2DCSV(const string& filename, int L, double J, double Tmin,
                                          double Tmax, int numTemps, int mcSteps) {
    auto res = ising2DMonteCarloSimulation(L, J, Tmin, Tmax, numTemps, mcSteps);
    ofstream out(filename);
    out << "T,E,M,Cv,chi\n";
    for (int i = 0; i < (int)res.temperature.size(); ++i)
        out << res.temperature[i] << "," << res.energy[i] << "," << res.magnetization[i] << ","
            << res.specificHeat[i] << "," << res.susceptibility[i] << "\n";
}

// ============================================================================
//  33. Antiferromagnetism
// ============================================================================

double SolidStateSystem::neelTemperature(double J, double z, double S) {
    return 2.0 * z * J * S * (S + 1.0) / (3.0 * kB);
}

double SolidStateSystem::afmSusceptibilityAboveTN(double T, double TN, double C) {
    return C / (T + TN);
}

double SolidStateSystem::afmSusceptibilityBelowTN_parallel(double T, double TN, double chi_TN) {
    if (TN < 1e-12) return 0.0;
    return chi_TN * T / TN;
}

double SolidStateSystem::afmSusceptibilityBelowTN_perp(double chi_TN) {
    return chi_TN;
}

double SolidStateSystem::afmSublatticeMagnetization(double T, double TN, double M0) {
    if (T >= TN) return 0.0;
    return M0 * sqrt(1.0 - (T / TN) * (T / TN));
}

void SolidStateSystem::exportAntiferromagnetismCSV(const string& filename, double J, double z,
                                                     double S, double Tmin, double Tmax, int numPoints) {
    double TN = neelTemperature(J, z, S);
    double C = mu0 * muB * muB * S * (S + 1.0) / (3.0 * kB);
    double chi_TN = afmSusceptibilityAboveTN(TN, TN, C);
    ofstream out(filename);
    out << "T,chi_parallel,chi_perp,M_sublattice\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 0.1) T = 0.1;
        double chi_par, chi_perp;
        if (T > TN) {
            chi_par = afmSusceptibilityAboveTN(T, TN, C);
            chi_perp = chi_par;
        } else {
            chi_par = afmSusceptibilityBelowTN_parallel(T, TN, chi_TN);
            chi_perp = afmSusceptibilityBelowTN_perp(chi_TN);
        }
        out << T << "," << chi_par << "," << chi_perp << ","
            << afmSublatticeMagnetization(T, TN, 1.0) << "\n";
    }
}

// ============================================================================
//  34. Ferrimagnetism
// ============================================================================

double SolidStateSystem::ferriNetMagnetization(double MA, double MB) {
    return abs(MA - MB);
}

double SolidStateSystem::ferriCompensationTemperature(double JA, double JB, double SA, double SB) {
    (void)JA; (void)JB;
    // Simplified estimate
    if (abs(SA - SB) < 1e-12) return 0.0;
    return abs(JA * SA - JB * SB) / (kB * abs(SA - SB));
}

MagnetizationResult SolidStateSystem::ferriMeanField(double JA, double JB, double JAB,
                                                       double SA, double SB, int NA, int NB,
                                                       double Tmin, double Tmax, int numPoints) {
    (void)JAB;
    MagnetizationResult res;
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 0.1) T = 0.1;
        res.temperature.push_back(T);
        double MA = NA * muB * SA * brillouinFunction(muB * SA * JA / (kB * T), SA);
        double MB = NB * muB * SB * brillouinFunction(muB * SB * JB / (kB * T), SB);
        res.magnetization.push_back(ferriNetMagnetization(MA, MB));
        res.susceptibility.push_back(0.0); // placeholder
    }
    return res;
}

void SolidStateSystem::exportFerrimagnetismCSV(const string& filename, double JA, double JB, double JAB,
                                                  double SA, double SB, int NA, int NB,
                                                  double Tmin, double Tmax, int numPoints) {
    auto res = ferriMeanField(JA, JB, JAB, SA, SB, NA, NB, Tmin, Tmax, numPoints);
    ofstream out(filename);
    out << "T,M_net,chi\n";
    for (int i = 0; i < (int)res.temperature.size(); ++i)
        out << res.temperature[i] << "," << res.magnetization[i] << "," << res.susceptibility[i] << "\n";
}

// ============================================================================
//  35. Magnetic Domains
// ============================================================================

double SolidStateSystem::domainWallWidth(double A_ex, double K_aniso) {
    return PI * sqrt(A_ex / K_aniso);
}

double SolidStateSystem::domainWallEnergy(double A_ex, double K_aniso) {
    return 4.0 * sqrt(A_ex * K_aniso);
}

double SolidStateSystem::singleDomainCriticalRadius(double A_ex, double K_aniso, double Ms) {
    return 9.0 * domainWallEnergy(A_ex, K_aniso) / (mu0 * Ms * Ms);
}

double SolidStateSystem::hysteresisLoop(double H, double Hc, double Ms, double Mr) {
    // Simplified Stoner-Wohlfarth-like model
    double x = H / Hc;
    if (abs(x) > 5.0) return (x > 0) ? Ms : -Ms;
    return Ms * tanh(x) + Mr * (2.0 / PI) * atan(x);
}

void SolidStateSystem::exportMagneticDomainsCSV(const string& filename, double A_ex, double K_aniso,
                                                   double Ms, double Hmax, int numPoints) {
    double Hc = 2.0 * K_aniso / (mu0 * Ms);
    double Mr = Ms * 0.7;
    ofstream out(filename);
    out << "H,M\n";
    for (int i = 0; i <= numPoints; ++i) {
        double H = -Hmax + 2.0 * Hmax * i / numPoints;
        out << H << "," << hysteresisLoop(H, Hc, Ms, Mr) << "\n";
    }
}

// ============================================================================
//  36. Spin-Orbit Coupling in Solids
// ============================================================================

double SolidStateSystem::rashbaEnergy(double k, double alpha_R, double mStar) {
    return hbar * hbar * k * k / (2.0 * mStar);
    // The Rashba splitting shifts this by +/- alpha_R * k
}

pair<double, double> SolidStateSystem::rashbaSplit(double k, double alpha_R, double mStar) {
    double E0 = hbar * hbar * k * k / (2.0 * mStar);
    return {E0 - alpha_R * abs(k), E0 + alpha_R * abs(k)};
}

double SolidStateSystem::dresselhausEnergy(double k, double beta_D, double mStar) {
    return hbar * hbar * k * k / (2.0 * mStar) + beta_D * k;
}

double SolidStateSystem::spinHallAngle(double sigma_SH, double sigma_xx) {
    if (abs(sigma_xx) < 1e-30) return 0.0;
    return sigma_SH / sigma_xx;
}

void SolidStateSystem::exportSpinOrbitCSV(const string& filename, double alpha_R, double beta_D,
                                            double mStar, int numPoints) {
    ofstream out(filename);
    out << "k,E_rashba_plus,E_rashba_minus,E_dresselhaus\n";
    double kMax = 1e10;
    for (int i = 0; i <= numPoints; ++i) {
        double k = -kMax + 2.0 * kMax * i / numPoints;
        auto [Em, Ep] = rashbaSplit(k, alpha_R, mStar);
        out << k << "," << Ep << "," << Em << "," << dresselhausEnergy(k, beta_D, mStar) << "\n";
    }
}

// ============================================================================
//  37. BCS Theory
// ============================================================================

double SolidStateSystem::bcsGap(double T, double Delta0, double Tc) {
    if (T >= Tc) return 0.0;
    if (T < 1e-12) return Delta0;
    return Delta0 * tanh(1.74 * sqrt(Tc / T - 1.0));
}

double SolidStateSystem::bcsCriticalTemperature(double omegaD, double N_EF, double V_eff) {
    return 1.13 * omegaD * exp(-1.0 / (N_EF * V_eff));
}

double SolidStateSystem::bcsCoherenceLength(double vF, double Delta0) {
    return hbar * vF / (PI * Delta0);
}

double SolidStateSystem::bcsSpecificHeatJump() {
    return 1.43; // Delta C / (gamma * Tc) = 1.43
}

double SolidStateSystem::bcsDOSQuasiparticle(double E, double Delta) {
    if (abs(E) < Delta) return 0.0;
    return abs(E) / sqrt(E * E - Delta * Delta);
}

BCSResult SolidStateSystem::computeBCSGapVsTemp(double omegaD, double N_EF, double V_eff, int numPoints) {
    BCSResult res;
    res.Tc = bcsCriticalTemperature(omegaD, N_EF, V_eff);
    res.Delta0 = 1.764 * kB * res.Tc;
    res.lambda_BCS = N_EF * V_eff;
    for (int i = 0; i <= numPoints; ++i) {
        double T = res.Tc * 1.2 * i / numPoints;
        res.temperature.push_back(T);
        res.gap.push_back(bcsGap(T, res.Delta0, res.Tc));
    }
    return res;
}

void SolidStateSystem::exportBCSTheoryCSV(const string& filename, double omegaD, double N_EF,
                                             double V_eff, int numPoints) {
    auto res = computeBCSGapVsTemp(omegaD, N_EF, V_eff, numPoints);
    ofstream out(filename);
    out << "T,Delta,T_over_Tc\n";
    for (int i = 0; i < (int)res.temperature.size(); ++i)
        out << res.temperature[i] << "," << res.gap[i] << ","
            << ((res.Tc > 0) ? res.temperature[i]/res.Tc : 0.0) << "\n";
}

// ============================================================================
//  38. Ginzburg-Landau Theory
// ============================================================================

double SolidStateSystem::glCoherenceLength(double xi0, double T, double Tc) {
    double t = T / Tc;
    if (t >= 1.0) return 1e30;
    return xi0 / sqrt(1.0 - t);
}

double SolidStateSystem::glPenetrationDepth(double lambda0, double T, double Tc) {
    double t = T / Tc;
    if (t >= 1.0) return 1e30;
    return lambda0 / sqrt(1.0 - t);
}

double SolidStateSystem::glKappaParameter(double lambda, double xi) {
    if (xi < 1e-30) return 0.0;
    return lambda / xi;
}

double SolidStateSystem::glCriticalField(double Hc0, double T, double Tc) {
    double t = T / Tc;
    if (t >= 1.0) return 0.0;
    return Hc0 * (1.0 - t * t);
}

double SolidStateSystem::glUpperCriticalField(double Hc, double kappa) {
    return sqrt(2.0) * kappa * Hc;
}

double SolidStateSystem::glLowerCriticalField(double Hc, double kappa) {
    if (kappa < 1e-30) return Hc;
    return Hc * log(kappa) / (sqrt(2.0) * kappa);
}

GLResult SolidStateSystem::computeGLProfile(double xi, double lambda, int numPoints) {
    GLResult res;
    res.xi = xi;
    res.lambda_L = lambda;
    res.kappa = glKappaParameter(lambda, xi);
    double xMax = 5.0 * max(xi, lambda);
    for (int i = 0; i <= numPoints; ++i) {
        double x = xMax * i / numPoints;
        res.x.push_back(x);
        res.orderParameter.push_back(tanh(x / (sqrt(2.0) * xi)));
        res.magneticField.push_back(exp(-x / lambda));
    }
    return res;
}

void SolidStateSystem::exportGLTheoryCSV(const string& filename, double xi0, double lambda0,
                                            double Tc, double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,xi,lambda,kappa,Hc\n";
    double Hc0 = 1e5; // placeholder
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        double xi = glCoherenceLength(xi0, T, Tc);
        double lam = glPenetrationDepth(lambda0, T, Tc);
        double kap = glKappaParameter(lam, xi);
        double Hc = glCriticalField(Hc0, T, Tc);
        out << T << "," << xi << "," << lam << "," << kap << "," << Hc << "\n";
    }
}

// ============================================================================
//  39. Type I & II Superconductors
// ============================================================================

double SolidStateSystem::thermodynamicCriticalField(double Hc0, double T, double Tc) {
    return glCriticalField(Hc0, T, Tc);
}

bool SolidStateSystem::isTypeII(double kappa) {
    return kappa > 1.0 / sqrt(2.0);
}

double SolidStateSystem::fluxQuantum() {
    return h_planck / (2.0 * e_charge);
}

double SolidStateSystem::vortexSpacing(double B) {
    double Phi0 = fluxQuantum();
    return sqrt(2.0 * Phi0 / (sqrt(3.0) * B));
}

double SolidStateSystem::mixedStateB(double H, double Hc1, double Hc2) {
    if (H <= Hc1) return 0.0;
    if (H >= Hc2) return mu0 * H;
    return mu0 * H * (H - Hc1) / (Hc2 - Hc1);
}

void SolidStateSystem::exportSCTypesCSV(const string& filename, double Hc0, double kappa,
                                           double Tc, int numPoints) {
    ofstream out(filename);
    out << "T,Hc,Hc1,Hc2\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tc * i / numPoints;
        double Hc = glCriticalField(Hc0, T, Tc);
        out << T << "," << Hc << "," << glLowerCriticalField(Hc, kappa) << ","
            << glUpperCriticalField(Hc, kappa) << "\n";
    }
}

// ============================================================================
//  40. Josephson Effect
// ============================================================================

double SolidStateSystem::josephsonCriticalCurrent(double Delta, double Rn) {
    return PI * Delta / (2.0 * e_charge * Rn);
}

double SolidStateSystem::josephsonDCCurrent(double Ic, double phi) {
    return Ic * sin(phi);
}

double SolidStateSystem::josephsonACFrequency(double V) {
    return 2.0 * e_charge * V / h_planck;
}

double SolidStateSystem::josephsonPlasmaFrequency(double Ic, double C) {
    return sqrt(2.0 * e_charge * Ic / (hbar * C));
}

double SolidStateSystem::shapiroStepVoltage(int n, double freq) {
    return n * h_planck * freq / (2.0 * e_charge);
}

double SolidStateSystem::squidFluxResponse(double Ic, double flux, double L) {
    double Phi0 = fluxQuantum();
    return 2.0 * Ic * abs(cos(PI * flux / Phi0));
}

void SolidStateSystem::exportJosephsonCSV(const string& filename, double Ic, double Rn,
                                             double Vmax, int numPoints) {
    ofstream out(filename);
    out << "phi,I_dc,V_ac_freq\n";
    for (int i = 0; i <= numPoints; ++i) {
        double phi = TWO_PI * i / numPoints;
        double V = Vmax * i / numPoints;
        out << phi << "," << josephsonDCCurrent(Ic, phi) << ","
            << josephsonACFrequency(V) << "\n";
    }
}

// ============================================================================
//  41. Electron-Phonon Interaction
// ============================================================================

double SolidStateSystem::eliashbergCouplingLambda(double N_EF, double g2, double omegaPh) {
    if (omegaPh < 1e-30) return 0.0;
    return 2.0 * N_EF * g2 / omegaPh;
}

double SolidStateSystem::mcMillanTc(double omegaLog, double lambda, double muStar) {
    double denom = lambda - muStar * (1.0 + 0.62 * lambda);
    if (denom < 1e-12) return 0.0;
    return (omegaLog / 1.20) * exp(-1.04 * (1.0 + lambda) / denom);
}

double SolidStateSystem::electronPhononScatteringRate(double T, double thetaD) {
    if (T < 1e-12) return 0.0;
    return (T / thetaD); // simplified
}

double SolidStateSystem::resistivityBlochGruneisen(double T, double thetaD, double rho0, double K) {
    if (T < 1e-12) return rho0;
    double x = thetaD / T;
    // Simplified BG: rho = rho0 + K * (T/thetaD)^5 * integral
    double integral = debyeFunction(x, 5) * pow(x, 5);
    return rho0 + K * pow(T / thetaD, 5.0) * integral;
}

void SolidStateSystem::exportElectronPhononCSV(const string& filename, double thetaD, double rho0,
                                                  double K, double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "T,rho\n";
    for (int i = 0; i <= numPoints; ++i) {
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 1.0) T = 1.0;
        out << T << "," << resistivityBlochGruneisen(T, thetaD, rho0, K) << "\n";
    }
}

// ============================================================================
//  42. Anderson Localization
// ============================================================================

double SolidStateSystem::andersonLocalizationLength1D(double mfp) {
    return 2.0 * mfp;
}

double SolidStateSystem::andersonLocalizationLength2D(double mfp) {
    return mfp * exp(PI * mfp / 2.0); // Thouless scaling
}

double SolidStateSystem::ioffRegelCriterion(double kF, double mfp) {
    return kF * mfp;
}

double SolidStateSystem::thoulessEnergy(double D, double L) {
    return hbar * D / (L * L);
}

double SolidStateSystem::conductanceDimensionless(double g0, double L, double xi) {
    return g0 * exp(-L / xi);
}

void SolidStateSystem::exportAndersonCSV(const string& filename, double mfp, double kF,
                                           double Lmin, double Lmax, int numPoints) {
    ofstream out(filename);
    out << "L,xi_1D,xi_2D,IR,g\n";
    double xi1 = andersonLocalizationLength1D(mfp);
    double xi2 = andersonLocalizationLength2D(mfp);
    double ir = ioffRegelCriterion(kF, mfp);
    for (int i = 0; i <= numPoints; ++i) {
        double L = Lmin + (Lmax - Lmin) * i / numPoints;
        out << L << "," << xi1 << "," << xi2 << "," << ir << ","
            << conductanceDimensionless(1.0, L, xi1) << "\n";
    }
}

// ============================================================================
//  43. Quantum Hall Effect (Extended)
// ============================================================================

double SolidStateSystem::landauLevel(int n, double B, double mStar) {
    double omega_c = e_charge * B / mStar;
    return hbar * omega_c * (n + 0.5);
}

double SolidStateSystem::magneticLength(double B) {
    return sqrt(hbar / (e_charge * B));
}

double SolidStateSystem::fillingFactor(double n2D, double B) {
    return n2D * h_planck / (e_charge * B);
}

double SolidStateSystem::hallConductance(double nu) {
    return nu * e_charge * e_charge / h_planck;
}

double SolidStateSystem::hallResistance(double nu) {
    if (abs(nu) < 1e-12) return 0.0;
    return h_planck / (nu * e_charge * e_charge);
}

double SolidStateSystem::longitudinalResistance(double B, double n2D, double mStar, double broadening) {
    double nu = fillingFactor(n2D, B);
    double nuRound = round(nu);
    double delta = nu - nuRound;
    return exp(-delta * delta / (2.0 * broadening * broadening));
}

double SolidStateSystem::fractionalFilling(int p, int q) {
    return static_cast<double>(p) / q;
}

double SolidStateSystem::laughlinGapEstimate(double B, double mStar) {
    double omega_c = e_charge * B / mStar;
    return 0.1 * hbar * omega_c; // rough estimate
}

void SolidStateSystem::exportQuantumHallCSV(const string& filename, double n2D, double mStar,
                                               double Bmin, double Bmax, int numPoints) {
    ofstream out(filename);
    out << "B,nu,sigma_xy,R_xy,R_xx\n";
    for (int i = 1; i <= numPoints; ++i) {
        double B = Bmin + (Bmax - Bmin) * i / numPoints;
        double nu = fillingFactor(n2D, B);
        out << B << "," << nu << "," << hallConductance(nu) << ","
            << hallResistance(nu) << "," << longitudinalResistance(B, n2D, mStar, 0.1) << "\n";
    }
}

// ============================================================================
//  44. Topological Insulators (Toy)
// ============================================================================

double SolidStateSystem::sshEnergy(double k, double v, double w) {
    double dx = v + w * cos(k);
    double dy = w * sin(k);
    return sqrt(dx*dx + dy*dy);
}

pair<double, double> SolidStateSystem::sshBands(double k, double v, double w) {
    double E = sshEnergy(k, v, w);
    return {-E, E};
}

double SolidStateSystem::windingNumber(double v, double w) {
    return (abs(w) > abs(v)) ? 1.0 : 0.0;
}

int SolidStateSystem::cherNumberTwoLevel(double m, double t) {
    if (m > 0 && m < 4.0 * t) return 1;
    if (m > 4.0 * t && m < 8.0 * t) return -1;
    return 0;
}

BandResult SolidStateSystem::computeSSHBands(double v, double w, int numPoints) {
    BandResult result;
    result.k_points.resize(numPoints);
    result.bands.resize(2);
    result.bands[0].resize(numPoints);
    result.bands[1].resize(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        double k = -PI + TWO_PI * i / (numPoints - 1);
        result.k_points[i] = k;
        auto [Em, Ep] = sshBands(k, v, w);
        result.bands[0][i] = Em;
        result.bands[1][i] = Ep;
    }
    return result;
}

void SolidStateSystem::exportTopologicalCSV(const string& filename, double v, double w, int numPoints) {
    auto result = computeSSHBands(v, w, numPoints);
    ofstream out(filename);
    out << "k,E_lower,E_upper\n";
    for (int i = 0; i < numPoints; ++i)
        out << result.k_points[i] << "," << result.bands[0][i] << "," << result.bands[1][i] << "\n";
}

// ============================================================================
//  45. Graphene & Dirac Materials
// ============================================================================

double SolidStateSystem::grapheneDispersion(double kx, double ky, double t, double a) {
    double fx = cos(kx * a) + 2.0 * cos(kx * a / 2.0) * cos(ky * a * sqrt(3.0) / 2.0);
    return t * sqrt(3.0 + fx);
}

pair<double, double> SolidStateSystem::grapheneBands(double kx, double ky, double t, double a) {
    double E = grapheneDispersion(kx, ky, t, a);
    return {-E, E};
}

double SolidStateSystem::grapheneDOS(double E, double t) {
    if (abs(E) > 3.0 * abs(t)) return 0.0;
    return 2.0 * abs(E) / (PI * t * t); // linear DOS near Dirac point
}

double SolidStateSystem::grapheneFermiVelocity(double t, double a) {
    return 3.0 * t * a / (2.0 * hbar);
}

double SolidStateSystem::grapheneLandauLevel(int n, double B) {
    double vF = 1.0e6; // m/s for graphene
    double lB = magneticLength(B);
    return vF * sqrt(2.0 * hbar * e_charge * B * abs(n));
}

double SolidStateSystem::grapheneCyclotronMass(double EF, double vF) {
    return abs(EF) / (vF * vF);
}

void SolidStateSystem::exportGrapheneCSV(const string& filename, double t, double a, int numPoints) {
    ofstream out(filename);
    out << "k,E_plus,E_minus\n";
    double kMax = 2.0 * PI / a;
    for (int i = 0; i <= numPoints; ++i) {
        double kx = kMax * (i - numPoints/2.0) / numPoints;
        auto [Em, Ep] = grapheneBands(kx, 0.0, t, a);
        out << kx << "," << Ep << "," << Em << "\n";
    }
}

// ============================================================================
//  46. Plasmons & Dielectric Function
// ============================================================================

double SolidStateSystem::plasmaFrequency(double n, double mStar, double epsilon_inf) {
    return sqrt(n * e_charge * e_charge / (mStar * epsilon_inf * epsilon0));
}

complex<double> SolidStateSystem::drudePermittivity(double omega, double omega_p, double gamma) {
    if (omega < 1e-30) return complex<double>(1e30, 0);
    return 1.0 - omega_p * omega_p / (omega * omega + complex<double>(0, gamma * omega));
}

double SolidStateSystem::lindhardDielectric(double q, double omega, double kF, double vF, double mStar) {
    (void)omega; (void)mStar;
    if (q < 1e-30) return 1.0;
    double x = q / (2.0 * kF);
    // Static Lindhard at T=0: simplified
    return 1.0 + (e_charge * e_charge * kF) / (PI * PI * epsilon0 * hbar * vF)
           * (0.5 + (1.0 - x*x) / (4.0*x) * log(abs((1.0+x)/(1.0-x+1e-30))));
}

double SolidStateSystem::thomasFermiScreeningLength(double n, double mStar) {
    double EF = fermiEnergy(n, mStar);
    double dos = freeElectronDOS3D(EF, mStar);
    return sqrt(epsilon0 / (e_charge * e_charge * dos));
}

double SolidStateSystem::plasmonDispersion(double q, double omega_p, double vF) {
    return sqrt(omega_p * omega_p + 3.0/5.0 * vF * vF * q * q);
}

void SolidStateSystem::exportPlasmonCSV(const string& filename, double n, double mStar,
                                           double gamma, double omegaMax, int numPoints) {
    double omega_p = plasmaFrequency(n, mStar, 1.0);
    ofstream out(filename);
    out << "omega,eps_real,eps_imag\n";
    for (int i = 1; i <= numPoints; ++i) {
        double omega = omegaMax * i / numPoints;
        auto eps = drudePermittivity(omega, omega_p, gamma);
        out << omega << "," << eps.real() << "," << eps.imag() << "\n";
    }
}

// ============================================================================
//  47. Excitons
// ============================================================================

double SolidStateSystem::excitonBindingEnergy(int n, double mStar_r, double epsilon_r) {
    double Ry_exc = mStar_r * e_charge * e_charge * e_charge * e_charge
                  / (2.0 * hbar * hbar * (4.0*PI*epsilon0*epsilon_r) * (4.0*PI*epsilon0*epsilon_r));
    return Ry_exc / (n * n);
}

double SolidStateSystem::excitonBohrRadius(double mStar_r, double epsilon_r) {
    return 4.0 * PI * epsilon0 * epsilon_r * hbar * hbar / (mStar_r * e_charge * e_charge);
}

double SolidStateSystem::excitonReducedMass(double mStar_e, double mStar_h) {
    return mStar_e * mStar_h / (mStar_e + mStar_h);
}

double SolidStateSystem::mottDensity(double aB_exc, int dim) {
    if (dim == 2) return 1.0 / (PI * aB_exc * aB_exc);
    return 1.0 / (aB_exc * aB_exc * aB_exc); // 3D
}

double SolidStateSystem::excitonAbsorption(double E, double Eg, double Eb, double broadening) {
    // Continuum + bound states (simplified)
    double alpha = 0.0;
    // Bound states
    for (int n = 1; n <= 5; ++n) {
        double En = Eg - Eb / (n * n);
        double x = (E - En) / broadening;
        alpha += (1.0 / (n * n * n)) * exp(-0.5 * x * x) / (broadening * sqrt(TWO_PI));
    }
    // Continuum
    if (E > Eg) alpha += sqrt(E - Eg);
    return alpha;
}

void SolidStateSystem::exportExcitonCSV(const string& filename, double mStar_e, double mStar_h,
                                           double epsilon_r, double Eg, int numPoints) {
    double mr = excitonReducedMass(mStar_e, mStar_h);
    double Eb = excitonBindingEnergy(1, mr, epsilon_r);
    ofstream out(filename);
    out << "E_eV,absorption\n";
    for (int i = 0; i <= numPoints; ++i) {
        double E = (Eg - 3.0*Eb/eV) + (6.0*Eb/eV) * i / numPoints;
        out << E << "," << excitonAbsorption(E*eV, Eg*eV, Eb, 0.01*eV) << "\n";
    }
}

// ============================================================================
//  48. Polaritons
// ============================================================================

pair<double, double> SolidStateSystem::polaritonDispersion(double k, double omega_c, double omega_exc,
                                                             double g_coupling) {
    double avg = (omega_c + omega_exc) / 2.0;
    double delta = omega_c - omega_exc;
    double disc = sqrt(delta * delta / 4.0 + g_coupling * g_coupling);
    return {avg - disc, avg + disc};
}

double SolidStateSystem::rabiSplitting(double g_coupling) {
    return 2.0 * g_coupling;
}

double SolidStateSystem::polaritonHopfieldCoeff(double detuning, double g_coupling) {
    double Omega = sqrt(detuning * detuning + 4.0 * g_coupling * g_coupling);
    return 0.5 * (1.0 + detuning / Omega);
}

double SolidStateSystem::cavityPhotonEnergy(double k, double omega0, double mStar_ph) {
    return sqrt(omega0 * omega0 + hbar * hbar * k * k / (mStar_ph * mStar_ph));
}

void SolidStateSystem::exportPolaritonCSV(const string& filename, double omega_c, double omega_exc,
                                             double g_coupling, int numPoints) {
    ofstream out(filename);
    out << "detuning,E_lower,E_upper\n";
    for (int i = 0; i <= numPoints; ++i) {
        double det = -5.0*g_coupling + 10.0*g_coupling*i/numPoints;
        auto [El, Eu] = polaritonDispersion(0, omega_c + det, omega_exc, g_coupling);
        out << det << "," << El << "," << Eu << "\n";
    }
}

// ============================================================================
//  49. Amorphous Solids
// ============================================================================

double SolidStateSystem::radialDistributionFunction(double r, double r0, double sigma, double n0) {
    if (r < 1e-30) return 0.0;
    double g = 1.0;
    // First coordination shell
    g += (n0 / (r * sqrt(TWO_PI) * sigma)) * exp(-0.5 * (r - r0) * (r - r0) / (sigma * sigma));
    // Second shell
    g += 0.5 * (n0 / (r * sqrt(TWO_PI) * sigma)) * exp(-0.5 * (r - 1.7*r0) * (r - 1.7*r0) / (sigma * sigma * 2.0));
    return g;
}

double SolidStateSystem::tunnelingSystems(double T, double P0) {
    return P0 * T; // linear in T (two-level systems)
}

double SolidStateSystem::amorphousSpecificHeat(double T, double alpha_c, double beta_c) {
    return alpha_c * T + beta_c * T * T * T; // linear + Debye
}

double SolidStateSystem::amorphousThermalConductivity(double T, double kMin, double kPhonon) {
    return kMin + kPhonon * T * T; // plateau + T^2
}

void SolidStateSystem::exportAmorphousSolidsCSV(const string& filename, double r0, double sigma,
                                                   double n0, double Tmin, double Tmax, int numPoints) {
    ofstream out(filename);
    out << "r,g_r,T,Cv_amorphous,kappa_amorphous\n";
    for (int i = 0; i <= numPoints; ++i) {
        double r = 0.5 * r0 + 3.0 * r0 * i / numPoints;
        double T = Tmin + (Tmax - Tmin) * i / numPoints;
        if (T < 0.1) T = 0.1;
        out << r << "," << radialDistributionFunction(r, r0, sigma, n0) << ","
            << T << "," << amorphousSpecificHeat(T, 1e-3, 1e-6) << ","
            << amorphousThermalConductivity(T, 0.01, 1e-4) << "\n";
    }
}

// ============================================================================
//  50. Mesoscopic Transport
// ============================================================================

double SolidStateSystem::landauerConductance(double T_transmission, int numChannels) {
    return numChannels * T_transmission * conductanceQuantum();
}

double SolidStateSystem::conductanceQuantum() {
    return 2.0 * e_charge * e_charge / h_planck;
}

double SolidStateSystem::coulombBlockadeChargingEnergy(double C) {
    return e_charge * e_charge / (2.0 * C);
}

double SolidStateSystem::coulombDiamondBoundary(double Vg, double Vsd, double Ec, int n) {
    (void)Vsd;
    return Ec * (2.0 * n + 1.0) - e_charge * Vg;
}

TransportResult SolidStateSystem::computeLandauerTransport(double barrierHeight, double barrierWidth,
                                                              double mStar, double Emin, double Emax,
                                                              int numPoints) {
    TransportResult res;
    double G0 = conductanceQuantum();
    for (int i = 0; i <= numPoints; ++i) {
        double E = Emin + (Emax - Emin) * i / numPoints;
        res.energy.push_back(E);
        double T;
        if (E >= barrierHeight) {
            T = 1.0;
        } else {
            double kappa = sqrt(2.0 * mStar * (barrierHeight - E)) / hbar;
            T = exp(-2.0 * kappa * barrierWidth);
        }
        res.transmission.push_back(T);
        res.conductance.push_back(G0 * T);
    }
    return res;
}

double SolidStateSystem::abFlux(double B, double area) {
    return B * area;
}

double SolidStateSystem::aharonovBohmOscillationPeriod(double area) {
    return fluxQuantum() / area;
}

void SolidStateSystem::exportMesoscopicCSV(const string& filename, double barrierH, double barrierW,
                                              double mStar, double Emin, double Emax, int numPoints) {
    auto res = computeLandauerTransport(barrierH, barrierW, mStar, Emin, Emax, numPoints);
    ofstream out(filename);
    out << "E,T,G\n";
    for (int i = 0; i < (int)res.energy.size(); ++i)
        out << res.energy[i] << "," << res.transmission[i] << "," << res.conductance[i] << "\n";
}
