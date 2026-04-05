#pragma once
// ============================================================================
//  SolidStateSystem.h  —  Condensed-matter physics engine (50 simulations)
// ============================================================================

#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <array>
#include <tuple>
#include <functional>

// ── Result Structs ──────────────────────────────────────────────────────────

struct LatticeVector {
    double x, y, z;
};

struct BravaisLattice {
    LatticeVector a1, a2, a3;
    std::string name;
};

struct DiffractionPeak {
    int h, k, l;
    double d_spacing;
    double twoTheta;
    double intensity;
};

struct PhononBranch {
    std::vector<double> k_points;
    std::vector<double> omega;
    std::string label;   // "acoustic" or "optical"
};

struct BandResult {
    std::vector<double> k_points;
    std::vector<std::vector<double>> bands;  // bands[bandIndex][kIndex]
};

struct FermiSurfacePoint {
    double kx, ky, kz;
    double energy;
};

struct PNJunctionResult {
    double V_bi;           // built-in voltage
    double W;              // depletion width
    double xn, xp;         // depletion edges
    std::vector<double> x;
    std::vector<double> potential;
    std::vector<double> efield;
    std::vector<double> charge;
};

struct MagnetizationResult {
    std::vector<double> temperature;
    std::vector<double> magnetization;
    std::vector<double> susceptibility;
};

struct IsingResult {
    std::vector<double> temperature;
    std::vector<double> energy;
    std::vector<double> magnetization;
    std::vector<double> specificHeat;
    std::vector<double> susceptibility;
};

struct BCSResult {
    double Delta0;          // gap at T=0
    double Tc;              // critical temperature
    double lambda_BCS;      // coupling constant
    std::vector<double> temperature;
    std::vector<double> gap;
};

struct GLResult {
    std::vector<double> x;
    std::vector<double> orderParameter;
    std::vector<double> magneticField;
    double xi;              // coherence length
    double lambda_L;        // penetration depth
    double kappa;           // GL parameter
};

struct TransportResult {
    std::vector<double> energy;
    std::vector<double> transmission;
    std::vector<double> conductance;
};

// ============================================================================
//  SolidStateSystem  —  Core physics class
// ============================================================================

class SolidStateSystem {
public:
    // Constructor
    SolidStateSystem(const std::string& materialName = "Generic",
                     double latticeConstant = 5.0e-10,
                     double effectiveMass = 9.1093837015e-31);

    // Accessors
    const std::string& getMaterialName() const { return m_materialName; }
    double getLatticeConstant() const { return m_a; }
    double getEffectiveMass() const { return m_mStar; }
    void setMaterial(const std::string& name, double a, double mStar);

    // ========================================================================
    //  1. Bravais Lattices
    // ========================================================================
    static BravaisLattice getBravaisLattice(int type, double a, double b = 0, double c = 0,
                                             double alpha = 90, double beta = 90, double gamma = 90);
    static LatticeVector computeReciprocalVector(const LatticeVector& a1,
                                                  const LatticeVector& a2,
                                                  const LatticeVector& a3, int which);
    static double latticeSpacing(int h, int k, int l, double a);
    void exportBravaisLatticeCSV(const std::string& filename, int type, double a);

    // ========================================================================
    //  2. Reciprocal Lattice
    // ========================================================================
    static LatticeVector reciprocalLatticePoint(const LatticeVector& b1,
                                                 const LatticeVector& b2,
                                                 const LatticeVector& b3,
                                                 int n1, int n2, int n3);
    static double brillouinZoneBoundary1D(double a);
    void exportReciprocalLatticeCSV(const std::string& filename, double a, int maxIndex);

    // ========================================================================
    //  3. X-Ray Diffraction
    // ========================================================================
    static double braggAngle(double d, double wavelength, int n = 1);
    static double structureFactor(int h, int k, int l, const std::string& latticeType);
    static std::vector<DiffractionPeak> computeDiffractionPattern(double a, double wavelength,
                                                                    const std::string& latticeType, int maxHKL);
    static double scatteringVector(double theta, double wavelength);
    static double debyeWallerFactor(double B, double theta, double wavelength);
    void exportDiffractionCSV(const std::string& filename, double a, double wavelength,
                               const std::string& latticeType, int maxHKL);

    // ========================================================================
    //  4. Lattice Vibrations — 1D Monatomic
    // ========================================================================
    static double monatomicDispersion(double k, double a, double springK, double mass);
    static double monatomicGroupVelocity(double k, double a, double springK, double mass);
    static double monatomicDOS1D(double omega, double omegaMax);
    static double debyeCutoffFrequency(double springK, double mass);
    void exportMonatomicDispersionCSV(const std::string& filename, double a,
                                       double springK, double mass, int numPoints);

    // ========================================================================
    //  5. Lattice Vibrations — 1D Diatomic
    // ========================================================================
    static std::pair<double, double> diatomicDispersion(double k, double a,
                                                         double springK, double m1, double m2);
    static double diatomicBandGap(double springK, double m1, double m2);
    void exportDiatomicDispersionCSV(const std::string& filename, double a,
                                      double springK, double m1, double m2, int numPoints);

    // ========================================================================
    //  6. Debye Model
    // ========================================================================
    static double debyeSpecificHeat(double T, double thetaD);
    static double debyeEnergy(double T, double thetaD, int N);
    static double debyeEntropy(double T, double thetaD, int N);
    static double debyeFunction(double x, int n);
    void exportDebyeModelCSV(const std::string& filename, double thetaD, int N,
                              double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  7. Einstein Model
    // ========================================================================
    static double einsteinSpecificHeat(double T, double thetaE);
    static double einsteinEnergy(double T, double thetaE, int N);
    static double einsteinEntropy(double T, double thetaE, int N);
    void exportEinsteinModelCSV(const std::string& filename, double thetaE, int N,
                                 double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  8. Phonon Density of States
    // ========================================================================
    static double debyeDOS3D(double omega, double omegaD);
    static double vanHoveSingularity1D(double omega, double omegaMax);
    static std::vector<double> computePhononDOS(const std::vector<double>& frequencies,
                                                  double broadening, int numBins);
    void exportPhononDOSCSV(const std::string& filename, double omegaD, int numPoints);

    // ========================================================================
    //  9. Thermal Expansion
    // ========================================================================
    static double grueneisenParameter(double bulkModulus, double Cv, double volume, double alpha);
    static double thermalExpansionCoeff(double gamma, double Cv, double bulkModulus, double volume);
    static double anharmonicPotential(double x, double k2, double k3, double k4);
    void exportThermalExpansionCSV(const std::string& filename, double thetaD, double gamma,
                                    double V0, double B, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  10. Lattice Thermal Conductivity
    // ========================================================================
    static double kineticTheoryThermalConductivity(double Cv, double v, double mfp);
    static double umklappScatteringRate(double T, double thetaD, double omega);
    static double phononMeanFreePath(double T, double thetaD);
    void exportThermalConductivityCSV(const std::string& filename, double thetaD,
                                       double v_sound, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  11. Free Electron Model
    // ========================================================================
    static double freeElectronEnergy(double k, double mStar);
    static double fermiEnergy(double n_density, double mStar);
    static double fermiWavevector(double EF, double mStar);
    static double fermiVelocity(double EF, double mStar);
    static double fermiTemperature(double EF);
    static double freeElectronDOS3D(double E, double mStar);
    static double freeElectronDOS2D(double mStar);
    static double freeElectronDOS1D(double E, double mStar);
    static double fermiDiracDistribution(double E, double EF, double T);
    static double electronicSpecificHeat(double T, double EF, double n_density);
    void exportFreeElectronCSV(const std::string& filename, double n_density, double mStar,
                                double Emax, int numPoints);

    // ========================================================================
    //  12. Nearly Free Electron Model
    // ========================================================================
    static std::pair<double, double> NFEBandEdges(double k, double G, double V_G, double mStar);
    static double NFEBandGap(double V_G);
    static BandResult computeNFEBands(double a, double V_G, double mStar, int numBands, int numPoints);
    void exportNFEBandsCSV(const std::string& filename, double a, double V_G,
                            double mStar, int numBands, int numPoints);

    // ========================================================================
    //  13. Tight-Binding (Extended)
    // ========================================================================
    static double tightBindingEnergy1D(double E0, double t, double k, double a);
    static double tightBindingEnergy2D(double E0, double t, double kx, double ky, double a);
    static double tightBindingEffectiveMass(double t, double a, double mStar);
    static double tightBindingDOS1D(double E, double E0, double t);
    static double tightBindingBandwidth(double t);
    static BandResult computeTBBands1D(double E0, double t, double a, int numPoints);
    static BandResult computeTBBands2D(double E0, double t, double a, int numKPoints);
    void exportTBBandsCSV(const std::string& filename, double E0, double t, double a, int numPoints);

    // ========================================================================
    //  14. Kronig-Penney (Extended)
    // ========================================================================
    static double kronigPenneyDispersion(double E, double V0, double a, double b, double mStar);
    static double kronigPenneyDeltaDispersion(double E, double P, double a);
    static BandResult computeKPBands(double V0, double a, double b, double mStar,
                                      int numBands, int numPoints);
    void exportKPBandsCSV(const std::string& filename, double V0, double a, double b,
                           double mStar, int numBands, int numPoints);

    // ========================================================================
    //  15. APW / OPW Methods
    // ========================================================================
    static double muffinTinPotential(double r, double rMT, double V0);
    static double apwPhaseShift(int l, double E, double rMT, double V0);
    static BandResult computeAPWBands(double a, double rMT, double V0, int lMax, int numPoints);
    void exportAPWBandsCSV(const std::string& filename, double a, double rMT,
                            double V0, int lMax, int numPoints);

    // ========================================================================
    //  16. DOS & Fermi Surface
    // ========================================================================
    static std::vector<double> computeDOSFromBands(const BandResult& bands, double broadening, int numBins);
    static std::vector<FermiSurfacePoint> computeFermiSurface2D(
        std::function<double(double, double)> Ek, double EF, int numPoints);
    void exportDOSFromBandsCSV(const std::string& filename, const BandResult& bands,
                                double broadening, int numBins);
    void exportFermiSurface2DCSV(const std::string& filename, double E0, double t,
                                   double a, double EF, int numPoints);

    // ========================================================================
    //  17. Effective Mass & k·p Theory
    // ========================================================================
    static double effectiveMassFromBand(double d2Edk2);
    static double kpBandEnergy(double k, double Eg, double P, double mStar);
    static std::pair<double, double> kpTwoBand(double k, double Eg, double P);
    static double cyclotronMass(double EF, double dAdE);
    void exportKPTheoryCSV(const std::string& filename, double Eg, double P, int numPoints);

    // ========================================================================
    //  18. Pseudopotential Method
    // ========================================================================
    static double pseudopotentialFormFactor(double q, double V_s, double r_c);
    static BandResult computePseudopotentialBands(double a, const std::vector<double>& V_G,
                                                    int numBands, int numPoints);
    void exportPseudopotentialCSV(const std::string& filename, double a,
                                   const std::vector<double>& V_G, int numBands, int numPoints);

    // ========================================================================
    //  19. Wannier Functions
    // ========================================================================
    static double wannierFunction1D(double x, double a, int nBands,
                                     const std::vector<double>& blochCoeffs);
    static double wannierSpread(double a, int nBands);
    void exportWannierCSV(const std::string& filename, double a, int nBands, int numPoints);

    // ========================================================================
    //  20. Band Structure Visualization
    // ========================================================================
    static BandResult computeModelBandStructure(const std::string& model, double a,
                                                  double param1, double param2, int numBands, int numPoints);
    void exportBandStructureCSV(const std::string& filename, const std::string& model,
                                 double a, double param1, double param2, int numBands, int numPoints);

    // ========================================================================
    //  21. Intrinsic Semiconductors
    // ========================================================================
    static double intrinsicCarrierConcentration(double Eg, double T, double mStar_e, double mStar_h);
    static double intrinsicFermiLevel(double Eg, double T, double mStar_e, double mStar_h);
    static double conductivityIntrinsic(double ni, double mu_e, double mu_h);
    void exportIntrinsicSemiCSV(const std::string& filename, double Eg, double mStar_e,
                                  double mStar_h, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  22. Extrinsic Semiconductors
    // ========================================================================
    static double donorIonizationFraction(double Nd, double Ed, double T, double mStar);
    static double electronConcentration(double Nd, double Na, double Eg, double T,
                                          double mStar_e, double mStar_h);
    static double holeConcentration(double n, double ni);
    static double extrinsicFermiLevel(double Eg, double T, double Nd, double Na,
                                        double mStar_e, double mStar_h);
    void exportExtrinsicSemiCSV(const std::string& filename, double Eg, double Nd, double Na,
                                  double mStar_e, double mStar_h, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  23. p-n Junction
    // ========================================================================
    static double builtInVoltage(double Na, double Nd, double ni, double T);
    static double depletionWidth(double Na, double Nd, double epsilon_r, double V_bi, double V_app);
    static PNJunctionResult computePNJunction(double Na, double Nd, double ni, double epsilon_r,
                                               double T, double V_app, int numPoints);
    static double diodeCurrent(double Is, double V, double T, double n_ideal);
    void exportPNJunctionCSV(const std::string& filename, double Na, double Nd, double ni,
                              double epsilon_r, double T, double V_app, int numPoints);

    // ========================================================================
    //  24. Semiconductor Optical Properties
    // ========================================================================
    static double directAbsorptionCoeff(double E_photon, double Eg, double mStar_r);
    static double indirectAbsorptionCoeff(double E_photon, double Eg, double E_phonon, double T);
    static double luminescenceSpectrum(double E_photon, double Eg, double T);
    static double jointDOS(double E, double Eg, double mStar_r);
    void exportOpticalPropertiesCSV(const std::string& filename, double Eg, double mStar_r,
                                     double T, double Emin, double Emax, int numPoints);

    // ========================================================================
    //  25. Quantum Wells & Heterostructures
    // ========================================================================
    double quantumWellEnergy(int n, double Lz, double mStar_w);
    double quantumWellWavefunction(int n, double z, double Lz);
    std::pair<double, double> finiteQWBoundState(double V0, double Lz, double mStar_w, double mStar_b, int n);
    static double superlatticeMiniBand(double t, double k, double d);
    void exportQuantumWellCSV(const std::string& filename, double Lz, double mStar_w,
                                int maxN, int numPoints);

    // ========================================================================
    //  26. Hall Effect & Magnetoresistance
    // ========================================================================
    static double hallCoefficient(double n_carrier, double charge);
    static double hallVoltage(double I, double B, double n_carrier, double t);
    static double hallMobility(double RH, double sigma);
    static double magnetoresistanceRatio(double mu, double B);
    static double cyclotronFrequency(double B, double mStar);
    void exportHallEffectCSV(const std::string& filename, double n_carrier, double mu,
                              double Bmin, double Bmax, int numPoints);

    // ========================================================================
    //  27. Boltzmann Transport
    // ========================================================================
    static double drudeCondutivity(double n, double tau, double mStar);
    static double drudeThermalConductivity(double n, double tau, double T, double mStar);
    static double wiedemannFranzRatio(double T);
    static double lorenzNumber();
    static double relaxationTimeApprox(double mu, double mStar);
    void exportBoltzmannTransportCSV(const std::string& filename, double n, double mStar,
                                       double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  28. Thermoelectric Effects
    // ========================================================================
    static double seebeckCoefficient(double T, double EF, double mStar);
    static double peltierCoefficient(double S, double T);
    static double thomsonCoefficient(double S, double T, double dSdT);
    static double thermoelectricZT(double S, double sigma, double kappa, double T);
    static double maxCoolingDeltaT(double ZT, double Tc);
    void exportThermoelectricCSV(const std::string& filename, double EF, double mStar,
                                   double sigma, double kappa, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  29. Diamagnetism & Paramagnetism
    // ========================================================================
    static double langevinFunction(double x);
    static double brillouinFunction(double x, double J);
    static double curieParaMagnetization(double B, double T, double J, double g, int N);
    static double curieParaSusceptibility(double T, double C);
    static double landauDiamagneticSusceptibility(double n, double mStar);
    static double pauliParamagneticSusceptibility(double dosEF);
    void exportParamagnetismCSV(const std::string& filename, double J, double g, int N,
                                  double Bmax, double T, int numPoints);

    // ========================================================================
    //  30. Ferromagnetism
    // ========================================================================
    static double weissFieldMagnetization(double M, double n, double J, double g,
                                            double T, double lambda_mf);
    static double curieWeissTemperature(double n, double J, double g, double lambda_mf);
    static MagnetizationResult solveMeanFieldFerromagnet(double n, double J, double g,
                                                          double lambda_mf, double Tmin,
                                                          double Tmax, int numPoints);
    static double spinWaveDispersion(double k, double a, double J, double S);
    void exportFerromagnetismCSV(const std::string& filename, double n, double J, double g,
                                   double lambda_mf, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  31. Heisenberg Model
    // ========================================================================
    static double heisenbergExchangeEnergy(double J, double S1_dot_S2);
    static double heisenbergGroundStateEnergy1D(int N, double J, double S);
    static double spinWaveEnergy(double k, double a, double J, double S);
    static double blochT32Law(double T, double Tc, double M0);
    void exportHeisenbergCSV(const std::string& filename, double J, double S, double a,
                               int numPoints);

    // ========================================================================
    //  32. Ising Model (1D/2D)
    // ========================================================================
    static double ising1DPartitionFunction(int N, double J, double B, double T);
    static double ising1DMagnetization(double J, double B, double T);
    static double ising1DCorrelationLength(double J, double T);
    static double ising2DCriticalTemperature(double J);
    static double ising2DExactMagnetization(double T, double Tc);
    static IsingResult ising2DMonteCarloSimulation(int L, double J, double Tmin, double Tmax,
                                                     int numTemps, int mcSteps);
    void exportIsing1DCSV(const std::string& filename, double J, double Bmax,
                            double Tmin, double Tmax, int numPoints);
    void exportIsing2DCSV(const std::string& filename, int L, double J, double Tmin,
                            double Tmax, int numTemps, int mcSteps);

    // ========================================================================
    //  33. Antiferromagnetism
    // ========================================================================
    static double neelTemperature(double J, double z, double S);
    static double afmSusceptibilityAboveTN(double T, double TN, double C);
    static double afmSusceptibilityBelowTN_parallel(double T, double TN, double chi_TN);
    static double afmSusceptibilityBelowTN_perp(double chi_TN);
    static double afmSublatticeMagnetization(double T, double TN, double M0);
    void exportAntiferromagnetismCSV(const std::string& filename, double J, double z,
                                      double S, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  34. Ferrimagnetism
    // ========================================================================
    static double ferriNetMagnetization(double MA, double MB);
    static double ferriCompensationTemperature(double JA, double JB, double SA, double SB);
    static MagnetizationResult ferriMeanField(double JA, double JB, double JAB,
                                                double SA, double SB, int NA, int NB,
                                                double Tmin, double Tmax, int numPoints);
    void exportFerrimagnetismCSV(const std::string& filename, double JA, double JB, double JAB,
                                   double SA, double SB, int NA, int NB,
                                   double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  35. Magnetic Domains
    // ========================================================================
    static double domainWallWidth(double A_ex, double K_aniso);
    static double domainWallEnergy(double A_ex, double K_aniso);
    static double singleDomainCriticalRadius(double A_ex, double K_aniso, double Ms);
    static double hysteresisLoop(double H, double Hc, double Ms, double Mr);
    void exportMagneticDomainsCSV(const std::string& filename, double A_ex, double K_aniso,
                                    double Ms, double Hmax, int numPoints);

    // ========================================================================
    //  36. Spin-Orbit Coupling in Solids
    // ========================================================================
    static double rashbaEnergy(double k, double alpha_R, double mStar);
    static std::pair<double, double> rashbaSplit(double k, double alpha_R, double mStar);
    static double dresselhausEnergy(double k, double beta_D, double mStar);
    static double spinHallAngle(double sigma_SH, double sigma_xx);
    void exportSpinOrbitCSV(const std::string& filename, double alpha_R, double beta_D,
                              double mStar, int numPoints);

    // ========================================================================
    //  37. BCS Theory
    // ========================================================================
    static double bcsGap(double T, double Delta0, double Tc);
    static double bcsCriticalTemperature(double omegaD, double N_EF, double V_eff);
    static double bcsCoherenceLength(double vF, double Delta0);
    static double bcsSpecificHeatJump();
    static double bcsDOSQuasiparticle(double E, double Delta);
    static BCSResult computeBCSGapVsTemp(double omegaD, double N_EF, double V_eff, int numPoints);
    void exportBCSTheoryCSV(const std::string& filename, double omegaD, double N_EF,
                              double V_eff, int numPoints);

    // ========================================================================
    //  38. Ginzburg-Landau Theory
    // ========================================================================
    static double glCoherenceLength(double xi0, double T, double Tc);
    static double glPenetrationDepth(double lambda0, double T, double Tc);
    static double glKappaParameter(double lambda, double xi);
    static double glCriticalField(double Hc0, double T, double Tc);
    static double glUpperCriticalField(double Hc, double kappa);
    static double glLowerCriticalField(double Hc, double kappa);
    static GLResult computeGLProfile(double xi, double lambda, int numPoints);
    void exportGLTheoryCSV(const std::string& filename, double xi0, double lambda0,
                             double Tc, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  39. Type I & II Superconductors
    // ========================================================================
    static double thermodynamicCriticalField(double Hc0, double T, double Tc);
    static bool isTypeII(double kappa);
    static double fluxQuantum();
    static double vortexSpacing(double B);
    static double mixedStateB(double H, double Hc1, double Hc2);
    void exportSCTypesCSV(const std::string& filename, double Hc0, double kappa,
                            double Tc, int numPoints);

    // ========================================================================
    //  40. Josephson Effect
    // ========================================================================
    static double josephsonCriticalCurrent(double Delta, double Rn);
    static double josephsonDCCurrent(double Ic, double phi);
    static double josephsonACFrequency(double V);
    static double josephsonPlasmaFrequency(double Ic, double C);
    static double shapiroStepVoltage(int n, double freq);
    static double squidFluxResponse(double Ic, double flux, double L);
    void exportJosephsonCSV(const std::string& filename, double Ic, double Rn,
                              double Vmax, int numPoints);

    // ========================================================================
    //  41. Electron-Phonon Interaction
    // ========================================================================
    static double eliashbergCouplingLambda(double N_EF, double g2, double omegaPh);
    static double mcMillanTc(double omegaLog, double lambda, double muStar);
    static double electronPhononScatteringRate(double T, double thetaD);
    static double resistivityBlochGruneisen(double T, double thetaD, double rho0, double K);
    void exportElectronPhononCSV(const std::string& filename, double thetaD, double rho0,
                                   double K, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  42. Anderson Localization
    // ========================================================================
    static double andersonLocalizationLength1D(double mfp);
    static double andersonLocalizationLength2D(double mfp);
    static double ioffRegelCriterion(double kF, double mfp);
    static double thoulessEnergy(double D, double L);
    static double conductanceDimensionless(double g0, double L, double xi);
    void exportAndersonCSV(const std::string& filename, double mfp, double kF,
                             double Lmin, double Lmax, int numPoints);

    // ========================================================================
    //  43. Quantum Hall Effect (Extended)
    // ========================================================================
    static double landauLevel(int n, double B, double mStar);
    static double magneticLength(double B);
    static double fillingFactor(double n2D, double B);
    static double hallConductance(double nu);
    static double hallResistance(double nu);
    static double longitudinalResistance(double B, double n2D, double mStar, double broadening);
    static double fractionalFilling(int p, int q);
    static double laughlinGapEstimate(double B, double mStar);
    void exportQuantumHallCSV(const std::string& filename, double n2D, double mStar,
                                double Bmin, double Bmax, int numPoints);

    // ========================================================================
    //  44. Topological Insulators (Toy)
    // ========================================================================
    static double sshEnergy(double k, double v, double w);
    static std::pair<double, double> sshBands(double k, double v, double w);
    static double windingNumber(double v, double w);
    static int cherNumberTwoLevel(double m, double t);
    static BandResult computeSSHBands(double v, double w, int numPoints);
    void exportTopologicalCSV(const std::string& filename, double v, double w, int numPoints);

    // ========================================================================
    //  45. Graphene & Dirac Materials
    // ========================================================================
    static double grapheneDispersion(double kx, double ky, double t, double a);
    static std::pair<double, double> grapheneBands(double kx, double ky, double t, double a);
    static double grapheneDOS(double E, double t);
    static double grapheneFermiVelocity(double t, double a);
    static double grapheneLandauLevel(int n, double B);
    static double grapheneCyclotronMass(double EF, double vF);
    void exportGrapheneCSV(const std::string& filename, double t, double a, int numPoints);

    // ========================================================================
    //  46. Plasmons & Dielectric Function
    // ========================================================================
    static double plasmaFrequency(double n, double mStar, double epsilon_inf);
    static std::complex<double> drudePermittivity(double omega, double omega_p, double gamma);
    static double lindhardDielectric(double q, double omega, double kF, double vF, double mStar);
    static double thomasFermiScreeningLength(double n, double mStar);
    static double plasmonDispersion(double q, double omega_p, double vF);
    void exportPlasmonCSV(const std::string& filename, double n, double mStar,
                            double gamma, double omegaMax, int numPoints);

    // ========================================================================
    //  47. Excitons
    // ========================================================================
    static double excitonBindingEnergy(int n, double mStar_r, double epsilon_r);
    static double excitonBohrRadius(double mStar_r, double epsilon_r);
    static double excitonReducedMass(double mStar_e, double mStar_h);
    static double mottDensity(double aB_exc, int dim);
    static double excitonAbsorption(double E, double Eg, double Eb, double broadening);
    void exportExcitonCSV(const std::string& filename, double mStar_e, double mStar_h,
                            double epsilon_r, double Eg, int numPoints);

    // ========================================================================
    //  48. Polaritons
    // ========================================================================
    static std::pair<double, double> polaritonDispersion(double k, double omega_c, double omega_exc,
                                                          double g_coupling);
    static double rabiSplitting(double g_coupling);
    static double polaritonHopfieldCoeff(double detuning, double g_coupling);
    static double cavityPhotonEnergy(double k, double omega0, double mStar_ph);
    void exportPolaritonCSV(const std::string& filename, double omega_c, double omega_exc,
                              double g_coupling, int numPoints);

    // ========================================================================
    //  49. Amorphous Solids
    // ========================================================================
    static double radialDistributionFunction(double r, double r0, double sigma, double n0);
    static double tunnelingSystems(double T, double P0);
    static double amorphousSpecificHeat(double T, double alpha, double beta);
    static double amorphousThermalConductivity(double T, double kMin, double kPhonon);
    void exportAmorphousSolidsCSV(const std::string& filename, double r0, double sigma,
                                    double n0, double Tmin, double Tmax, int numPoints);

    // ========================================================================
    //  50. Mesoscopic Transport
    // ========================================================================
    static double landauerConductance(double T_transmission, int numChannels);
    static double conductanceQuantum();
    static double coulombBlockadeChargingEnergy(double C);
    static double coulombDiamondBoundary(double Vg, double Vsd, double Ec, int n);
    static TransportResult computeLandauerTransport(double barrierHeight, double barrierWidth,
                                                      double mStar, double Emin, double Emax,
                                                      int numPoints);
    static double abFlux(double B, double area);
    static double aharonovBohmOscillationPeriod(double area);
    void exportMesoscopicCSV(const std::string& filename, double barrierH, double barrierW,
                               double mStar, double Emin, double Emax, int numPoints);

private:
    std::string m_materialName;
    double m_a;       // lattice constant (m)
    double m_mStar;   // effective mass (kg)
};
