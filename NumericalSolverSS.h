#pragma once
// ============================================================================
//  NumericalSolverSS.h  —  Numerical methods for solid-state physics
// ============================================================================

#include <vector>
#include <complex>
#include <string>
#include <functional>

class NumericalSolverSS {
public:
    // ── Finite-Difference Band Solver ──────────────────────────────────────
    static void solveBandStructureFDM(
        double mStar, double xMin, double xMax, int numPoints,
        const std::vector<double>& potential, int numEigenstates,
        const std::string& outputFilename);

    // ── Crank-Nicolson Time Evolution ──────────────────────────────────────
    static void timeEvolveCrankNicolson(
        double mStar, double xMin, double xMax, int numPoints,
        const std::vector<double>& V,
        const std::vector<std::complex<double>>& psi0,
        double dt, int numSteps, const std::string& outCsv,
        int snapshotEvery = 10);

    // ── Gaussian Initial State ─────────────────────────────────────────────
    static std::vector<std::complex<double>> makeGaussianWavepacket(
        int numPoints, double xMin, double xMax,
        double x0, double sigma, double k0);

    // ── Self-Consistent Poisson-Schrodinger ────────────────────────────────
    static void solvePoissonSchrodinger(
        double mStar, double epsilon_r, double xMin, double xMax,
        int numPoints, const std::vector<double>& dopingProfile,
        int maxIterations, double convergence,
        const std::string& outputFilename);

    // ── Transfer Matrix Method ─────────────────────────────────────────────
    static double transferMatrixTransmission(
        double E, const std::vector<double>& widths,
        const std::vector<double>& heights, double mStar);

    // ── Monte Carlo Integration ────────────────────────────────────────────
    static double monteCarloIntegrate(
        std::function<double(double)> f, double a, double b, int numSamples);
};
