// ============================================================================
//  NumericalSolverSS.cpp  —  Numerical solvers implementation
// ============================================================================

#include "NumericalSolverSS.h"
#include "SolidStateConstants.h"

#include <Eigen/Dense>

#include <cmath>
#include <fstream>
#include <algorithm>
#include <random>
#include <cassert>

using namespace std;
using namespace SSConst;

void NumericalSolverSS::solveBandStructureFDM(
    double mStar, double xMin, double xMax, int numPoints,
    const vector<double>& potential, int numEigenstates,
    const string& outputFilename)
{
    if (numPoints < 3) return;
    int N = numPoints - 2; // interior points (Dirichlet BCs)
    double dx = (xMax - xMin) / (numPoints - 1);
    double coeff = hbar * hbar / (2.0 * mStar * dx * dx);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        double Vi = (i + 1 < (int)potential.size()) ? potential[i + 1] : 0.0;
        H(i, i) = 2.0 * coeff + Vi;
        if (i > 0)     H(i, i - 1) = -coeff;
        if (i < N - 1) H(i, i + 1) = -coeff;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
    if (solver.info() != Eigen::Success) return;

    int nStates = min(numEigenstates, N);
    ofstream out(outputFilename);
    out << "x";
    for (int n = 0; n < nStates; ++n)
        out << ",E" << n << ",psi" << n;
    out << "\n";

    for (int i = 0; i < numPoints; ++i) {
        double x = xMin + i * dx;
        out << x;
        for (int n = 0; n < nStates; ++n) {
            double E = solver.eigenvalues()(n);
            double psi = 0.0;
            if (i > 0 && i < numPoints - 1)
                psi = solver.eigenvectors()(i - 1, n);
            out << "," << E << "," << psi;
        }
        out << "\n";
    }
}

void NumericalSolverSS::timeEvolveCrankNicolson(
    double mStar, double xMin, double xMax, int numPoints,
    const vector<double>& V, const vector<complex<double>>& psi0,
    double dt, int numSteps, const string& outCsv, int snapshotEvery)
{
    if (numPoints < 3) return;
    int N = numPoints;
    double dx = (xMax - xMin) / (N - 1);
    double r = hbar * dt / (4.0 * mStar * dx * dx);

    // Tridiagonal Crank-Nicolson: (I + i*H*dt/2) psi^{n+1} = (I - i*H*dt/2) psi^n
    // Build tridiagonal coefficients
    vector<complex<double>> psi(psi0.begin(), psi0.end());
    psi.resize(N, {0.0, 0.0});

    vector<complex<double>> diagL(N), offL(N), diagR(N), offR(N);
    complex<double> ir(0.0, r);
    for (int i = 0; i < N; ++i) {
        double Vi = (i < (int)V.size()) ? V[i] : 0.0;
        complex<double> vTerm(0.0, dt * Vi / (2.0 * hbar));
        diagL[i] = 1.0 + 2.0 * ir + vTerm;
        diagR[i] = 1.0 - 2.0 * ir - vTerm;
    }
    for (int i = 0; i < N; ++i) {
        offL[i] = -ir;
        offR[i] = ir;
    }

    // Tridiagonal solver (Thomas algorithm)
    auto solveTridiag = [&](vector<complex<double>>& a, vector<complex<double>>& b,
                            vector<complex<double>>& c, vector<complex<double>>& d, int n) {
        // a=lower, b=diag, c=upper, d=rhs; result in d
        vector<complex<double>> cp(n), dp(n);
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        for (int i = 1; i < n; ++i) {
            complex<double> m = b[i] - a[i] * cp[i - 1];
            cp[i] = c[i] / m;
            dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
        }
        d[n - 1] = dp[n - 1];
        for (int i = n - 2; i >= 0; --i)
            d[i] = dp[i] - cp[i] * d[i + 1];
    };

    ofstream out(outCsv);
    out << "step,x,Re_psi,Im_psi,prob\n";

    auto writeSnapshot = [&](int step) {
        for (int i = 0; i < N; ++i) {
            double x = xMin + i * dx;
            out << step << "," << x << "," << psi[i].real() << ","
                << psi[i].imag() << "," << norm(psi[i]) << "\n";
        }
    };
    writeSnapshot(0);

    for (int s = 1; s <= numSteps; ++s) {
        // Build RHS: (I - i*H*dt/2) * psi
        vector<complex<double>> rhs(N);
        for (int i = 0; i < N; ++i) {
            rhs[i] = diagR[i] * psi[i];
            if (i > 0)     rhs[i] += offR[i] * psi[i - 1];
            if (i < N - 1) rhs[i] += offR[i] * psi[i + 1];
        }
        // Dirichlet BCs
        rhs[0] = {0.0, 0.0};
        rhs[N - 1] = {0.0, 0.0};

        vector<complex<double>> lower(N), diag(diagL), upper(N);
        for (int i = 0; i < N; ++i) { lower[i] = offL[i]; upper[i] = offL[i]; }
        solveTridiag(lower, diag, upper, rhs, N);
        psi = rhs;
        psi[0] = {0.0, 0.0};
        psi[N - 1] = {0.0, 0.0};

        if (s % snapshotEvery == 0)
            writeSnapshot(s);
    }
}

vector<complex<double>> NumericalSolverSS::makeGaussianWavepacket(
    int numPoints, double xMin, double xMax,
    double x0, double sigma, double k0)
{
    vector<complex<double>> psi(numPoints);
    double dx = (xMax - xMin) / (numPoints - 1);
    double norm = 0.0;
    for (int i = 0; i < numPoints; ++i) {
        double x = xMin + i * dx;
        double env = exp(-0.5 * (x - x0) * (x - x0) / (sigma * sigma));
        psi[i] = env * complex<double>(cos(k0 * x), sin(k0 * x));
        norm += abs(psi[i]) * abs(psi[i]) * dx;
    }
    double invNorm = 1.0 / sqrt(norm);
    for (auto& p : psi) p *= invNorm;
    return psi;
}

void NumericalSolverSS::solvePoissonSchrodinger(
    double mStar, double epsilon_r, double xMin, double xMax,
    int numPoints, const vector<double>& dopingProfile,
    int maxIterations, double convergenceTol,
    const string& outputFilename)
{
    if (numPoints < 4) return;
    int N = numPoints;
    double dx = (xMax - xMin) / (N - 1);
    double eps = epsilon0 * epsilon_r;

    vector<double> V(N, 0.0);       // electrostatic potential energy
    vector<double> density(N, 0.0);  // electron density
    vector<double> Vold(N, 0.0);

    for (int iter = 0; iter < maxIterations; ++iter) {
        // --- Solve Schrodinger with current potential ---
        int Ni = N - 2;
        double coeff = hbar * hbar / (2.0 * mStar * dx * dx);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(Ni, Ni);
        for (int i = 0; i < Ni; ++i) {
            H(i, i) = 2.0 * coeff + V[i + 1];
            if (i > 0)      H(i, i - 1) = -coeff;
            if (i < Ni - 1) H(i, i + 1) = -coeff;
        }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
        if (solver.info() != Eigen::Success) return;

        // Compute electron density from lowest states
        fill(density.begin(), density.end(), 0.0);
        int nOccupied = min(3, Ni);
        for (int n = 0; n < nOccupied; ++n) {
            for (int i = 0; i < Ni; ++i) {
                double psi = solver.eigenvectors()(i, n);
                density[i + 1] += psi * psi / dx; // normalized
            }
        }

        // --- Solve Poisson: d^2V/dx^2 = -e^2*(N_d - n)/(eps) ---
        Vold = V;
        vector<double> rhs(N, 0.0);
        for (int i = 1; i < N - 1; ++i) {
            double Nd = (i < (int)dopingProfile.size()) ? dopingProfile[i] : 0.0;
            rhs[i] = e_charge * e_charge * (Nd - density[i]) / eps;
        }
        // Simple iterative Poisson solve (Jacobi)
        for (int pIter = 0; pIter < 200; ++pIter) {
            vector<double> Vnew(N, 0.0);
            for (int i = 1; i < N - 1; ++i)
                Vnew[i] = 0.5 * (V[i - 1] + V[i + 1] + dx * dx * rhs[i]);
            V = Vnew;
        }

        // Mix old and new for stability
        double alpha = 0.3;
        for (int i = 0; i < N; ++i)
            V[i] = alpha * V[i] + (1.0 - alpha) * Vold[i];

        // Check convergence
        double maxDiff = 0.0;
        for (int i = 0; i < N; ++i)
            maxDiff = max(maxDiff, abs(V[i] - Vold[i]));
        if (maxDiff < convergenceTol) break;
    }

    // Write output
    ofstream out(outputFilename);
    out << "x,V,density\n";
    for (int i = 0; i < N; ++i) {
        double x = xMin + i * dx;
        out << x << "," << V[i] << "," << density[i] << "\n";
    }
}

double NumericalSolverSS::transferMatrixTransmission(
    double E, const vector<double>& widths,
    const vector<double>& heights, double mStar)
{
    // 2x2 transfer matrix for multilayer potential
    double M[2][2] = {{1,0},{0,1}};
    double k0 = sqrt(2.0 * mStar * E) / hbar;

    for (size_t i = 0; i < widths.size(); ++i) {
        double V = heights[i];
        double w = widths[i];
        if (E > V) {
            double k = sqrt(2.0 * mStar * (E - V)) / hbar;
            double c = cos(k * w), s = sin(k * w);
            double r = k0 / k;
            double T[2][2] = {
                {c, s / k},
                {-k * s, c}
            };
            double Mnew[2][2];
            Mnew[0][0] = T[0][0]*M[0][0] + T[0][1]*M[1][0];
            Mnew[0][1] = T[0][0]*M[0][1] + T[0][1]*M[1][1];
            Mnew[1][0] = T[1][0]*M[0][0] + T[1][1]*M[1][0];
            Mnew[1][1] = T[1][0]*M[0][1] + T[1][1]*M[1][1];
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    M[a][b] = Mnew[a][b];
        } else {
            double kappa = sqrt(2.0 * mStar * (V - E)) / hbar;
            double ch = cosh(kappa * w), sh = sinh(kappa * w);
            double T[2][2] = {
                {ch, sh / kappa},
                {kappa * sh, ch}
            };
            double Mnew[2][2];
            Mnew[0][0] = T[0][0]*M[0][0] + T[0][1]*M[1][0];
            Mnew[0][1] = T[0][0]*M[0][1] + T[0][1]*M[1][1];
            Mnew[1][0] = T[1][0]*M[0][0] + T[1][1]*M[1][0];
            Mnew[1][1] = T[1][0]*M[0][1] + T[1][1]*M[1][1];
            for (int a = 0; a < 2; ++a)
                for (int b = 0; b < 2; ++b)
                    M[a][b] = Mnew[a][b];
        }
    }
    double denom = M[0][0]*M[0][0] + (M[0][1]*k0)*(M[0][1]*k0);
    if (denom < 1e-30) return 0.0;
    return 4.0 * k0 * k0 / (denom + (M[1][0]/k0)*(M[1][0]/k0) + M[1][1]*M[1][1]);
}

double NumericalSolverSS::monteCarloIntegrate(
    function<double(double)> f, double a, double b, int numSamples)
{
    mt19937 rng(42);
    uniform_real_distribution<double> dist(a, b);
    double sum = 0.0;
    for (int i = 0; i < numSamples; ++i)
        sum += f(dist(rng));
    return (b - a) * sum / numSamples;
}
