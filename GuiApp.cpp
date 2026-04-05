// ============================================================================
//  GuiApp.cpp  —  ImGui + ImPlot GUI for SolidStateCore (50 simulations)
// ============================================================================

#include "GuiApp.h"
#include "SolidStateSystem.h"
#include "SolidStateConstants.h"
#include "NumericalSolverSS.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include <GLFW/glfw3.h>
#include <cmath>
#include <sstream>
#include <vector>

// ============================================================================
//  Plot helpers
// ============================================================================

void GuiApp::clearPlot() {
    m_plotX.clear();
    m_plotY.clear();
    m_plotLabels.clear();
}

void GuiApp::addCurve(const std::string& label,
                      const std::vector<double>& x,
                      const std::vector<double>& y) {
    m_plotLabels.push_back(label);
    m_plotX.push_back(x);
    m_plotY.push_back(y);
}

// ============================================================================
//  Init / Run / Shutdown
// ============================================================================

GuiApp::GuiApp() {}
GuiApp::~GuiApp() { shutdown(); }

bool GuiApp::init(int width, int height, const char* title) {
    if (!glfwInit()) return false;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    m_window = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!m_window) { glfwTerminate(); return false; }
    glfwMakeContextCurrent(m_window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(m_window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    return true;
}

void GuiApp::run() {
    static const char* names[] = {
        " 1  Bravais Lattices",
        " 2  Reciprocal Lattice",
        " 3  X-Ray Diffraction",
        " 4  Monatomic Lattice Vibrations",
        " 5  Diatomic Lattice Vibrations",
        " 6  Debye Model",
        " 7  Einstein Model",
        " 8  Phonon DOS",
        " 9  Thermal Expansion",
        "10  Lattice Thermal Conductivity",
        "11  Free Electron Model",
        "12  Nearly Free Electron",
        "13  Tight-Binding (Extended)",
        "14  Kronig-Penney (Extended)",
        "15  APW / OPW Methods",
        "16  DOS & Fermi Surface",
        "17  Effective Mass & k.p",
        "18  Pseudopotential Method",
        "19  Wannier Functions",
        "20  Band Structure Visualization",
        "21  Intrinsic Semiconductors",
        "22  Extrinsic Semiconductors",
        "23  p-n Junction",
        "24  Semiconductor Optics",
        "25  Quantum Wells & Heterostructures",
        "26  Hall Effect & Magnetoresistance",
        "27  Boltzmann Transport",
        "28  Thermoelectric Effects",
        "29  Dia- & Paramagnetism",
        "30  Ferromagnetism",
        "31  Heisenberg Model",
        "32  Ising Model (1D/2D)",
        "33  Antiferromagnetism",
        "34  Ferrimagnetism",
        "35  Magnetic Domains",
        "36  Spin-Orbit Coupling in Solids",
        "37  BCS Theory",
        "38  Ginzburg-Landau Theory",
        "39  Type I & II Superconductors",
        "40  Josephson Effect",
        "41  Electron-Phonon Interaction",
        "42  Anderson Localization",
        "43  Quantum Hall Effect (Extended)",
        "44  Topological Insulators (Toy)",
        "45  Graphene & Dirac Materials",
        "46  Plasmons & Dielectric Function",
        "47  Excitons",
        "48  Polaritons",
        "49  Amorphous Solids",
        "50  Mesoscopic Transport"
    };

    while (!glfwWindowShouldClose(m_window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // ── Left sidebar ───────────────────────────────────────────────
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(320, (float)ImGui::GetIO().DisplaySize.y));
        ImGui::Begin("Simulations", nullptr,
                     ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse);
        for (int i = 0; i < 50; ++i) {
            bool selected = (m_selectedSim == i);
            if (ImGui::Selectable(names[i], selected))
                m_selectedSim = i;
        }
        ImGui::End();

        // ── Center panel ───────────────────────────────────────────────
        ImGui::SetNextWindowPos(ImVec2(320, 0));
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetIO().DisplaySize.x - 320,
                                         ImGui::GetIO().DisplaySize.y));
        ImGui::Begin("##MainPanel", nullptr,
                     ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                     ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar);

        switch (m_selectedSim) {
            case  0: renderSim01_BravaisLattices();       break;
            case  1: renderSim02_ReciprocalLattice();     break;
            case  2: renderSim03_XRayDiffraction();       break;
            case  3: renderSim04_MonatomicLattice();      break;
            case  4: renderSim05_DiatomicLattice();       break;
            case  5: renderSim06_DebyeModel();            break;
            case  6: renderSim07_EinsteinModel();         break;
            case  7: renderSim08_PhononDOS();             break;
            case  8: renderSim09_ThermalExpansion();       break;
            case  9: renderSim10_LatticeThermalConductivity(); break;
            case 10: renderSim11_FreeElectronModel();     break;
            case 11: renderSim12_NearlyFreeElectron();    break;
            case 12: renderSim13_TightBinding();          break;
            case 13: renderSim14_KronigPenney();          break;
            case 14: renderSim15_APW_OPW();               break;
            case 15: renderSim16_DOS_FermiSurface();      break;
            case 16: renderSim17_EffectiveMass_KP();      break;
            case 17: renderSim18_Pseudopotential();       break;
            case 18: renderSim19_WannierFunctions();      break;
            case 19: renderSim20_BandStructureVis();      break;
            case 20: renderSim21_IntrinsicSemiconductors(); break;
            case 21: renderSim22_ExtrinsicSemiconductors(); break;
            case 22: renderSim23_PNJunction();            break;
            case 23: renderSim24_SemiOpticalProperties(); break;
            case 24: renderSim25_QuantumWells();          break;
            case 25: renderSim26_HallEffect();            break;
            case 26: renderSim27_BoltzmannTransport();    break;
            case 27: renderSim28_ThermoelectricEffects(); break;
            case 28: renderSim29_DiaParamagnetism();      break;
            case 29: renderSim30_Ferromagnetism();        break;
            case 30: renderSim31_HeisenbergModel();       break;
            case 31: renderSim32_IsingModel();            break;
            case 32: renderSim33_Antiferromagnetism();    break;
            case 33: renderSim34_Ferrimagnetism();        break;
            case 34: renderSim35_MagneticDomains();       break;
            case 35: renderSim36_SpinOrbitCoupling();     break;
            case 36: renderSim37_BCSTheory();             break;
            case 37: renderSim38_GinzburgLandau();        break;
            case 38: renderSim39_TypeI_II_SC();           break;
            case 39: renderSim40_JosephsonEffect();       break;
            case 40: renderSim41_ElectronPhonon();        break;
            case 41: renderSim42_AndersonLocalization();  break;
            case 42: renderSim43_QuantumHallEffect();     break;
            case 43: renderSim44_TopologicalInsulators(); break;
            case 44: renderSim45_GrapheneDirac();         break;
            case 45: renderSim46_PlasmonsDielectric();    break;
            case 46: renderSim47_Excitons();              break;
            case 47: renderSim48_Polaritons();            break;
            case 48: renderSim49_AmorphousSolids();       break;
            case 49: renderSim50_MesoscopicTransport();   break;
        }

        // ── Plot area ──────────────────────────────────────────────────
        if (!m_plotLabels.empty()) {
            if (ImPlot::BeginPlot("##Plot", ImVec2(-1, 350))) {
                for (size_t c = 0; c < m_plotLabels.size(); ++c) {
                    if (!m_plotX[c].empty())
                        ImPlot::PlotLine(m_plotLabels[c].c_str(),
                                         m_plotX[c].data(), m_plotY[c].data(),
                                         (int)m_plotX[c].size());
                }
                ImPlot::EndPlot();
            }
        }

        // ── Result text ────────────────────────────────────────────────
        if (!m_resultText.empty()) {
            ImGui::Separator();
            ImGui::TextWrapped("%s", m_resultText.c_str());
        }

        ImGui::End();

        // ── Render ─────────────────────────────────────────────────────
        ImGui::Render();
        int dw, dh;
        glfwGetFramebufferSize(m_window, &dw, &dh);
        glViewport(0, 0, dw, dh);
        glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(m_window);
    }
}

void GuiApp::shutdown() {
    if (m_window) {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        glfwDestroyWindow(m_window);
        glfwTerminate();
        m_window = nullptr;
    }
}

// ============================================================================
//  50 Simulation Render Stubs
//  Each follows the pattern:
//    - CollapsingHeader("Theory: ...")
//    - Combo for mode (if needed)
//    - Parameter inputs
//    - Compute button  ->  clearPlot(), compute, addCurve(), build resultText
//    - Export CSV button
// ============================================================================

void GuiApp::renderSim01_BravaisLattices() {
    ImGui::Text("1. Bravais Lattices");
    if (ImGui::CollapsingHeader("Theory: Bravais Lattices"))
        ImGui::TextWrapped("A Bravais lattice is an infinite set of points generated by a set of discrete translation operations. In 3D there are 14 distinct Bravais lattices grouped into 7 crystal systems.");

    static int latticeType = 0;
    static double a = 5.0e-10;
    static const char* types[] = {"SC","BCC","FCC","Hexagonal","Tetragonal","Orthorhombic","Monoclinic"};
    ImGui::Combo("Lattice Type##Sim01", &latticeType, types, IM_ARRAYSIZE(types));
    ImGui::InputDouble("a (m)##Sim01", &a, 0, 0, "%.3e");

    if (ImGui::Button("Compute##Sim01")) {
        clearPlot();
        auto lat = SolidStateSystem::getBravaisLattice(latticeType, a);
        std::ostringstream oss;
        oss << "Lattice: " << lat.name << "\n"
            << "a1 = (" << lat.a1.x << ", " << lat.a1.y << ", " << lat.a1.z << ")\n"
            << "a2 = (" << lat.a2.x << ", " << lat.a2.y << ", " << lat.a2.z << ")\n"
            << "a3 = (" << lat.a3.x << ", " << lat.a3.y << ", " << lat.a3.z << ")\n"
            << "d(100) = " << SolidStateSystem::latticeSpacing(1,0,0,a) << " m\n"
            << "d(110) = " << SolidStateSystem::latticeSpacing(1,1,0,a) << " m\n"
            << "d(111) = " << SolidStateSystem::latticeSpacing(1,1,1,a) << " m";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim01")) {
        SolidStateSystem sys;
        sys.exportBravaisLatticeCSV("sim01_bravais.csv", latticeType, a);
    }
}

void GuiApp::renderSim02_ReciprocalLattice() {
    ImGui::Text("2. Reciprocal Lattice");
    if (ImGui::CollapsingHeader("Theory: Reciprocal Lattice"))
        ImGui::TextWrapped("The reciprocal lattice is the Fourier transform of the direct lattice. Reciprocal lattice vectors b_i satisfy a_i . b_j = 2*pi*delta_ij. The first Brillouin zone is the Wigner-Seitz cell of the reciprocal lattice.");

    static double a = 5.0e-10;
    static int maxIndex = 3;
    ImGui::InputDouble("a (m)##Sim02", &a, 0, 0, "%.3e");
    ImGui::SliderInt("Max Index##Sim02", &maxIndex, 1, 5);

    if (ImGui::Button("Compute##Sim02")) {
        clearPlot();
        auto lat = SolidStateSystem::getBravaisLattice(0, a);
        auto b1 = SolidStateSystem::computeReciprocalVector(lat.a1, lat.a2, lat.a3, 0);
        auto b2 = SolidStateSystem::computeReciprocalVector(lat.a1, lat.a2, lat.a3, 1);
        auto b3 = SolidStateSystem::computeReciprocalVector(lat.a1, lat.a2, lat.a3, 2);
        double bzBound = SolidStateSystem::brillouinZoneBoundary1D(a);
        std::ostringstream oss;
        oss << "b1 = (" << b1.x << ", " << b1.y << ", " << b1.z << ")\n"
            << "b2 = (" << b2.x << ", " << b2.y << ", " << b2.z << ")\n"
            << "b3 = (" << b3.x << ", " << b3.y << ", " << b3.z << ")\n"
            << "1st BZ boundary (1D) = +/-" << bzBound << " m^-1";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim02")) {
        SolidStateSystem sys;
        sys.exportReciprocalLatticeCSV("sim02_reciprocal.csv", a, maxIndex);
    }
}

void GuiApp::renderSim03_XRayDiffraction() {
    ImGui::Text("3. X-Ray Diffraction");
    if (ImGui::CollapsingHeader("Theory: X-Ray Diffraction"))
        ImGui::TextWrapped("Bragg's law: 2d*sin(theta) = n*lambda. X-ray diffraction probes crystal structure via constructive interference from lattice planes. The structure factor determines which reflections are allowed.");

    static double a = 5.43e-10;
    static double wavelength = 1.54e-10;
    static int latticeIdx = 2;
    static int maxHKL = 3;
    static const char* latticeNames[] = {"SC", "BCC", "FCC"};
    ImGui::InputDouble("a (m)##Sim03", &a, 0, 0, "%.3e");
    ImGui::InputDouble("Wavelength (m)##Sim03", &wavelength, 0, 0, "%.3e");
    ImGui::Combo("Lattice##Sim03", &latticeIdx, latticeNames, IM_ARRAYSIZE(latticeNames));
    ImGui::SliderInt("Max HKL##Sim03", &maxHKL, 1, 6);

    if (ImGui::Button("Compute##Sim03")) {
        clearPlot();
        std::string lt = latticeNames[latticeIdx];
        auto peaks = SolidStateSystem::computeDiffractionPattern(a, wavelength, lt, maxHKL);
        std::vector<double> xvals, yvals;
        for (auto& p : peaks) {
            xvals.push_back(p.twoTheta * 180.0 / SSConst::PI);
            yvals.push_back(p.intensity);
        }
        addCurve("Diffraction", xvals, yvals);
        std::ostringstream oss;
        oss << "Found " << peaks.size() << " peaks\n";
        for (auto& p : peaks)
            oss << "(" << p.h << p.k << p.l << ") 2theta=" << p.twoTheta*180.0/SSConst::PI << " deg  I=" << p.intensity << "\n";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim03")) {
        SolidStateSystem sys;
        sys.exportDiffractionCSV("sim03_xrd.csv", a, wavelength, latticeNames[latticeIdx], maxHKL);
    }
}

void GuiApp::renderSim04_MonatomicLattice() {
    ImGui::Text("4. Monatomic Lattice Vibrations");
    if (ImGui::CollapsingHeader("Theory: 1D Monatomic Chain"))
        ImGui::TextWrapped("omega(k) = 2*sqrt(K/M)*|sin(ka/2)|. The dispersion is periodic in k with period 2*pi/a. Group velocity vanishes at the zone boundary.");

    static double a = 3.0e-10;
    static double springK = 10.0;
    static double mass = 4.0e-26;
    static int numPts = 200;
    ImGui::InputDouble("a (m)##Sim04", &a, 0, 0, "%.3e");
    ImGui::InputDouble("Spring K (N/m)##Sim04", &springK);
    ImGui::InputDouble("Mass (kg)##Sim04", &mass, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim04", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim04")) {
        clearPlot();
        std::vector<double> kv, wv;
        double kmax = SSConst::PI / a;
        for (int i = 0; i < numPts; ++i) {
            double k = -kmax + 2.0 * kmax * i / (numPts - 1);
            kv.push_back(k);
            wv.push_back(SolidStateSystem::monatomicDispersion(k, a, springK, mass));
        }
        addCurve("omega(k)", kv, wv);
        double omMax = SolidStateSystem::debyeCutoffFrequency(springK, mass);
        std::ostringstream oss;
        oss << "omega_max = " << omMax << " rad/s";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim04")) {
        SolidStateSystem sys;
        sys.exportMonatomicDispersionCSV("sim04_monatomic.csv", a, springK, mass, numPts);
    }
}

void GuiApp::renderSim05_DiatomicLattice() {
    ImGui::Text("5. Diatomic Lattice Vibrations");
    if (ImGui::CollapsingHeader("Theory: 1D Diatomic Chain"))
        ImGui::TextWrapped("Two atoms per unit cell yield acoustic and optical branches. A band gap opens between them when m1 != m2. The optical branch has atoms moving out of phase.");

    static double a = 5.0e-10;
    static double springK = 10.0;
    static double m1 = 4.0e-26;
    static double m2 = 6.0e-26;
    static int numPts = 200;
    ImGui::InputDouble("a (m)##Sim05", &a, 0, 0, "%.3e");
    ImGui::InputDouble("Spring K (N/m)##Sim05", &springK);
    ImGui::InputDouble("m1 (kg)##Sim05", &m1, 0, 0, "%.3e");
    ImGui::InputDouble("m2 (kg)##Sim05", &m2, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim05", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim05")) {
        clearPlot();
        std::vector<double> kv, ac, op;
        double kmax = SSConst::PI / a;
        for (int i = 0; i < numPts; ++i) {
            double k = -kmax + 2.0 * kmax * i / (numPts - 1);
            kv.push_back(k);
            auto [wa, wo] = SolidStateSystem::diatomicDispersion(k, a, springK, m1, m2);
            ac.push_back(wa);
            op.push_back(wo);
        }
        addCurve("Acoustic", kv, ac);
        addCurve("Optical", kv, op);
        double gap = SolidStateSystem::diatomicBandGap(springK, m1, m2);
        std::ostringstream oss;
        oss << "Band gap = " << gap << " rad/s";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim05")) {
        SolidStateSystem sys;
        sys.exportDiatomicDispersionCSV("sim05_diatomic.csv", a, springK, m1, m2, numPts);
    }
}

void GuiApp::renderSim06_DebyeModel() {
    ImGui::Text("6. Debye Model");
    if (ImGui::CollapsingHeader("Theory: Debye Model"))
        ImGui::TextWrapped("The Debye model assumes a linear dispersion omega = v*k up to a cutoff frequency omega_D. It correctly predicts the T^3 low-temperature specific heat and approaches the Dulong-Petit limit at high T.");

    static double thetaD = 345.0;
    static int N = 1000;
    static double Tmin = 1.0, Tmax = 500.0;
    static int numPts = 200;
    ImGui::InputDouble("Theta_D (K)##Sim06", &thetaD);
    ImGui::InputInt("N atoms##Sim06", &N);
    ImGui::InputDouble("T min (K)##Sim06", &Tmin);
    ImGui::InputDouble("T max (K)##Sim06", &Tmax);
    ImGui::SliderInt("Points##Sim06", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim06")) {
        clearPlot();
        std::vector<double> Tv, Cv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            Tv.push_back(T);
            Cv.push_back(SolidStateSystem::debyeSpecificHeat(T, thetaD));
        }
        addCurve("Cv/3NkB (Debye)", Tv, Cv);
        std::ostringstream oss;
        oss << "Debye temperature = " << thetaD << " K\n"
            << "Cv(300K) = " << SolidStateSystem::debyeSpecificHeat(300.0, thetaD) << " *3NkB";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim06")) {
        SolidStateSystem sys;
        sys.exportDebyeModelCSV("sim06_debye.csv", thetaD, N, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim07_EinsteinModel() {
    ImGui::Text("7. Einstein Model");
    if (ImGui::CollapsingHeader("Theory: Einstein Model"))
        ImGui::TextWrapped("The Einstein model treats each atom as an independent quantum harmonic oscillator with the same frequency omega_E. It gives the correct high-T limit but falls off too rapidly at low T (exponential vs T^3).");

    static double thetaE = 200.0;
    static int N = 1000;
    static double Tmin = 1.0, Tmax = 500.0;
    static int numPts = 200;
    ImGui::InputDouble("Theta_E (K)##Sim07", &thetaE);
    ImGui::InputInt("N atoms##Sim07", &N);
    ImGui::InputDouble("T min (K)##Sim07", &Tmin);
    ImGui::InputDouble("T max (K)##Sim07", &Tmax);
    ImGui::SliderInt("Points##Sim07", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim07")) {
        clearPlot();
        std::vector<double> Tv, Cv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            Tv.push_back(T);
            Cv.push_back(SolidStateSystem::einsteinSpecificHeat(T, thetaE));
        }
        addCurve("Cv/3NkB (Einstein)", Tv, Cv);
        std::ostringstream oss;
        oss << "Einstein temperature = " << thetaE << " K\n"
            << "Cv(300K) = " << SolidStateSystem::einsteinSpecificHeat(300.0, thetaE) << " *3NkB";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim07")) {
        SolidStateSystem sys;
        sys.exportEinsteinModelCSV("sim07_einstein.csv", thetaE, N, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim08_PhononDOS() {
    ImGui::Text("8. Phonon Density of States");
    if (ImGui::CollapsingHeader("Theory: Phonon DOS"))
        ImGui::TextWrapped("The phonon density of states g(omega) counts the number of modes per unit frequency. In the Debye model g ~ omega^2 in 3D. Van Hove singularities appear where the group velocity vanishes.");

    static double omegaD = 5.0e13;
    static int numPts = 200;
    ImGui::InputDouble("omega_D (rad/s)##Sim08", &omegaD, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim08", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim08")) {
        clearPlot();
        std::vector<double> wv, gv;
        for (int i = 0; i < numPts; ++i) {
            double w = omegaD * i / (numPts - 1);
            wv.push_back(w);
            gv.push_back(SolidStateSystem::debyeDOS3D(w, omegaD));
        }
        addCurve("g(omega) Debye 3D", wv, gv);
        std::ostringstream oss;
        oss << "Debye cutoff omega_D = " << omegaD << " rad/s";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim08")) {
        SolidStateSystem sys;
        sys.exportPhononDOSCSV("sim08_phonon_dos.csv", omegaD, numPts);
    }
}

void GuiApp::renderSim09_ThermalExpansion() {
    ImGui::Text("9. Thermal Expansion");
    if (ImGui::CollapsingHeader("Theory: Thermal Expansion"))
        ImGui::TextWrapped("Thermal expansion arises from anharmonicity in the interatomic potential. The Gruneisen parameter gamma relates the thermal expansion coefficient to the specific heat: alpha = gamma*Cv/(B*V).");

    static double thetaD = 345.0;
    static double gamma_g = 2.0;
    static double V0 = 1.2e-29;
    static double B = 100.0e9;
    static double Tmin = 1.0, Tmax = 500.0;
    static int numPts = 200;
    ImGui::InputDouble("Theta_D (K)##Sim09", &thetaD);
    ImGui::InputDouble("Gruneisen gamma##Sim09", &gamma_g);
    ImGui::InputDouble("V0 (m^3)##Sim09", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("Bulk modulus B (Pa)##Sim09", &B, 0, 0, "%.3e");
    ImGui::InputDouble("T min (K)##Sim09", &Tmin);
    ImGui::InputDouble("T max (K)##Sim09", &Tmax);
    ImGui::SliderInt("Points##Sim09", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim09")) {
        clearPlot();
        std::vector<double> Tv, av;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            double Cv = SolidStateSystem::debyeSpecificHeat(T, thetaD) * 3.0 * SSConst::kB;
            double alpha = SolidStateSystem::thermalExpansionCoeff(gamma_g, Cv, B, V0);
            Tv.push_back(T);
            av.push_back(alpha);
        }
        addCurve("alpha(T)", Tv, av);
        double Cv300 = SolidStateSystem::debyeSpecificHeat(300.0, thetaD) * 3.0 * SSConst::kB;
        double alpha300 = SolidStateSystem::thermalExpansionCoeff(gamma_g, Cv300, B, V0);
        std::ostringstream oss;
        oss << "Gruneisen gamma = " << gamma_g << "\n"
            << "alpha(300K) = " << alpha300 << " K^-1";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim09")) {
        SolidStateSystem sys;
        sys.exportThermalExpansionCSV("sim09_thermal_expansion.csv", thetaD, gamma_g, V0, B, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim10_LatticeThermalConductivity() {
    ImGui::Text("10. Lattice Thermal Conductivity");
    if (ImGui::CollapsingHeader("Theory: Lattice Thermal Conductivity"))
        ImGui::TextWrapped("kappa = (1/3)*Cv*v*l where Cv is specific heat, v is sound velocity, and l is phonon mean free path. At low T, boundary scattering dominates (kappa ~ T^3). At high T, Umklapp scattering gives kappa ~ 1/T.");

    static double thetaD = 345.0;
    static double v_sound = 5000.0;
    static double Tmin = 1.0, Tmax = 500.0;
    static int numPts = 200;
    ImGui::InputDouble("Theta_D (K)##Sim10", &thetaD);
    ImGui::InputDouble("v_sound (m/s)##Sim10", &v_sound);
    ImGui::InputDouble("T min (K)##Sim10", &Tmin);
    ImGui::InputDouble("T max (K)##Sim10", &Tmax);
    ImGui::SliderInt("Points##Sim10", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim10")) {
        clearPlot();
        std::vector<double> Tv, kv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            double Cv = SolidStateSystem::debyeSpecificHeat(T, thetaD) * 3.0 * SSConst::kB;
            double mfp = SolidStateSystem::phononMeanFreePath(T, thetaD);
            double kappa = SolidStateSystem::kineticTheoryThermalConductivity(Cv, v_sound, mfp);
            Tv.push_back(T);
            kv.push_back(kappa);
        }
        addCurve("kappa(T)", Tv, kv);
        std::ostringstream oss;
        oss << "Theta_D = " << thetaD << " K\n"
            << "v_sound = " << v_sound << " m/s";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim10")) {
        SolidStateSystem sys;
        sys.exportThermalConductivityCSV("sim10_thermal_cond.csv", thetaD, v_sound, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim11_FreeElectronModel() {
    ImGui::Text("11. Free Electron Model");
    if (ImGui::CollapsingHeader("Theory: Free Electron Model"))
        ImGui::TextWrapped("Electrons in a metal treated as free particles in a box. E(k) = hbar^2*k^2/(2m). The Fermi-Dirac distribution governs occupation. DOS ~ sqrt(E) in 3D. Electronic specific heat is linear in T.");

    static double n_density = 8.5e28;
    static double mStar = SSConst::me;
    static double Emax = 15.0 * SSConst::eV;
    static int numPts = 200;
    ImGui::InputDouble("n (m^-3)##Sim11", &n_density, 0, 0, "%.3e");
    ImGui::InputDouble("m* (kg)##Sim11", &mStar, 0, 0, "%.3e");
    ImGui::InputDouble("E_max (J)##Sim11", &Emax, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim11", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim11")) {
        clearPlot();
        double EF = SolidStateSystem::fermiEnergy(n_density, mStar);
        std::vector<double> Ev, dosv;
        for (int i = 0; i < numPts; ++i) {
            double E = Emax * (i + 1) / numPts;
            Ev.push_back(E / SSConst::eV);
            dosv.push_back(SolidStateSystem::freeElectronDOS3D(E, mStar));
        }
        addCurve("DOS(E)", Ev, dosv);
        double kF = SolidStateSystem::fermiWavevector(EF, mStar);
        double vF = SolidStateSystem::fermiVelocity(EF, mStar);
        double TF = SolidStateSystem::fermiTemperature(EF);
        std::ostringstream oss;
        oss << "E_F = " << EF/SSConst::eV << " eV\n"
            << "k_F = " << kF << " m^-1\n"
            << "v_F = " << vF << " m/s\n"
            << "T_F = " << TF << " K";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim11")) {
        SolidStateSystem sys;
        sys.exportFreeElectronCSV("sim11_free_electron.csv", n_density, mStar, Emax, numPts);
    }
}

void GuiApp::renderSim12_NearlyFreeElectron() {
    ImGui::Text("12. Nearly Free Electron Model");
    if (ImGui::CollapsingHeader("Theory: Nearly Free Electron"))
        ImGui::TextWrapped("A weak periodic potential opens band gaps at Brillouin zone boundaries. Near a zone boundary, degenerate perturbation theory gives E = E_avg +/- sqrt((E1-E2)^2/4 + |V_G|^2). The band gap = 2|V_G|.");

    static double a = 5.0e-10;
    static double V_G = 1.0 * SSConst::eV;
    static double mStar = SSConst::me;
    static int numBands = 3;
    static int numPts = 200;
    ImGui::InputDouble("a (m)##Sim12", &a, 0, 0, "%.3e");
    ImGui::InputDouble("V_G (J)##Sim12", &V_G, 0, 0, "%.3e");
    ImGui::InputDouble("m* (kg)##Sim12", &mStar, 0, 0, "%.3e");
    ImGui::SliderInt("Bands##Sim12", &numBands, 1, 6);
    ImGui::SliderInt("Points##Sim12", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim12")) {
        clearPlot();
        auto bands = SolidStateSystem::computeNFEBands(a, V_G, mStar, numBands, numPts);
        for (int b = 0; b < (int)bands.bands.size(); ++b) {
            std::vector<double> Ev(bands.bands[b].size());
            for (size_t i = 0; i < bands.bands[b].size(); ++i)
                Ev[i] = bands.bands[b][i] / SSConst::eV;
            addCurve("Band " + std::to_string(b), bands.k_points, Ev);
        }
        double gap = SolidStateSystem::NFEBandGap(V_G);
        std::ostringstream oss;
        oss << "Band gap = 2|V_G| = " << gap/SSConst::eV << " eV";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim12")) {
        SolidStateSystem sys;
        sys.exportNFEBandsCSV("sim12_nfe.csv", a, V_G, mStar, numBands, numPts);
    }
}

void GuiApp::renderSim13_TightBinding() {
    ImGui::Text("13. Tight-Binding (Extended)");
    if (ImGui::CollapsingHeader("Theory: Tight-Binding"))
        ImGui::TextWrapped("Starting from atomic orbitals, the tight-binding model gives E(k) = E0 - 2t*cos(ka) in 1D. The bandwidth is 4t and the effective mass is hbar^2/(2ta^2). Extends to 2D/3D with nearest-neighbor hopping.");

    static double E0 = 0.0;
    static double t = 1.0 * SSConst::eV;
    static double a = 3.0e-10;
    static int numPts = 200;
    ImGui::InputDouble("E0 (J)##Sim13", &E0, 0, 0, "%.3e");
    ImGui::InputDouble("t (J)##Sim13", &t, 0, 0, "%.3e");
    ImGui::InputDouble("a (m)##Sim13", &a, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim13", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim13")) {
        clearPlot();
        auto bands = SolidStateSystem::computeTBBands1D(E0, t, a, numPts);
        for (int b = 0; b < (int)bands.bands.size(); ++b) {
            std::vector<double> Ev(bands.bands[b].size());
            for (size_t i = 0; i < bands.bands[b].size(); ++i)
                Ev[i] = bands.bands[b][i] / SSConst::eV;
            addCurve("Band " + std::to_string(b), bands.k_points, Ev);
        }
        double bw = SolidStateSystem::tightBindingBandwidth(t);
        std::ostringstream oss;
        oss << "Bandwidth = 4t = " << bw/SSConst::eV << " eV";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim13")) {
        SolidStateSystem sys;
        sys.exportTBBandsCSV("sim13_tight_binding.csv", E0, t, a, numPts);
    }
}

void GuiApp::renderSim14_KronigPenney() {
    ImGui::Text("14. Kronig-Penney (Extended)");
    if (ImGui::CollapsingHeader("Theory: Kronig-Penney"))
        ImGui::TextWrapped("Periodic square well potential: exact solution via transfer matrix gives cos(kd) = f(E). Allowed bands where |f(E)| <= 1, gaps where |f(E)| > 1. Delta function limit: cos(ka) + P*sin(ka)/(ka) = cos(Ka).");

    static double V0 = 5.0 * SSConst::eV;
    static double a = 5.0e-10;
    static double b = 1.0e-10;
    static double mStar = SSConst::me;
    static int numBands = 4;
    static int numPts = 200;
    ImGui::InputDouble("V0 (J)##Sim14", &V0, 0, 0, "%.3e");
    ImGui::InputDouble("a (m)##Sim14", &a, 0, 0, "%.3e");
    ImGui::InputDouble("b (m)##Sim14", &b, 0, 0, "%.3e");
    ImGui::InputDouble("m* (kg)##Sim14", &mStar, 0, 0, "%.3e");
    ImGui::SliderInt("Bands##Sim14", &numBands, 1, 8);
    ImGui::SliderInt("Points##Sim14", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim14")) {
        clearPlot();
        auto bands = SolidStateSystem::computeKPBands(V0, a, b, mStar, numBands, numPts);
        for (int bn = 0; bn < (int)bands.bands.size(); ++bn) {
            std::vector<double> Ev(bands.bands[bn].size());
            for (size_t i = 0; i < bands.bands[bn].size(); ++i)
                Ev[i] = bands.bands[bn][i] / SSConst::eV;
            addCurve("Band " + std::to_string(bn), bands.k_points, Ev);
        }
        std::ostringstream oss;
        oss << "V0 = " << V0/SSConst::eV << " eV\n"
            << "Period = a+b = " << (a+b) << " m\n"
            << "Bands computed: " << bands.bands.size();
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim14")) {
        SolidStateSystem sys;
        sys.exportKPBandsCSV("sim14_kronig_penney.csv", V0, a, b, mStar, numBands, numPts);
    }
}

void GuiApp::renderSim15_APW_OPW() {
    ImGui::Text("15. APW / OPW Methods");
    if (ImGui::CollapsingHeader("Theory: APW/OPW"))
        ImGui::TextWrapped("Augmented Plane Wave (APW): muffin-tin potential with spherical solutions inside spheres matched to plane waves outside. OPW orthogonalizes plane waves to core states. Both give realistic band structures.");

    static double a = 5.0e-10;
    static double rMT = 1.5e-10;
    static double V0 = 10.0 * SSConst::eV;
    static int lMax = 3;
    static int numPts = 200;
    ImGui::InputDouble("a (m)##Sim15", &a, 0, 0, "%.3e");
    ImGui::InputDouble("r_MT (m)##Sim15", &rMT, 0, 0, "%.3e");
    ImGui::InputDouble("V0 (J)##Sim15", &V0, 0, 0, "%.3e");
    ImGui::SliderInt("l_max##Sim15", &lMax, 0, 6);
    ImGui::SliderInt("Points##Sim15", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim15")) {
        clearPlot();
        auto bands = SolidStateSystem::computeAPWBands(a, rMT, V0, lMax, numPts);
        for (int bn = 0; bn < (int)bands.bands.size(); ++bn) {
            std::vector<double> Ev(bands.bands[bn].size());
            for (size_t i = 0; i < bands.bands[bn].size(); ++i)
                Ev[i] = bands.bands[bn][i] / SSConst::eV;
            addCurve("Band " + std::to_string(bn), bands.k_points, Ev);
        }
        std::ostringstream oss;
        oss << "APW bands computed: " << bands.bands.size() << "\n"
            << "r_MT/a = " << rMT/a;
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim15")) {
        SolidStateSystem sys;
        sys.exportAPWBandsCSV("sim15_apw.csv", a, rMT, V0, lMax, numPts);
    }
}

void GuiApp::renderSim16_DOS_FermiSurface() {
    ImGui::Text("16. DOS & Fermi Surface");
    if (ImGui::CollapsingHeader("Theory: DOS & Fermi Surface"))
        ImGui::TextWrapped("The density of states g(E) = sum_n integral dk delta(E - E_n(k)). The Fermi surface is the constant-energy surface at E = E_F in k-space. Its topology determines many electronic properties.");

    static double E0 = 0.0;
    static double t = 1.0 * SSConst::eV;
    static double a = 3.0e-10;
    static double broadening = 0.1 * SSConst::eV;
    static int numBins = 200;
    ImGui::InputDouble("E0 (J)##Sim16", &E0, 0, 0, "%.3e");
    ImGui::InputDouble("t (J)##Sim16", &t, 0, 0, "%.3e");
    ImGui::InputDouble("a (m)##Sim16", &a, 0, 0, "%.3e");
    ImGui::InputDouble("Broadening (J)##Sim16", &broadening, 0, 0, "%.3e");
    ImGui::SliderInt("Bins##Sim16", &numBins, 50, 500);

    if (ImGui::Button("Compute##Sim16")) {
        clearPlot();
        auto bands = SolidStateSystem::computeTBBands1D(E0, t, a, numBins);
        auto dos = SolidStateSystem::computeDOSFromBands(bands, broadening, numBins);
        double Emin = E0 - 2.5*std::abs(t), Emax = E0 + 2.5*std::abs(t);
        std::vector<double> Ev(numBins);
        for (int i = 0; i < numBins; ++i)
            Ev[i] = (Emin + (Emax - Emin) * i / (numBins - 1)) / SSConst::eV;
        std::vector<double> dosV(dos.begin(), dos.begin() + std::min((int)dos.size(), numBins));
        dosV.resize(numBins, 0.0);
        addCurve("DOS(E)", Ev, dosV);
        std::ostringstream oss;
        oss << "TB DOS computed with broadening = " << broadening/SSConst::eV << " eV";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim16")) {
        auto bands = SolidStateSystem::computeTBBands1D(E0, t, a, numBins);
        SolidStateSystem sys;
        sys.exportDOSFromBandsCSV("sim16_dos.csv", bands, broadening, numBins);
    }
}

void GuiApp::renderSim17_EffectiveMass_KP() {
    ImGui::Text("17. Effective Mass & k.p Theory");
    if (ImGui::CollapsingHeader("Theory: Effective Mass & k.p"))
        ImGui::TextWrapped("m* = hbar^2 / (d^2E/dk^2). The k.p method expands the Hamiltonian near band extrema, coupling nearby bands via momentum matrix elements. The Kane model gives non-parabolic bands in narrow-gap semiconductors.");

    static double Eg = 1.42 * SSConst::eV;
    static double P = 1.0e-9;
    static int numPts = 200;
    ImGui::InputDouble("Eg (J)##Sim17", &Eg, 0, 0, "%.3e");
    ImGui::InputDouble("Kane P (J.m)##Sim17", &P, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim17", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim17")) {
        clearPlot();
        std::vector<double> kv, Ecb, Evb;
        double kmax = 5.0e9;
        for (int i = 0; i < numPts; ++i) {
            double k = -kmax + 2.0 * kmax * i / (numPts - 1);
            kv.push_back(k);
            auto [Ec, Ev] = SolidStateSystem::kpTwoBand(k, Eg, P);
            Ecb.push_back(Ec / SSConst::eV);
            Evb.push_back(Ev / SSConst::eV);
        }
        addCurve("CB (k.p)", kv, Ecb);
        addCurve("VB (k.p)", kv, Evb);
        std::ostringstream oss;
        oss << "Eg = " << Eg/SSConst::eV << " eV\n"
            << "Kane parameter P = " << P << " J.m";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim17")) {
        SolidStateSystem sys;
        sys.exportKPTheoryCSV("sim17_kp.csv", Eg, P, numPts);
    }
}

void GuiApp::renderSim18_Pseudopotential() {
    ImGui::Text("18. Pseudopotential Method");
    if (ImGui::CollapsingHeader("Theory: Pseudopotential"))
        ImGui::TextWrapped("Pseudopotentials replace the strong ionic potential and core electrons with a weaker effective potential that reproduces the correct valence electron behavior. Only a few Fourier components V_G are needed.");

    static double a = 5.43e-10;
    static double V3 = -0.21 * SSConst::eV;
    static double V8 = 0.04 * SSConst::eV;
    static double V11 = 0.08 * SSConst::eV;
    static int numBands = 4;
    static int numPts = 200;
    ImGui::InputDouble("a (m)##Sim18", &a, 0, 0, "%.3e");
    ImGui::InputDouble("V3 (J)##Sim18", &V3, 0, 0, "%.3e");
    ImGui::InputDouble("V8 (J)##Sim18", &V8, 0, 0, "%.3e");
    ImGui::InputDouble("V11 (J)##Sim18", &V11, 0, 0, "%.3e");
    ImGui::SliderInt("Bands##Sim18", &numBands, 1, 8);
    ImGui::SliderInt("Points##Sim18", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim18")) {
        clearPlot();
        std::vector<double> VG = {V3, V8, V11};
        auto bands = SolidStateSystem::computePseudopotentialBands(a, VG, numBands, numPts);
        for (int bn = 0; bn < (int)bands.bands.size(); ++bn) {
            std::vector<double> Ev(bands.bands[bn].size());
            for (size_t i = 0; i < bands.bands[bn].size(); ++i)
                Ev[i] = bands.bands[bn][i] / SSConst::eV;
            addCurve("Band " + std::to_string(bn), bands.k_points, Ev);
        }
        std::ostringstream oss;
        oss << "Pseudopotential bands: " << bands.bands.size() << "\n"
            << "V3=" << V3/SSConst::eV << ", V8=" << V8/SSConst::eV << ", V11=" << V11/SSConst::eV << " eV";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim18")) {
        std::vector<double> VG = {V3, V8, V11};
        SolidStateSystem sys;
        sys.exportPseudopotentialCSV("sim18_pseudopotential.csv", a, VG, numBands, numPts);
    }
}

void GuiApp::renderSim19_WannierFunctions() {
    ImGui::Text("19. Wannier Functions");
    if (ImGui::CollapsingHeader("Theory: Wannier Functions"))
        ImGui::TextWrapped("Wannier functions are localized real-space counterparts of Bloch functions: w(r-R) = (1/sqrt(N)) * sum_k e^{-ikR} psi_k(r). They form a complete orthonormal set and are useful for tight-binding parameterization.");

    static double a = 5.0e-10;
    static int nBands = 3;
    static int numPts = 200;
    ImGui::InputDouble("a (m)##Sim19", &a, 0, 0, "%.3e");
    ImGui::SliderInt("nBands##Sim19", &nBands, 1, 8);
    ImGui::SliderInt("Points##Sim19", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim19")) {
        clearPlot();
        std::vector<double> xv, wv;
        std::vector<double> coeffs(nBands, 1.0);
        for (int i = 0; i < numPts; ++i) {
            double x = -3.0 * a + 6.0 * a * i / (numPts - 1);
            xv.push_back(x);
            wv.push_back(SolidStateSystem::wannierFunction1D(x, a, nBands, coeffs));
        }
        addCurve("w(x)", xv, wv);
        double spread = SolidStateSystem::wannierSpread(a, nBands);
        std::ostringstream oss;
        oss << "Wannier spread = " << spread << " m\n"
            << "nBands = " << nBands;
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim19")) {
        SolidStateSystem sys;
        sys.exportWannierCSV("sim19_wannier.csv", a, nBands, numPts);
    }
}

void GuiApp::renderSim20_BandStructureVis() {
    ImGui::Text("20. Band Structure Visualization");
    if (ImGui::CollapsingHeader("Theory: Band Structure"))
        ImGui::TextWrapped("Band structure E_n(k) plotted along high-symmetry paths in the Brillouin zone. Compare free electron, tight-binding, NFE, and Kronig-Penney models.");

    static int modelIdx = 0;
    static double a = 5.0e-10;
    static double param1 = 1.0 * SSConst::eV;
    static double param2 = 0.0;
    static int numBands = 4;
    static int numPts = 200;
    static const char* models[] = {"free_electron","tight_binding","nfe","kronig_penney"};
    ImGui::Combo("Model##Sim20", &modelIdx, models, IM_ARRAYSIZE(models));
    ImGui::InputDouble("a (m)##Sim20", &a, 0, 0, "%.3e");
    ImGui::InputDouble("param1 (J)##Sim20", &param1, 0, 0, "%.3e");
    ImGui::InputDouble("param2 (J)##Sim20", &param2, 0, 0, "%.3e");
    ImGui::SliderInt("Bands##Sim20", &numBands, 1, 8);
    ImGui::SliderInt("Points##Sim20", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim20")) {
        clearPlot();
        auto bands = SolidStateSystem::computeModelBandStructure(models[modelIdx], a, param1, param2, numBands, numPts);
        for (int bn = 0; bn < (int)bands.bands.size(); ++bn) {
            std::vector<double> Ev(bands.bands[bn].size());
            for (size_t i = 0; i < bands.bands[bn].size(); ++i)
                Ev[i] = bands.bands[bn][i] / SSConst::eV;
            addCurve("Band " + std::to_string(bn), bands.k_points, Ev);
        }
        std::ostringstream oss;
        oss << "Model: " << models[modelIdx] << "\n"
            << "Bands: " << bands.bands.size();
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim20")) {
        SolidStateSystem sys;
        sys.exportBandStructureCSV("sim20_bands.csv", models[modelIdx], a, param1, param2, numBands, numPts);
    }
}

void GuiApp::renderSim21_IntrinsicSemiconductors() {
    ImGui::Text("21. Intrinsic Semiconductors");
    if (ImGui::CollapsingHeader("Theory: Intrinsic Semiconductors"))
        ImGui::TextWrapped("In an intrinsic semiconductor n = p = n_i = sqrt(Nc*Nv)*exp(-Eg/2kT). The Fermi level lies near mid-gap. Carrier concentration is exponentially activated with temperature.");

    static double Eg = 1.12 * SSConst::eV;
    static double mStar_e = 0.26 * SSConst::me;
    static double mStar_h = 0.36 * SSConst::me;
    static double Tmin = 100.0, Tmax = 600.0;
    static int numPts = 200;
    ImGui::InputDouble("Eg (J)##Sim21", &Eg, 0, 0, "%.3e");
    ImGui::InputDouble("m*_e (kg)##Sim21", &mStar_e, 0, 0, "%.3e");
    ImGui::InputDouble("m*_h (kg)##Sim21", &mStar_h, 0, 0, "%.3e");
    ImGui::InputDouble("T min (K)##Sim21", &Tmin);
    ImGui::InputDouble("T max (K)##Sim21", &Tmax);
    ImGui::SliderInt("Points##Sim21", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim21")) {
        clearPlot();
        std::vector<double> Tv, niv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            Tv.push_back(T);
            niv.push_back(SolidStateSystem::intrinsicCarrierConcentration(Eg, T, mStar_e, mStar_h));
        }
        addCurve("n_i(T)", Tv, niv);
        double ni300 = SolidStateSystem::intrinsicCarrierConcentration(Eg, 300.0, mStar_e, mStar_h);
        double EF300 = SolidStateSystem::intrinsicFermiLevel(Eg, 300.0, mStar_e, mStar_h);
        std::ostringstream oss;
        oss << "n_i(300K) = " << ni300 << " m^-3\n"
            << "E_F(300K) = " << EF300/SSConst::eV << " eV";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim21")) {
        SolidStateSystem sys;
        sys.exportIntrinsicSemiCSV("sim21_intrinsic.csv", Eg, mStar_e, mStar_h, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim22_ExtrinsicSemiconductors() {
    ImGui::Text("22. Extrinsic Semiconductors");
    if (ImGui::CollapsingHeader("Theory: Extrinsic Semiconductors"))
        ImGui::TextWrapped("Doping with donors (n-type) or acceptors (p-type) shifts the Fermi level and dramatically increases carrier concentration. At low T: freeze-out. At intermediate T: extrinsic regime. At high T: intrinsic regime.");

    static double Eg = 1.12 * SSConst::eV;
    static double Nd = 1.0e22;
    static double Na = 0.0;
    static double mStar_e = 0.26 * SSConst::me;
    static double mStar_h = 0.36 * SSConst::me;
    static double Tmin = 50.0, Tmax = 600.0;
    static int numPts = 200;
    ImGui::InputDouble("Eg (J)##Sim22", &Eg, 0, 0, "%.3e");
    ImGui::InputDouble("Nd (m^-3)##Sim22", &Nd, 0, 0, "%.3e");
    ImGui::InputDouble("Na (m^-3)##Sim22", &Na, 0, 0, "%.3e");
    ImGui::InputDouble("m*_e (kg)##Sim22", &mStar_e, 0, 0, "%.3e");
    ImGui::InputDouble("m*_h (kg)##Sim22", &mStar_h, 0, 0, "%.3e");
    ImGui::InputDouble("T min (K)##Sim22", &Tmin);
    ImGui::InputDouble("T max (K)##Sim22", &Tmax);
    ImGui::SliderInt("Points##Sim22", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim22")) {
        clearPlot();
        std::vector<double> Tv, nv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            Tv.push_back(T);
            nv.push_back(SolidStateSystem::electronConcentration(Nd, Na, Eg, T, mStar_e, mStar_h));
        }
        addCurve("n(T)", Tv, nv);
        double n300 = SolidStateSystem::electronConcentration(Nd, Na, Eg, 300.0, mStar_e, mStar_h);
        double EF300 = SolidStateSystem::extrinsicFermiLevel(Eg, 300.0, Nd, Na, mStar_e, mStar_h);
        std::ostringstream oss;
        oss << "n(300K) = " << n300 << " m^-3\n"
            << "E_F(300K) = " << EF300/SSConst::eV << " eV";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim22")) {
        SolidStateSystem sys;
        sys.exportExtrinsicSemiCSV("sim22_extrinsic.csv", Eg, Nd, Na, mStar_e, mStar_h, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim23_PNJunction() {
    ImGui::Text("23. p-n Junction");
    if (ImGui::CollapsingHeader("Theory: p-n Junction"))
        ImGui::TextWrapped("A p-n junction forms a depletion region with built-in voltage V_bi = (kT/e)*ln(Na*Nd/ni^2). The depletion width W = sqrt(2*eps*(V_bi - V_app)*(Na+Nd)/(e*Na*Nd)). Current follows the Shockley diode equation.");

    static double Na = 1.0e23;
    static double Nd = 1.0e22;
    static double ni = 1.5e16;
    static double epsilon_r = 11.7;
    static double T = 300.0;
    static double V_app = 0.0;
    static int numPts = 200;
    ImGui::InputDouble("Na (m^-3)##Sim23", &Na, 0, 0, "%.3e");
    ImGui::InputDouble("Nd (m^-3)##Sim23", &Nd, 0, 0, "%.3e");
    ImGui::InputDouble("ni (m^-3)##Sim23", &ni, 0, 0, "%.3e");
    ImGui::InputDouble("epsilon_r##Sim23", &epsilon_r);
    ImGui::InputDouble("T (K)##Sim23", &T);
    ImGui::InputDouble("V_app (V)##Sim23", &V_app);
    ImGui::SliderInt("Points##Sim23", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim23")) {
        clearPlot();
        auto pn = SolidStateSystem::computePNJunction(Na, Nd, ni, epsilon_r, T, V_app, numPts);
        addCurve("Potential(x)", pn.x, pn.potential);
        addCurve("E-field(x)", pn.x, pn.efield);
        std::ostringstream oss;
        oss << "V_bi = " << pn.V_bi << " V\n"
            << "W = " << pn.W << " m\n"
            << "xn = " << pn.xn << " m, xp = " << pn.xp << " m";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim23")) {
        SolidStateSystem sys;
        sys.exportPNJunctionCSV("sim23_pn_junction.csv", Na, Nd, ni, epsilon_r, T, V_app, numPts);
    }
}

void GuiApp::renderSim24_SemiOpticalProperties() {
    ImGui::Text("24. Semiconductor Optical Properties");
    if (ImGui::CollapsingHeader("Theory: Semiconductor Optics"))
        ImGui::TextWrapped("Direct transitions: alpha ~ sqrt(E - Eg). Indirect transitions require phonon assistance: alpha ~ (E - Eg +/- E_ph)^2. Luminescence spectrum peaks near the band edge and decays as exp(-(E-Eg)/kT).");

    static double Eg = 1.42 * SSConst::eV;
    static double mStar_r = 0.05 * SSConst::me;
    static double T = 300.0;
    static double Emin = 1.0 * SSConst::eV;
    static double Emax = 3.0 * SSConst::eV;
    static int numPts = 200;
    ImGui::InputDouble("Eg (J)##Sim24", &Eg, 0, 0, "%.3e");
    ImGui::InputDouble("m*_r (kg)##Sim24", &mStar_r, 0, 0, "%.3e");
    ImGui::InputDouble("T (K)##Sim24", &T);
    ImGui::InputDouble("E min (J)##Sim24", &Emin, 0, 0, "%.3e");
    ImGui::InputDouble("E max (J)##Sim24", &Emax, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim24", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim24")) {
        clearPlot();
        std::vector<double> Ev, absv, lumv;
        for (int i = 0; i < numPts; ++i) {
            double E = Emin + (Emax - Emin) * i / (numPts - 1);
            Ev.push_back(E / SSConst::eV);
            absv.push_back(SolidStateSystem::directAbsorptionCoeff(E, Eg, mStar_r));
            lumv.push_back(SolidStateSystem::luminescenceSpectrum(E, Eg, T));
        }
        addCurve("Absorption", Ev, absv);
        addCurve("Luminescence", Ev, lumv);
        std::ostringstream oss;
        oss << "Eg = " << Eg/SSConst::eV << " eV\n"
            << "T = " << T << " K";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim24")) {
        SolidStateSystem sys;
        sys.exportOpticalPropertiesCSV("sim24_optical.csv", Eg, mStar_r, T, Emin, Emax, numPts);
    }
}

void GuiApp::renderSim25_QuantumWells() {
    ImGui::Text("25. Quantum Wells & Heterostructures");
    if (ImGui::CollapsingHeader("Theory: Quantum Wells"))
        ImGui::TextWrapped("Quantum confinement in thin layers: E_n = n^2*pi^2*hbar^2/(2m*L^2). Finite wells have fewer bound states. Superlattices create minibands with width ~4t where t is the inter-well coupling.");

    static double Lz = 10.0e-9;
    static double mStar_w = 0.067 * SSConst::me;
    static int maxN = 5;
    static int numPts = 200;
    ImGui::InputDouble("Lz (m)##Sim25", &Lz, 0, 0, "%.3e");
    ImGui::InputDouble("m*_w (kg)##Sim25", &mStar_w, 0, 0, "%.3e");
    ImGui::SliderInt("Max n##Sim25", &maxN, 1, 10);
    ImGui::SliderInt("Points##Sim25", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim25")) {
        clearPlot();
        SolidStateSystem sys("QW", 5.0e-10, mStar_w);
        for (int n = 1; n <= maxN; ++n) {
            std::vector<double> zv, pv;
            for (int i = 0; i < numPts; ++i) {
                double z = Lz * i / (numPts - 1);
                zv.push_back(z * 1e9);
                pv.push_back(sys.quantumWellWavefunction(n, z, Lz));
            }
            addCurve("psi_" + std::to_string(n), zv, pv);
        }
        std::ostringstream oss;
        oss << "Quantum Well Energies:\n";
        for (int n = 1; n <= maxN; ++n) {
            double En = sys.quantumWellEnergy(n, Lz, mStar_w);
            oss << "E_" << n << " = " << En/SSConst::eV << " eV\n";
        }
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim25")) {
        SolidStateSystem sys("QW", 5.0e-10, mStar_w);
        sys.exportQuantumWellCSV("sim25_quantum_well.csv", Lz, mStar_w, maxN, numPts);
    }
}

void GuiApp::renderSim26_HallEffect() {
    ImGui::Text("26. Hall Effect & Magnetoresistance");
    if (ImGui::CollapsingHeader("Theory: Hall Effect"))
        ImGui::TextWrapped("Hall coefficient R_H = -1/(ne) for electrons. The Hall voltage V_H = IB/(net). Magnetoresistance: Delta_rho/rho ~ (mu*B)^2 for classical transport. Used to determine carrier type, density, and mobility.");

    static double n_carrier = 8.5e28;
    static double mu = 0.01;
    static double Bmin = 0.0, Bmax = 10.0;
    static int numPts = 200;
    ImGui::InputDouble("n (m^-3)##Sim26", &n_carrier, 0, 0, "%.3e");
    ImGui::InputDouble("mu (m^2/Vs)##Sim26", &mu);
    ImGui::InputDouble("B min (T)##Sim26", &Bmin);
    ImGui::InputDouble("B max (T)##Sim26", &Bmax);
    ImGui::SliderInt("Points##Sim26", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim26")) {
        clearPlot();
        std::vector<double> Bv, MRv;
        for (int i = 0; i < numPts; ++i) {
            double B = Bmin + (Bmax - Bmin) * i / (numPts - 1);
            Bv.push_back(B);
            MRv.push_back(SolidStateSystem::magnetoresistanceRatio(mu, B));
        }
        addCurve("MR ratio", Bv, MRv);
        double RH = SolidStateSystem::hallCoefficient(n_carrier, SSConst::e_charge);
        std::ostringstream oss;
        oss << "R_H = " << RH << " m^3/C\n"
            << "mu = " << mu << " m^2/Vs";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim26")) {
        SolidStateSystem sys;
        sys.exportHallEffectCSV("sim26_hall.csv", n_carrier, mu, Bmin, Bmax, numPts);
    }
}

void GuiApp::renderSim27_BoltzmannTransport() {
    ImGui::Text("27. Boltzmann Transport");
    if (ImGui::CollapsingHeader("Theory: Boltzmann Transport"))
        ImGui::TextWrapped("Drude conductivity sigma = ne^2*tau/m. Wiedemann-Franz law: kappa/(sigma*T) = L = pi^2/3*(k_B/e)^2. The Boltzmann equation with relaxation time approximation describes carrier dynamics.");

    static double n = 8.5e28;
    static double mStar = SSConst::me;
    static double Tmin = 10.0, Tmax = 500.0;
    static int numPts = 200;
    ImGui::InputDouble("n (m^-3)##Sim27", &n, 0, 0, "%.3e");
    ImGui::InputDouble("m* (kg)##Sim27", &mStar, 0, 0, "%.3e");
    ImGui::InputDouble("T min (K)##Sim27", &Tmin);
    ImGui::InputDouble("T max (K)##Sim27", &Tmax);
    ImGui::SliderInt("Points##Sim27", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim27")) {
        clearPlot();
        double mu_val = 0.01;
        double tau = SolidStateSystem::relaxationTimeApprox(mu_val, mStar);
        std::vector<double> Tv, WFv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            Tv.push_back(T);
            WFv.push_back(SolidStateSystem::wiedemannFranzRatio(T));
        }
        addCurve("WF ratio", Tv, WFv);
        double sigma = SolidStateSystem::drudeCondutivity(n, tau, mStar);
        double L0 = SolidStateSystem::lorenzNumber();
        std::ostringstream oss;
        oss << "sigma = " << sigma << " S/m\n"
            << "tau = " << tau << " s\n"
            << "Lorenz number = " << L0 << " W*Ohm/K^2";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim27")) {
        SolidStateSystem sys;
        sys.exportBoltzmannTransportCSV("sim27_boltzmann.csv", n, mStar, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim28_ThermoelectricEffects() {
    ImGui::Text("28. Thermoelectric Effects");
    if (ImGui::CollapsingHeader("Theory: Thermoelectrics"))
        ImGui::TextWrapped("Seebeck effect: voltage from temperature gradient, S = -pi^2*k_B^2*T/(3eE_F). Peltier effect: heat transport by current, Pi = ST. Figure of merit ZT = S^2*sigma*T/kappa determines thermoelectric efficiency.");

    static double EF = 5.0 * SSConst::eV;
    static double mStar = SSConst::me;
    static double sigma = 1.0e6;
    static double kappa = 10.0;
    static double Tmin = 100.0, Tmax = 800.0;
    static int numPts = 200;
    ImGui::InputDouble("E_F (J)##Sim28", &EF, 0, 0, "%.3e");
    ImGui::InputDouble("m* (kg)##Sim28", &mStar, 0, 0, "%.3e");
    ImGui::InputDouble("sigma (S/m)##Sim28", &sigma, 0, 0, "%.3e");
    ImGui::InputDouble("kappa (W/mK)##Sim28", &kappa);
    ImGui::InputDouble("T min (K)##Sim28", &Tmin);
    ImGui::InputDouble("T max (K)##Sim28", &Tmax);
    ImGui::SliderInt("Points##Sim28", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim28")) {
        clearPlot();
        std::vector<double> Tv, ZTv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            double S = SolidStateSystem::seebeckCoefficient(T, EF, mStar);
            double ZT = SolidStateSystem::thermoelectricZT(S, sigma, kappa, T);
            Tv.push_back(T);
            ZTv.push_back(ZT);
        }
        addCurve("ZT(T)", Tv, ZTv);
        double S300 = SolidStateSystem::seebeckCoefficient(300.0, EF, mStar);
        std::ostringstream oss;
        oss << "S(300K) = " << S300*1e6 << " uV/K\n"
            << "Pi(300K) = " << SolidStateSystem::peltierCoefficient(S300, 300.0) << " V";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim28")) {
        SolidStateSystem sys;
        sys.exportThermoelectricCSV("sim28_thermoelectric.csv", EF, mStar, sigma, kappa, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim29_DiaParamagnetism() {
    ImGui::Text("29. Diamagnetism & Paramagnetism");
    if (ImGui::CollapsingHeader("Theory: Dia- & Paramagnetism"))
        ImGui::TextWrapped("Diamagnetism: all materials, chi < 0 (Langevin, Landau). Paramagnetism: unpaired spins, chi > 0, Curie law chi = C/T. Brillouin function B_J(x) describes quantum paramagnetism. Pauli paramagnetism for conduction electrons.");

    static double J = 3.5;
    static double g = 2.0;
    static int N = 1000;
    static double Bmax = 10.0;
    static double T = 300.0;
    static int numPts = 200;
    ImGui::InputDouble("J (quantum)##Sim29", &J);
    ImGui::InputDouble("g-factor##Sim29", &g);
    ImGui::InputInt("N spins##Sim29", &N);
    ImGui::InputDouble("B max (T)##Sim29", &Bmax);
    ImGui::InputDouble("T (K)##Sim29", &T);
    ImGui::SliderInt("Points##Sim29", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim29")) {
        clearPlot();
        std::vector<double> Bv, Mv;
        for (int i = 0; i < numPts; ++i) {
            double B = Bmax * i / (numPts - 1);
            Bv.push_back(B);
            Mv.push_back(SolidStateSystem::curieParaMagnetization(B, T, J, g, N));
        }
        addCurve("M(B)", Bv, Mv);
        std::ostringstream oss;
        oss << "J = " << J << ", g = " << g << "\n"
            << "T = " << T << " K, N = " << N;
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim29")) {
        SolidStateSystem sys;
        sys.exportParamagnetismCSV("sim29_paramag.csv", J, g, N, Bmax, T, numPts);
    }
}

void GuiApp::renderSim30_Ferromagnetism() {
    ImGui::Text("30. Ferromagnetism");
    if (ImGui::CollapsingHeader("Theory: Ferromagnetism"))
        ImGui::TextWrapped("Spontaneous magnetization below Curie temperature T_C. Mean-field (Weiss) theory: self-consistent equation M = M_sat*B_J(g*mu_B*J*lambda*M/(kT)). Spin waves (magnons) with dispersion omega ~ k^2.");

    static double n = 8.5e28;
    static double J = 0.5;
    static double g = 2.0;
    static double lambda_mf = 1000.0;
    static double Tmin = 1.0, Tmax = 1200.0;
    static int numPts = 200;
    ImGui::InputDouble("n (m^-3)##Sim30", &n, 0, 0, "%.3e");
    ImGui::InputDouble("J (quantum)##Sim30", &J);
    ImGui::InputDouble("g-factor##Sim30", &g);
    ImGui::InputDouble("lambda_mf##Sim30", &lambda_mf);
    ImGui::InputDouble("T min (K)##Sim30", &Tmin);
    ImGui::InputDouble("T max (K)##Sim30", &Tmax);
    ImGui::SliderInt("Points##Sim30", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim30")) {
        clearPlot();
        auto res = SolidStateSystem::solveMeanFieldFerromagnet(n, J, g, lambda_mf, Tmin, Tmax, numPts);
        addCurve("M(T)", res.temperature, res.magnetization);
        addCurve("chi(T)", res.temperature, res.susceptibility);
        double Tc = SolidStateSystem::curieWeissTemperature(n, J, g, lambda_mf);
        std::ostringstream oss;
        oss << "T_C (Weiss) = " << Tc << " K";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim30")) {
        SolidStateSystem sys;
        sys.exportFerromagnetismCSV("sim30_ferromagnet.csv", n, J, g, lambda_mf, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim31_HeisenbergModel() {
    ImGui::Text("31. Heisenberg Model");
    if (ImGui::CollapsingHeader("Theory: Heisenberg Model"))
        ImGui::TextWrapped("H = -2J * sum S_i . S_j. Isotropic exchange coupling. Spin wave excitations E(k) = 4JS(1 - cos(ka)). Bloch T^(3/2) law for magnetization reduction due to magnons.");

    static double J = 1.0e-21;
    static double S = 0.5;
    static double a = 3.0e-10;
    static int numPts = 200;
    ImGui::InputDouble("J (J)##Sim31", &J, 0, 0, "%.3e");
    ImGui::InputDouble("S##Sim31", &S);
    ImGui::InputDouble("a (m)##Sim31", &a, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim31", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim31")) {
        clearPlot();
        std::vector<double> kv, Ev;
        double kmax = SSConst::PI / a;
        for (int i = 0; i < numPts; ++i) {
            double k = -kmax + 2.0 * kmax * i / (numPts - 1);
            kv.push_back(k);
            Ev.push_back(SolidStateSystem::spinWaveEnergy(k, a, J, S));
        }
        addCurve("E_magnon(k)", kv, Ev);
        double E0 = SolidStateSystem::heisenbergGroundStateEnergy1D(100, J, S);
        std::ostringstream oss;
        oss << "Ground state E (N=100) = " << E0 << " J\n"
            << "Bandwidth = " << SolidStateSystem::spinWaveEnergy(SSConst::PI/a, a, J, S) << " J";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim31")) {
        SolidStateSystem sys;
        sys.exportHeisenbergCSV("sim31_heisenberg.csv", J, S, a, numPts);
    }
}

void GuiApp::renderSim32_IsingModel() {
    ImGui::Text("32. Ising Model (1D/2D)");
    if (ImGui::CollapsingHeader("Theory: Ising Model"))
        ImGui::TextWrapped("H = -J * sum s_i*s_j with s = +/-1. 1D: exact solution, no phase transition. 2D: Onsager exact solution, T_C = 2J/(k_B*ln(1+sqrt(2))). Monte Carlo simulation with Metropolis algorithm.");

    static double J = 1.0e-21;
    static int L = 16;
    static double Tmin = 0.5, Tmax = 5.0;
    static int numTemps = 30;
    static int mcSteps = 5000;
    ImGui::InputDouble("J (J)##Sim32", &J, 0, 0, "%.3e");
    ImGui::SliderInt("L (lattice)##Sim32", &L, 4, 32);
    ImGui::InputDouble("T min (K)##Sim32", &Tmin);
    ImGui::InputDouble("T max (K)##Sim32", &Tmax);
    ImGui::SliderInt("Temps##Sim32", &numTemps, 5, 50);
    ImGui::SliderInt("MC steps##Sim32", &mcSteps, 1000, 20000);

    if (ImGui::Button("Compute##Sim32")) {
        clearPlot();
        auto res = SolidStateSystem::ising2DMonteCarloSimulation(L, J, Tmin, Tmax, numTemps, mcSteps);
        addCurve("|M|(T)", res.temperature, res.magnetization);
        addCurve("Cv(T)", res.temperature, res.specificHeat);
        double Tc = SolidStateSystem::ising2DCriticalTemperature(J);
        std::ostringstream oss;
        oss << "T_c (Onsager) = " << Tc << " K\n"
            << "L = " << L << ", MC steps = " << mcSteps;
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim32")) {
        SolidStateSystem sys;
        sys.exportIsing2DCSV("sim32_ising.csv", L, J, Tmin, Tmax, numTemps, mcSteps);
    }
}

void GuiApp::renderSim33_Antiferromagnetism() {
    ImGui::Text("33. Antiferromagnetism");
    if (ImGui::CollapsingHeader("Theory: Antiferromagnetism"))
        ImGui::TextWrapped("Neighboring spins align antiparallel below the Neel temperature T_N. Above T_N: chi = C/(T + T_N). Below T_N: chi_parallel decreases linearly, chi_perpendicular is constant.");

    static double J = 1.0e-21;
    static double z = 6.0;
    static double S = 0.5;
    static double Tmin = 1.0, Tmax = 500.0;
    static int numPts = 200;
    ImGui::InputDouble("J (J)##Sim33", &J, 0, 0, "%.3e");
    ImGui::InputDouble("z (neighbors)##Sim33", &z);
    ImGui::InputDouble("S##Sim33", &S);
    ImGui::InputDouble("T min (K)##Sim33", &Tmin);
    ImGui::InputDouble("T max (K)##Sim33", &Tmax);
    ImGui::SliderInt("Points##Sim33", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim33")) {
        clearPlot();
        double TN = SolidStateSystem::neelTemperature(J, z, S);
        double C = 1.0;
        double chi_TN = SolidStateSystem::afmSusceptibilityAboveTN(TN + 0.01, TN, C);
        std::vector<double> Tv, chiv;
        for (int i = 0; i < numPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (numPts - 1);
            Tv.push_back(T);
            if (T > TN)
                chiv.push_back(SolidStateSystem::afmSusceptibilityAboveTN(T, TN, C));
            else
                chiv.push_back(SolidStateSystem::afmSusceptibilityBelowTN_parallel(T, TN, chi_TN));
        }
        addCurve("chi(T)", Tv, chiv);
        std::ostringstream oss;
        oss << "T_N = " << TN << " K";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim33")) {
        SolidStateSystem sys;
        sys.exportAntiferromagnetismCSV("sim33_afm.csv", J, z, S, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim34_Ferrimagnetism() {
    ImGui::Text("34. Ferrimagnetism");
    if (ImGui::CollapsingHeader("Theory: Ferrimagnetism"))
        ImGui::TextWrapped("Two magnetic sublattices with unequal moments. Net magnetization M = |M_A - M_B|. Compensation temperature where M = 0. Ferrites are technologically important ferrimagnets.");

    static double JA = 1.0e-21, JB = 0.5e-21, JAB = -0.8e-21;
    static double SA = 2.5, SB = 2.0;
    static int NA = 1000, NB = 1000;
    static double Tmin = 1.0, Tmax = 800.0;
    static int numPts = 200;
    ImGui::InputDouble("J_A (J)##Sim34", &JA, 0, 0, "%.3e");
    ImGui::InputDouble("J_B (J)##Sim34", &JB, 0, 0, "%.3e");
    ImGui::InputDouble("J_AB (J)##Sim34", &JAB, 0, 0, "%.3e");
    ImGui::InputDouble("S_A##Sim34", &SA);
    ImGui::InputDouble("S_B##Sim34", &SB);
    ImGui::InputInt("N_A##Sim34", &NA);
    ImGui::InputInt("N_B##Sim34", &NB);
    ImGui::InputDouble("T min (K)##Sim34", &Tmin);
    ImGui::InputDouble("T max (K)##Sim34", &Tmax);
    ImGui::SliderInt("Points##Sim34", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim34")) {
        clearPlot();
        auto res = SolidStateSystem::ferriMeanField(JA, JB, JAB, SA, SB, NA, NB, Tmin, Tmax, numPts);
        addCurve("M(T)", res.temperature, res.magnetization);
        double Tcomp = SolidStateSystem::ferriCompensationTemperature(JA, JB, SA, SB);
        std::ostringstream oss;
        oss << "Compensation T = " << Tcomp << " K";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim34")) {
        SolidStateSystem sys;
        sys.exportFerrimagnetismCSV("sim34_ferri.csv", JA, JB, JAB, SA, SB, NA, NB, Tmin, Tmax, numPts);
    }
}

void GuiApp::renderSim35_MagneticDomains() {
    ImGui::Text("35. Magnetic Domains");
    if (ImGui::CollapsingHeader("Theory: Magnetic Domains"))
        ImGui::TextWrapped("Domains minimize total energy (exchange + anisotropy + magnetostatic). Domain wall width delta = pi*sqrt(A/K). Domain wall energy sigma = 4*sqrt(AK). Hysteresis from domain wall pinning and nucleation.");

    static double A_ex = 1.0e-11;
    static double K_aniso = 5.0e4;
    static double Ms = 1.0e6;
    static double Hmax = 1.0e5;
    static int numPts = 200;
    ImGui::InputDouble("A_ex (J/m)##Sim35", &A_ex, 0, 0, "%.3e");
    ImGui::InputDouble("K_aniso (J/m^3)##Sim35", &K_aniso, 0, 0, "%.3e");
    ImGui::InputDouble("Ms (A/m)##Sim35", &Ms, 0, 0, "%.3e");
    ImGui::InputDouble("H_max (A/m)##Sim35", &Hmax, 0, 0, "%.3e");
    ImGui::SliderInt("Points##Sim35", &numPts, 50, 500);

    if (ImGui::Button("Compute##Sim35")) {
        clearPlot();
        std::vector<double> Hv, Mv;
        double Hc = K_aniso / Ms;
        double Mr = 0.8 * Ms;
        for (int i = 0; i < numPts; ++i) {
            double H = -Hmax + 2.0 * Hmax * i / (numPts - 1);
            Hv.push_back(H);
            Mv.push_back(SolidStateSystem::hysteresisLoop(H, Hc, Ms, Mr));
        }
        addCurve("M(H)", Hv, Mv);
        double dw = SolidStateSystem::domainWallWidth(A_ex, K_aniso);
        double dwe = SolidStateSystem::domainWallEnergy(A_ex, K_aniso);
        double rc = SolidStateSystem::singleDomainCriticalRadius(A_ex, K_aniso, Ms);
        std::ostringstream oss;
        oss << "Wall width = " << dw << " m\n"
            << "Wall energy = " << dwe << " J/m^2\n"
            << "Critical radius = " << rc << " m";
        m_resultText = oss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim35")) {
        SolidStateSystem sys;
        sys.exportMagneticDomainsCSV("sim35_domains.csv", A_ex, K_aniso, Ms, Hmax, numPts);
    }
}

void GuiApp::renderSim36_SpinOrbitCoupling() {
    ImGui::Text("36. Spin-Orbit Coupling in Solids");
    if (ImGui::CollapsingHeader("Theory: Spin-Orbit Coupling"))
        ImGui::TextWrapped("Rashba effect: H_R = alpha_R*(sigma x k).z from structural inversion asymmetry. Dresselhaus: from bulk inversion asymmetry. Spin-orbit splits bands and enables spin Hall effect and topological states.");

    static double alpha_R = 1e-11;
    static double beta_D = 0.5e-11;
    static double mStar = SSConst::mStar_GaAs * SSConst::me;
    static int nPts = 200;

    ImGui::InputDouble("Rashba alpha (J*m)##Sim36", &alpha_R, 0, 0, "%.3e");
    ImGui::InputDouble("Dresselhaus beta (J*m)##Sim36", &beta_D, 0, 0, "%.3e");
    ImGui::InputDouble("Effective mass (kg)##Sim36", &mStar, 0, 0, "%.4e");
    ImGui::InputInt("Points##Sim36", &nPts);

    if (ImGui::Button("Compute##Sim36")) {
        clearPlot();
        double kMax = 5e9;
        std::vector<double> kv(nPts), Eup(nPts), Edn(nPts), Edress(nPts);
        for (int i = 0; i < nPts; ++i) {
            double k = -kMax + 2.0 * kMax * i / (nPts - 1);
            kv[i] = k * 1e-9;
            auto [Ep, Em] = SolidStateSystem::rashbaSplit(k, alpha_R, mStar);
            Eup[i] = Ep / SSConst::eV;
            Edn[i] = Em / SSConst::eV;
            Edress[i] = SolidStateSystem::dresselhausEnergy(k, beta_D, mStar) / SSConst::eV;
        }
        addCurve("Rashba +", kv, Eup);
        addCurve("Rashba -", kv, Edn);
        addCurve("Dresselhaus", kv, Edress);
        std::ostringstream ss;
        ss << "Rashba alpha = " << alpha_R << " J*m\n"
           << "Dresselhaus beta = " << beta_D << " J*m";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim36")) {
        SolidStateSystem sys;
        sys.exportSpinOrbitCSV("sim36_spinorbit.csv", alpha_R, beta_D, mStar, nPts);
    }
}

void GuiApp::renderSim37_BCSTheory() {
    ImGui::Text("37. BCS Theory");
    if (ImGui::CollapsingHeader("Theory: BCS"))
        ImGui::TextWrapped("Cooper pairs form via electron-phonon interaction below T_C. Gap equation gives Delta(T). T_C = 1.13*omega_D*exp(-1/(N(E_F)*V)). Specific heat jump Delta_C/(gamma*T_C) = 1.43. Quasiparticle DOS has gap.");

    static double omegaD = 3e13;
    static double N_EF = 1e47;
    static double V_eff = 0.3e-47;
    static int nPts = 200;

    ImGui::InputDouble("Debye freq omega_D (rad/s)##Sim37", &omegaD, 0, 0, "%.3e");
    ImGui::InputDouble("DOS at E_F N(E_F) (1/J)##Sim37", &N_EF, 0, 0, "%.3e");
    ImGui::InputDouble("Coupling V_eff (J)##Sim37", &V_eff, 0, 0, "%.3e");
    ImGui::InputInt("Points##Sim37", &nPts);

    if (ImGui::Button("Compute##Sim37")) {
        clearPlot();
        BCSResult res = SolidStateSystem::computeBCSGapVsTemp(omegaD, N_EF, V_eff, nPts);
        std::vector<double> Tnorm(res.temperature.size()), gapEv(res.gap.size());
        for (size_t i = 0; i < res.temperature.size(); ++i) {
            Tnorm[i] = (res.Tc > 0) ? res.temperature[i] / res.Tc : res.temperature[i];
            gapEv[i] = res.gap[i] / SSConst::eV;
        }
        addCurve("Delta(T/Tc)", Tnorm, gapEv);
        std::ostringstream ss;
        ss << "Tc = " << res.Tc << " K\n"
           << "Delta(0) = " << res.Delta0 / SSConst::eV << " eV\n"
           << "Specific heat jump ratio = " << SolidStateSystem::bcsSpecificHeatJump();
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim37")) {
        SolidStateSystem sys;
        sys.exportBCSTheoryCSV("sim37_bcs.csv", omegaD, N_EF, V_eff, nPts);
    }
}

void GuiApp::renderSim38_GinzburgLandau() {
    ImGui::Text("38. Ginzburg-Landau Theory");
    if (ImGui::CollapsingHeader("Theory: Ginzburg-Landau"))
        ImGui::TextWrapped("Order parameter psi(r). Coherence length xi = xi_0/sqrt(1-T/Tc). Penetration depth lambda = lambda_0/sqrt(1-T/Tc). GL parameter kappa = lambda/xi. Type I: kappa < 1/sqrt(2), Type II: kappa > 1/sqrt(2).");

    static double xi0 = 80e-9;
    static double lambda0 = 40e-9;
    static double Tc = 9.2;
    static double Tmin = 0.1;
    static double Tmax = 9.0;
    static int nPts = 200;

    ImGui::InputDouble("xi_0 (m)##Sim38", &xi0, 0, 0, "%.3e");
    ImGui::InputDouble("lambda_0 (m)##Sim38", &lambda0, 0, 0, "%.3e");
    ImGui::InputDouble("Tc (K)##Sim38", &Tc);
    ImGui::InputDouble("T min (K)##Sim38", &Tmin);
    ImGui::InputDouble("T max (K)##Sim38", &Tmax);
    ImGui::InputInt("Points##Sim38", &nPts);

    if (ImGui::Button("Compute##Sim38")) {
        clearPlot();
        std::vector<double> Tv(nPts), xiv(nPts), lamv(nPts);
        for (int i = 0; i < nPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (nPts - 1);
            Tv[i] = T;
            xiv[i] = SolidStateSystem::glCoherenceLength(xi0, T, Tc) * 1e9;
            lamv[i] = SolidStateSystem::glPenetrationDepth(lambda0, T, Tc) * 1e9;
        }
        addCurve("xi (nm)", Tv, xiv);
        addCurve("lambda (nm)", Tv, lamv);
        std::ostringstream ss;
        double kappa0 = SolidStateSystem::glKappaParameter(lambda0, xi0);
        ss << "kappa(0) = " << kappa0 << "\n"
           << (kappa0 < 1.0 / std::sqrt(2.0) ? "Type I" : "Type II") << " superconductor";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim38")) {
        SolidStateSystem sys;
        sys.exportGLTheoryCSV("sim38_gl.csv", xi0, lambda0, Tc, Tmin, Tmax, nPts);
    }
}

void GuiApp::renderSim39_TypeI_II_SC() {
    ImGui::Text("39. Type I & II Superconductors");
    if (ImGui::CollapsingHeader("Theory: Type I & II"))
        ImGui::TextWrapped("Type I: complete Meissner effect, single H_c. Type II: mixed state between H_c1 and H_c2 with Abrikosov vortex lattice. Each vortex carries one flux quantum Phi_0 = h/(2e).");

    static double Hc0 = 1e5;
    static double kappa = 0.5;
    static double Tc = 7.2;
    static int nPts = 200;

    ImGui::InputDouble("Hc0 (A/m)##Sim39", &Hc0, 0, 0, "%.3e");
    ImGui::InputDouble("kappa##Sim39", &kappa);
    ImGui::InputDouble("Tc (K)##Sim39", &Tc);
    ImGui::InputInt("Points##Sim39", &nPts);

    if (ImGui::Button("Compute##Sim39")) {
        clearPlot();
        std::vector<double> Tv(nPts), Hcv(nPts), Hc1v(nPts), Hc2v(nPts);
        for (int i = 0; i < nPts; ++i) {
            double T = 0.01 + (Tc - 0.02) * i / (nPts - 1);
            Tv[i] = T;
            double Hc = SolidStateSystem::thermodynamicCriticalField(Hc0, T, Tc);
            Hcv[i] = Hc;
            Hc1v[i] = SolidStateSystem::glLowerCriticalField(Hc, kappa);
            Hc2v[i] = SolidStateSystem::glUpperCriticalField(Hc, kappa);
        }
        addCurve("Hc (A/m)", Tv, Hcv);
        addCurve("Hc1 (A/m)", Tv, Hc1v);
        addCurve("Hc2 (A/m)", Tv, Hc2v);
        std::ostringstream ss;
        ss << "Flux quantum = " << SolidStateSystem::fluxQuantum() << " Wb\n"
           << (SolidStateSystem::isTypeII(kappa) ? "Type II" : "Type I") << " (kappa = " << kappa << ")";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim39")) {
        SolidStateSystem sys;
        sys.exportSCTypesCSV("sim39_sctypes.csv", Hc0, kappa, Tc, nPts);
    }
}

void GuiApp::renderSim40_JosephsonEffect() {
    ImGui::Text("40. Josephson Effect");
    if (ImGui::CollapsingHeader("Theory: Josephson Effect"))
        ImGui::TextWrapped("DC Josephson: I = I_c*sin(phi). AC Josephson: dphi/dt = 2eV/hbar, frequency f = 2eV/h. SQUID: quantum interference of two junctions. Shapiro steps at V_n = n*h*f/(2e).");

    static double Ic = 1e-6;
    static double Rn = 10.0;
    static double freq = 10e9;
    static int nPts = 300;

    ImGui::InputDouble("Ic (A)##Sim40", &Ic, 0, 0, "%.3e");
    ImGui::InputDouble("Rn (Ohm)##Sim40", &Rn);
    ImGui::InputDouble("Microwave freq (Hz)##Sim40", &freq, 0, 0, "%.3e");
    ImGui::InputInt("Points##Sim40", &nPts);

    if (ImGui::Button("Compute##Sim40")) {
        clearPlot();
        std::vector<double> phi(nPts), Idc(nPts);
        for (int i = 0; i < nPts; ++i) {
            double p = -SSConst::PI + SSConst::TWO_PI * i / (nPts - 1);
            phi[i] = p;
            Idc[i] = SolidStateSystem::josephsonDCCurrent(Ic, p);
        }
        addCurve("I(phi) (A)", phi, Idc);
        std::ostringstream ss;
        double fAC = SolidStateSystem::josephsonACFrequency(1e-4);
        ss << "AC freq at 0.1 mV = " << fAC / 1e9 << " GHz\n"
           << "Shapiro V_1 = " << SolidStateSystem::shapiroStepVoltage(1, freq) * 1e6 << " uV\n"
           << "Shapiro V_2 = " << SolidStateSystem::shapiroStepVoltage(2, freq) * 1e6 << " uV";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim40")) {
        SolidStateSystem sys;
        sys.exportJosephsonCSV("sim40_josephson.csv", Ic, Rn, Ic * Rn, nPts);
    }
}

void GuiApp::renderSim41_ElectronPhonon() {
    ImGui::Text("41. Electron-Phonon Interaction");
    if (ImGui::CollapsingHeader("Theory: Electron-Phonon"))
        ImGui::TextWrapped("Electron-phonon coupling constant lambda = 2*N(E_F)*<g^2>/omega_ph. McMillan formula for T_c. Bloch-Gruneisen resistivity: rho ~ T^5 at low T, rho ~ T at high T.");

    static double thetaD = 343.0;
    static double rho0 = 1e-8;
    static double K = 5e-7;
    static double Tmin = 1.0;
    static double Tmax = 500.0;
    static int nPts = 200;

    ImGui::InputDouble("Debye temp (K)##Sim41", &thetaD);
    ImGui::InputDouble("Residual rho_0 (Ohm*m)##Sim41", &rho0, 0, 0, "%.3e");
    ImGui::InputDouble("BG constant K (Ohm*m)##Sim41", &K, 0, 0, "%.3e");
    ImGui::InputDouble("T min (K)##Sim41", &Tmin);
    ImGui::InputDouble("T max (K)##Sim41", &Tmax);
    ImGui::InputInt("Points##Sim41", &nPts);

    if (ImGui::Button("Compute##Sim41")) {
        clearPlot();
        std::vector<double> Tv(nPts), rhov(nPts);
        for (int i = 0; i < nPts; ++i) {
            double T = Tmin + (Tmax - Tmin) * i / (nPts - 1);
            Tv[i] = T;
            rhov[i] = SolidStateSystem::resistivityBlochGruneisen(T, thetaD, rho0, K) * 1e8;
        }
        addCurve("rho (uOhm*cm)", Tv, rhov);
        std::ostringstream ss;
        ss << "Theta_D = " << thetaD << " K\n"
           << "Residual rho_0 = " << rho0 * 1e8 << " uOhm*cm";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim41")) {
        SolidStateSystem sys;
        sys.exportElectronPhononCSV("sim41_ephonon.csv", thetaD, rho0, K, Tmin, Tmax, nPts);
    }
}

void GuiApp::renderSim42_AndersonLocalization() {
    ImGui::Text("42. Anderson Localization");
    if (ImGui::CollapsingHeader("Theory: Anderson Localization"))
        ImGui::TextWrapped("Disorder-induced localization. 1D/2D: all states localized for any disorder. 3D: mobility edge separates extended/localized states. Ioffe-Regel criterion k_F*l ~ 1. Thouless energy E_Th = hbar*D/L^2.");

    static double mfp = 50e-9;
    static double kF = 1e10;
    static double Lmin = 1e-9;
    static double Lmax = 1e-6;
    static int nPts = 200;

    ImGui::InputDouble("Mean free path (m)##Sim42", &mfp, 0, 0, "%.3e");
    ImGui::InputDouble("Fermi wavevector kF (1/m)##Sim42", &kF, 0, 0, "%.3e");
    ImGui::InputDouble("L min (m)##Sim42", &Lmin, 0, 0, "%.3e");
    ImGui::InputDouble("L max (m)##Sim42", &Lmax, 0, 0, "%.3e");
    ImGui::InputInt("Points##Sim42", &nPts);

    if (ImGui::Button("Compute##Sim42")) {
        clearPlot();
        double xi1D = SolidStateSystem::andersonLocalizationLength1D(mfp);
        double g0 = kF * mfp;
        std::vector<double> Lv(nPts), gv(nPts);
        for (int i = 0; i < nPts; ++i) {
            double L = Lmin + (Lmax - Lmin) * i / (nPts - 1);
            Lv[i] = L * 1e9;
            gv[i] = SolidStateSystem::conductanceDimensionless(g0, L, xi1D);
        }
        addCurve("g(L)", Lv, gv);
        std::ostringstream ss;
        ss << "xi_1D = " << xi1D * 1e9 << " nm\n"
           << "xi_2D = " << SolidStateSystem::andersonLocalizationLength2D(mfp) * 1e9 << " nm\n"
           << "Ioffe-Regel kF*l = " << SolidStateSystem::ioffRegelCriterion(kF, mfp);
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim42")) {
        SolidStateSystem sys;
        sys.exportAndersonCSV("sim42_anderson.csv", mfp, kF, Lmin, Lmax, nPts);
    }
}

void GuiApp::renderSim43_QuantumHallEffect() {
    ImGui::Text("43. Quantum Hall Effect (Extended)");
    if (ImGui::CollapsingHeader("Theory: QHE"))
        ImGui::TextWrapped("Landau levels E_n = hbar*omega_c*(n+1/2). Filling factor nu = n_2D*h/(eB). Integer QHE: sigma_xy = nu*e^2/h exactly. Fractional QHE: Laughlin states at nu = 1/3, 1/5, etc.");

    static double n2D = 3e15;
    static double mStar = SSConst::mStar_GaAs * SSConst::me;
    static double Bmin = 0.5;
    static double Bmax = 15.0;
    static double broadening = 0.5;
    static int nPts = 500;

    ImGui::InputDouble("n_2D (1/m^2)##Sim43", &n2D, 0, 0, "%.3e");
    ImGui::InputDouble("Effective mass (kg)##Sim43", &mStar, 0, 0, "%.4e");
    ImGui::InputDouble("B min (T)##Sim43", &Bmin);
    ImGui::InputDouble("B max (T)##Sim43", &Bmax);
    ImGui::InputDouble("Broadening##Sim43", &broadening);
    ImGui::InputInt("Points##Sim43", &nPts);

    if (ImGui::Button("Compute##Sim43")) {
        clearPlot();
        std::vector<double> Bv(nPts), Rxy(nPts), Rxx(nPts);
        for (int i = 0; i < nPts; ++i) {
            double B = Bmin + (Bmax - Bmin) * i / (nPts - 1);
            Bv[i] = B;
            double nu = SolidStateSystem::fillingFactor(n2D, B);
            Rxy[i] = SolidStateSystem::hallResistance(nu);
            Rxx[i] = SolidStateSystem::longitudinalResistance(B, n2D, mStar, broadening);
        }
        addCurve("R_xy (Ohm)", Bv, Rxy);
        addCurve("R_xx (Ohm)", Bv, Rxx);
        std::ostringstream ss;
        ss << "Magnetic length at 5T = " << SolidStateSystem::magneticLength(5.0) * 1e9 << " nm\n";
        for (int n = 0; n < 4; ++n)
            ss << "E_" << n << " at 5T = " << SolidStateSystem::landauLevel(n, 5.0, mStar) / SSConst::eV * 1e3 << " meV\n";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim43")) {
        SolidStateSystem sys;
        sys.exportQuantumHallCSV("sim43_qhe.csv", n2D, mStar, Bmin, Bmax, nPts);
    }
}

void GuiApp::renderSim44_TopologicalInsulators() {
    ImGui::Text("44. Topological Insulators (Toy)");
    if (ImGui::CollapsingHeader("Theory: Topological Insulators"))
        ImGui::TextWrapped("SSH model: H = sum (v*c_A^dag*c_B + w*c_B^dag*c_A'). Topological phase when |w| > |v| (winding number = 1). Edge states appear in the gap. Chern insulator in 2D.");

    static double v = 1.0;
    static double w = 1.5;
    static int nPts = 300;

    ImGui::InputDouble("Intra-cell hopping v##Sim44", &v);
    ImGui::InputDouble("Inter-cell hopping w##Sim44", &w);
    ImGui::InputInt("Points##Sim44", &nPts);

    if (ImGui::Button("Compute##Sim44")) {
        clearPlot();
        BandResult res = SolidStateSystem::computeSSHBands(v, w, nPts);
        for (size_t b = 0; b < res.bands.size(); ++b)
            addCurve(b == 0 ? "Lower" : "Upper", res.k_points, res.bands[b]);
        std::ostringstream ss;
        double wn = SolidStateSystem::windingNumber(v, w);
        ss << "Winding number = " << wn << "\n"
           << (wn > 0.5 ? "Topological phase (edge states)" : "Trivial phase");
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim44")) {
        SolidStateSystem sys;
        sys.exportTopologicalCSV("sim44_topo.csv", v, w, nPts);
    }
}

void GuiApp::renderSim45_GrapheneDirac() {
    ImGui::Text("45. Graphene & Dirac Materials");
    if (ImGui::CollapsingHeader("Theory: Graphene"))
        ImGui::TextWrapped("Honeycomb lattice with E(k) = +/- t*|f(k)|. Linear dispersion near K points: E = hbar*v_F*|k|. Dirac fermions with v_F ~ 10^6 m/s. Landau levels: E_n = v_F*sqrt(2*e*hbar*B*|n|). Anomalous QHE.");

    static double t = 2.8;
    static double a = 2.46e-10;
    static double B_ll = 10.0;
    static int nPts = 300;

    ImGui::InputDouble("Hopping t (eV)##Sim45", &t);
    ImGui::InputDouble("Lattice constant a (m)##Sim45", &a, 0, 0, "%.3e");
    ImGui::InputDouble("B for Landau levels (T)##Sim45", &B_ll);
    ImGui::InputInt("Points##Sim45", &nPts);

    if (ImGui::Button("Compute##Sim45")) {
        clearPlot();
        double K_pt = 4.0 * SSConst::PI / (3.0 * a);
        std::vector<double> kv(nPts), Eup(nPts), Edn(nPts);
        for (int i = 0; i < nPts; ++i) {
            double kx = -K_pt + 2.0 * K_pt * i / (nPts - 1);
            kv[i] = kx * a;
            auto [Ep, Em] = SolidStateSystem::grapheneBands(kx, 0.0, t, a);
            Eup[i] = Ep;
            Edn[i] = Em;
        }
        addCurve("Upper band (eV)", kv, Eup);
        addCurve("Lower band (eV)", kv, Edn);
        std::ostringstream ss;
        double vF = SolidStateSystem::grapheneFermiVelocity(t * SSConst::eV, a);
        ss << "v_F = " << vF / 1e6 << " x 10^6 m/s\n";
        for (int n = 0; n < 4; ++n)
            ss << "LL(" << n << ") = " << SolidStateSystem::grapheneLandauLevel(n, B_ll) / SSConst::eV * 1e3 << " meV\n";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim45")) {
        SolidStateSystem sys;
        sys.exportGrapheneCSV("sim45_graphene.csv", t, a, nPts);
    }
}

void GuiApp::renderSim46_PlasmonsDielectric() {
    ImGui::Text("46. Plasmons & Dielectric Function");
    if (ImGui::CollapsingHeader("Theory: Plasmons"))
        ImGui::TextWrapped("Plasma frequency omega_p = sqrt(ne^2/(m*eps_0)). Drude: eps(omega) = 1 - omega_p^2/(omega^2 + i*gamma*omega). Lindhard: quantum dielectric function. Thomas-Fermi screening. Plasmon dispersion.");

    static double n_el = 8.5e28;
    static double mStar = SSConst::me;
    static double gamma_d = 1e13;
    static double omegaMax = 3e16;
    static int nPts = 500;

    ImGui::InputDouble("Electron density (1/m^3)##Sim46", &n_el, 0, 0, "%.3e");
    ImGui::InputDouble("Effective mass (kg)##Sim46", &mStar, 0, 0, "%.4e");
    ImGui::InputDouble("Damping gamma (rad/s)##Sim46", &gamma_d, 0, 0, "%.3e");
    ImGui::InputDouble("omega max (rad/s)##Sim46", &omegaMax, 0, 0, "%.3e");
    ImGui::InputInt("Points##Sim46", &nPts);

    if (ImGui::Button("Compute##Sim46")) {
        clearPlot();
        double wp = SolidStateSystem::plasmaFrequency(n_el, mStar, 1.0);
        std::vector<double> wv(nPts), epsR(nPts), epsI(nPts);
        for (int i = 0; i < nPts; ++i) {
            double w = 1e12 + (omegaMax - 1e12) * i / (nPts - 1);
            wv[i] = w / wp;
            auto eps = SolidStateSystem::drudePermittivity(w, wp, gamma_d);
            epsR[i] = eps.real();
            epsI[i] = eps.imag();
        }
        addCurve("Re(eps)", wv, epsR);
        addCurve("Im(eps)", wv, epsI);
        std::ostringstream ss;
        ss << "omega_p = " << wp / 1e15 << " x 10^15 rad/s\n"
           << "Thomas-Fermi l_TF = " << SolidStateSystem::thomasFermiScreeningLength(n_el, mStar) * 1e10 << " A";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim46")) {
        SolidStateSystem sys;
        sys.exportPlasmonCSV("sim46_plasmon.csv", n_el, mStar, gamma_d, omegaMax, nPts);
    }
}

void GuiApp::renderSim47_Excitons() {
    ImGui::Text("47. Excitons");
    if (ImGui::CollapsingHeader("Theory: Excitons"))
        ImGui::TextWrapped("Bound electron-hole pairs: E_n = Eg - Ry*/(n^2) where Ry* = (m_r/m_e)*(1/eps_r^2)*13.6 eV. Wannier excitons in semiconductors (large radius). Mott density marks exciton ionization.");

    static double mStar_e = 0.067;
    static double mStar_h = 0.45;
    static double epsilon_r = 12.9;
    static double Eg = 1.42;
    static double broadening = 0.005;
    static int nPts = 500;

    ImGui::InputDouble("m*_e / m_e##Sim47", &mStar_e);
    ImGui::InputDouble("m*_h / m_e##Sim47", &mStar_h);
    ImGui::InputDouble("epsilon_r##Sim47", &epsilon_r);
    ImGui::InputDouble("Eg (eV)##Sim47", &Eg);
    ImGui::InputDouble("Broadening (eV)##Sim47", &broadening);
    ImGui::InputInt("Points##Sim47", &nPts);

    if (ImGui::Button("Compute##Sim47")) {
        clearPlot();
        double mr = SolidStateSystem::excitonReducedMass(mStar_e * SSConst::me, mStar_h * SSConst::me);
        double Eb = SolidStateSystem::excitonBindingEnergy(1, mr, epsilon_r);
        std::vector<double> Ev(nPts), absv(nPts);
        double Emin = Eg - 0.1, Emax = Eg + 0.5;
        for (int i = 0; i < nPts; ++i) {
            double E = Emin + (Emax - Emin) * i / (nPts - 1);
            Ev[i] = E;
            absv[i] = SolidStateSystem::excitonAbsorption(E * SSConst::eV, Eg * SSConst::eV,
                                                           Eb, broadening * SSConst::eV);
        }
        addCurve("Absorption", Ev, absv);
        std::ostringstream ss;
        double aB = SolidStateSystem::excitonBohrRadius(mr, epsilon_r);
        ss << "Exciton binding E_1 = " << Eb / SSConst::eV * 1e3 << " meV\n"
           << "Bohr radius = " << aB * 1e9 << " nm\n"
           << "Mott density (3D) = " << SolidStateSystem::mottDensity(aB, 3) << " m^-3";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim47")) {
        SolidStateSystem sys;
        sys.exportExcitonCSV("sim47_exciton.csv", mStar_e * SSConst::me, mStar_h * SSConst::me,
                              epsilon_r, Eg * SSConst::eV, nPts);
    }
}

void GuiApp::renderSim48_Polaritons() {
    ImGui::Text("48. Polaritons");
    if (ImGui::CollapsingHeader("Theory: Polaritons"))
        ImGui::TextWrapped("Coupled light-matter modes. Upper and lower polariton branches with anticrossing. Rabi splitting = 2g. Hopfield coefficients describe photon/exciton content. Exciton-polariton condensation.");

    static double omega_c = 2.5e15;
    static double omega_exc = 2.5e15;
    static double g_coupling = 5e13;
    static int nPts = 300;

    ImGui::InputDouble("Cavity omega_c (rad/s)##Sim48", &omega_c, 0, 0, "%.3e");
    ImGui::InputDouble("Exciton omega_exc (rad/s)##Sim48", &omega_exc, 0, 0, "%.3e");
    ImGui::InputDouble("Coupling g (rad/s)##Sim48", &g_coupling, 0, 0, "%.3e");
    ImGui::InputInt("Points##Sim48", &nPts);

    if (ImGui::Button("Compute##Sim48")) {
        clearPlot();
        double kMax = 1e7;
        std::vector<double> kv(nPts), Eup(nPts), Elow(nPts);
        for (int i = 0; i < nPts; ++i) {
            double k = -kMax + 2.0 * kMax * i / (nPts - 1);
            kv[i] = k * 1e-6;
            auto [wUp, wLow] = SolidStateSystem::polaritonDispersion(k, omega_c, omega_exc, g_coupling);
            Eup[i] = SSConst::hbar * wUp / SSConst::eV * 1e3;
            Elow[i] = SSConst::hbar * wLow / SSConst::eV * 1e3;
        }
        addCurve("Upper polariton (meV)", kv, Eup);
        addCurve("Lower polariton (meV)", kv, Elow);
        std::ostringstream ss;
        ss << "Rabi splitting = " << SolidStateSystem::rabiSplitting(g_coupling) / 1e12 << " THz\n"
           << "              = " << SSConst::hbar * SolidStateSystem::rabiSplitting(g_coupling) / SSConst::eV * 1e3 << " meV";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim48")) {
        SolidStateSystem sys;
        sys.exportPolaritonCSV("sim48_polariton.csv", omega_c, omega_exc, g_coupling, nPts);
    }
}

void GuiApp::renderSim49_AmorphousSolids() {
    ImGui::Text("49. Amorphous Solids");
    if (ImGui::CollapsingHeader("Theory: Amorphous Solids"))
        ImGui::TextWrapped("No long-range order. Radial distribution function g(r) shows short-range order. Two-level systems give linear specific heat C ~ T and thermal conductivity plateau. Anderson localization of phonons.");

    static double r0 = 2.5e-10;
    static double sigma = 0.3e-10;
    static double n0 = 5e28;
    static double alpha = 1e-3;
    static double beta = 1e-5;
    static double Tmin = 0.1;
    static double Tmax = 100.0;
    static int nPts = 200;

    ImGui::InputDouble("r0 (m)##Sim49", &r0, 0, 0, "%.3e");
    ImGui::InputDouble("sigma (m)##Sim49", &sigma, 0, 0, "%.3e");
    ImGui::InputDouble("n0 (1/m^3)##Sim49", &n0, 0, 0, "%.3e");
    ImGui::InputDouble("TLS alpha##Sim49", &alpha, 0, 0, "%.3e");
    ImGui::InputDouble("Phonon beta##Sim49", &beta, 0, 0, "%.3e");
    ImGui::InputDouble("T min (K)##Sim49", &Tmin);
    ImGui::InputDouble("T max (K)##Sim49", &Tmax);
    ImGui::InputInt("Points##Sim49", &nPts);

    if (ImGui::Button("Compute##Sim49")) {
        clearPlot();
        std::vector<double> rv(nPts), gv(nPts), Tv(nPts), Cv(nPts);
        double rMin = 1e-10, rMax = 10e-10;
        for (int i = 0; i < nPts; ++i) {
            double r = rMin + (rMax - rMin) * i / (nPts - 1);
            rv[i] = r * 1e10;
            gv[i] = SolidStateSystem::radialDistributionFunction(r, r0, sigma, n0);
            double T = Tmin + (Tmax - Tmin) * i / (nPts - 1);
            Tv[i] = T;
            Cv[i] = SolidStateSystem::amorphousSpecificHeat(T, alpha, beta);
        }
        addCurve("g(r)", rv, gv);
        addCurve("C(T)", Tv, Cv);
        std::ostringstream ss;
        ss << "RDF peak at r0 = " << r0 * 1e10 << " A\n"
           << "Low-T specific heat: C ~ " << alpha << "*T + " << beta << "*T^3";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim49")) {
        SolidStateSystem sys;
        sys.exportAmorphousSolidsCSV("sim49_amorphous.csv", r0, sigma, n0, Tmin, Tmax, nPts);
    }
}

void GuiApp::renderSim50_MesoscopicTransport() {
    ImGui::Text("50. Mesoscopic Transport");
    if (ImGui::CollapsingHeader("Theory: Mesoscopic Transport"))
        ImGui::TextWrapped("Landauer formula: G = (2e^2/h)*T. Conductance quantization in QPCs. Coulomb blockade: E_c = e^2/(2C). Aharonov-Bohm oscillations in rings. Universal conductance fluctuations.");

    static double barrierH = 0.3 * SSConst::eV;
    static double barrierW = 2e-9;
    static double mStar = SSConst::me;
    static double Emin = 0.01 * SSConst::eV;
    static double Emax = 0.6 * SSConst::eV;
    static double C_dot = 1e-18;
    static int nPts = 300;

    ImGui::InputDouble("Barrier height (J)##Sim50", &barrierH, 0, 0, "%.3e");
    ImGui::InputDouble("Barrier width (m)##Sim50", &barrierW, 0, 0, "%.3e");
    ImGui::InputDouble("Effective mass (kg)##Sim50", &mStar, 0, 0, "%.4e");
    ImGui::InputDouble("E min (J)##Sim50", &Emin, 0, 0, "%.3e");
    ImGui::InputDouble("E max (J)##Sim50", &Emax, 0, 0, "%.3e");
    ImGui::InputDouble("QD capacitance C (F)##Sim50", &C_dot, 0, 0, "%.3e");
    ImGui::InputInt("Points##Sim50", &nPts);

    if (ImGui::Button("Compute##Sim50")) {
        clearPlot();
        TransportResult res = SolidStateSystem::computeLandauerTransport(
            barrierH, barrierW, mStar, Emin, Emax, nPts);
        std::vector<double> EeV(res.energy.size());
        for (size_t i = 0; i < res.energy.size(); ++i)
            EeV[i] = res.energy[i] / SSConst::eV;
        addCurve("Transmission", EeV, res.transmission);
        std::vector<double> Gnorm(res.conductance.size());
        double G0 = SolidStateSystem::conductanceQuantum();
        for (size_t i = 0; i < res.conductance.size(); ++i)
            Gnorm[i] = res.conductance[i] / G0;
        addCurve("G / G_0", EeV, Gnorm);
        std::ostringstream ss;
        ss << "G_0 = " << G0 * 1e6 << " uS\n"
           << "Coulomb E_c = " << SolidStateSystem::coulombBlockadeChargingEnergy(C_dot) / SSConst::eV * 1e3 << " meV\n"
           << "AB period (100 nm ring) = " << SolidStateSystem::aharonovBohmOscillationPeriod(SSConst::PI * 1e-14) << " T";
        m_resultText = ss.str();
    }
    ImGui::SameLine();
    if (ImGui::Button("Export CSV##Sim50")) {
        SolidStateSystem sys;
        sys.exportMesoscopicCSV("sim50_mesoscopic.csv", barrierH, barrierW, mStar, Emin, Emax, nPts);
    }
}
