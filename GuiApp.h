#pragma once
// ============================================================================
//  GuiApp.h  —  ImGui + ImPlot GUI for SolidStateCore (50 simulations)
// ============================================================================

#include <string>
#include <vector>

struct GLFWwindow;

class GuiApp {
public:
    GuiApp();
    ~GuiApp();

    bool init(int width = 1600, int height = 900, const char* title = "SolidStateCore");
    void run();
    void shutdown();

private:
    GLFWwindow* m_window = nullptr;
    int m_selectedSim = 0;
    std::string m_resultText;

    // Plot helpers
    std::vector<std::vector<double>> m_plotX;
    std::vector<std::vector<double>> m_plotY;
    std::vector<std::string>         m_plotLabels;
    void clearPlot();
    void addCurve(const std::string& label,
                  const std::vector<double>& x,
                  const std::vector<double>& y);

    // ── 50 render functions ────────────────────────────────────────────────
    void renderSim01_BravaisLattices();
    void renderSim02_ReciprocalLattice();
    void renderSim03_XRayDiffraction();
    void renderSim04_MonatomicLattice();
    void renderSim05_DiatomicLattice();
    void renderSim06_DebyeModel();
    void renderSim07_EinsteinModel();
    void renderSim08_PhononDOS();
    void renderSim09_ThermalExpansion();
    void renderSim10_LatticeThermalConductivity();
    void renderSim11_FreeElectronModel();
    void renderSim12_NearlyFreeElectron();
    void renderSim13_TightBinding();
    void renderSim14_KronigPenney();
    void renderSim15_APW_OPW();
    void renderSim16_DOS_FermiSurface();
    void renderSim17_EffectiveMass_KP();
    void renderSim18_Pseudopotential();
    void renderSim19_WannierFunctions();
    void renderSim20_BandStructureVis();
    void renderSim21_IntrinsicSemiconductors();
    void renderSim22_ExtrinsicSemiconductors();
    void renderSim23_PNJunction();
    void renderSim24_SemiOpticalProperties();
    void renderSim25_QuantumWells();
    void renderSim26_HallEffect();
    void renderSim27_BoltzmannTransport();
    void renderSim28_ThermoelectricEffects();
    void renderSim29_DiaParamagnetism();
    void renderSim30_Ferromagnetism();
    void renderSim31_HeisenbergModel();
    void renderSim32_IsingModel();
    void renderSim33_Antiferromagnetism();
    void renderSim34_Ferrimagnetism();
    void renderSim35_MagneticDomains();
    void renderSim36_SpinOrbitCoupling();
    void renderSim37_BCSTheory();
    void renderSim38_GinzburgLandau();
    void renderSim39_TypeI_II_SC();
    void renderSim40_JosephsonEffect();
    void renderSim41_ElectronPhonon();
    void renderSim42_AndersonLocalization();
    void renderSim43_QuantumHallEffect();
    void renderSim44_TopologicalInsulators();
    void renderSim45_GrapheneDirac();
    void renderSim46_PlasmonsDielectric();
    void renderSim47_Excitons();
    void renderSim48_Polaritons();
    void renderSim49_AmorphousSolids();
    void renderSim50_MesoscopicTransport();
};
