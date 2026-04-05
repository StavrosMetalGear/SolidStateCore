# SolidStateCore

An interactive condensed-matter physics simulator with 50 simulations covering crystal structure, electronic bands, magnetism, superconductivity, and mesoscopic transport. Built with C++17, ImGui, and ImPlot.

![C++17](https://img.shields.io/badge/C%2B%2B-17-blue)
![CMake](https://img.shields.io/badge/CMake-%E2%89%A53.15-green)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

---

## Features

| #   | Simulation                        | #   | Simulation                        |
|-----|-----------------------------------|-----|-----------------------------------|
| 1   | Bravais Lattices                  | 26  | Hall Effect & Magnetoresistance   |
| 2   | Reciprocal Lattice                | 27  | Boltzmann Transport               |
| 3   | X-Ray Diffraction                 | 28  | Thermoelectric Effects            |
| 4   | Monatomic Lattice Vibrations      | 29  | Dia- & Paramagnetism              |
| 5   | Diatomic Lattice Vibrations       | 30  | Ferromagnetism                    |
| 6   | Debye Model                       | 31  | Heisenberg Model                  |
| 7   | Einstein Model                    | 32  | Ising Model (1D/2D)              |
| 8   | Phonon DOS                        | 33  | Antiferromagnetism                |
| 9   | Thermal Expansion                 | 34  | Ferrimagnetism                    |
| 10  | Lattice Thermal Conductivity      | 35  | Magnetic Domains                  |
| 11  | Free Electron Model               | 36  | Spin-Orbit Coupling               |
| 12  | Nearly Free Electron              | 37  | BCS Theory                        |
| 13  | Tight-Binding (Extended)          | 38  | Ginzburg-Landau Theory            |
| 14  | Kronig-Penney (Extended)          | 39  | Type I & II Superconductors       |
| 15  | APW / OPW Methods                 | 40  | Josephson Effect                  |
| 16  | DOS & Fermi Surface               | 41  | Electron-Phonon Interaction       |
| 17  | Effective Mass & k·p              | 42  | Anderson Localization             |
| 18  | Pseudopotential Method            | 43  | Quantum Hall Effect (Extended)    |
| 19  | Wannier Functions                 | 44  | Topological Insulators (Toy)      |
| 20  | Band Structure Visualization      | 45  | Graphene & Dirac Materials        |
| 21  | Intrinsic Semiconductors          | 46  | Plasmons & Dielectric Function    |
| 22  | Extrinsic Semiconductors          | 47  | Excitons                          |
| 23  | p-n Junction                      | 48  | Polaritons                        |
| 24  | Semiconductor Optics              | 49  | Amorphous Solids                  |
| 25  | Quantum Wells & Heterostructures  | 50  | Mesoscopic Transport              |

Every simulation provides:
- **Theory panel** — collapsible header with a concise physics summary.
- **Adjustable parameters** — sliders and input fields for real-time exploration.
- **Live plots** — rendered with [ImPlot](https://github.com/epezent/implot).
- **CSV export** — one-click export for post-processing in Python, MATLAB, etc.

## Architecture

```
SolidStateCore/
├── CMakeLists.txt              # Build configuration (FetchContent for deps)
├── main.cpp                    # Entry point
├── GuiApp.h / GuiApp.cpp       # ImGui/ImPlot GUI (50 render functions)
├── SolidStateSystem.h / .cpp   # Core physics engine (static + instance API)
├── SolidStateConstants.h       # Physical constants (SI units)
├── NumericalSolverSS.h / .cpp  # FDM, Crank-Nicolson, Poisson-Schrödinger, TMM
├── solidstate_c_api.h / .cpp   # Plain-C wrapper for FFI / scripting bindings
├── pch.h                       # Precompiled header (MSVC)
└── Eigen/                      # Bundled Eigen headers (linear algebra)
```

The build produces two targets:

| Target                | Type              | Description                         |
|-----------------------|-------------------|-------------------------------------|
| `SolidStatePhysics`  | Static library    | Physics engine — no GUI dependency  |
| `SolidStateCore`     | Executable        | Full GUI application                |

## Requirements

| Dependency | Version      | Acquired via        |
|------------|--------------|---------------------|
| C++ compiler | C++17      | —                   |
| CMake      | ≥ 3.15       | —                   |
| OpenGL     | ≥ 3.3        | System              |
| GLFW       | 3.4          | FetchContent        |
| Dear ImGui | 1.91.8       | FetchContent        |
| ImPlot     | 0.16         | FetchContent        |
| Eigen      | 3.x          | Bundled in `Eigen/` |

## Building

```bash
# Clone
git clone https://github.com/<your-username>/SolidStateCore.git
cd SolidStateCore

# Configure & build (Ninja recommended)
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build

# Run
./build/SolidStateCore
```

On Windows with Visual Studio:

```powershell
cmake -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release
.\build\Release\SolidStateCore.exe
```

> **Note:** GLFW, ImGui, and ImPlot are downloaded automatically via CMake `FetchContent` — no manual dependency installation is required.

## Usage

1. Launch the executable. A window opens with a **sidebar** listing all 50 simulations.
2. Click a simulation name to load its panel in the center area.
3. Expand the **Theory** header for a physics overview.
4. Adjust parameters and press **Compute** to generate plots and numerical results.
5. Press **Export CSV** to save data to the working directory.

## C API

`solidstate_c_api.h` exposes a subset of the physics engine as plain C functions, making it straightforward to call from Python (ctypes/cffi), C#, or other languages:

```c
double ss_fermi_energy(double n_density, double mStar);
double ss_debye_specific_heat(double T, double thetaD);
double ss_bcs_tc(double omegaD, double N_EF, double V_eff);
// ... and more
```

## Contributing

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/my-simulation`).
3. Commit your changes and open a pull request.

## License

This project is provided under the [MIT License](LICENSE).
