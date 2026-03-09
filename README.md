# BL4S Electromagnetic Shower Simulation

> **Beams for Schools 2026** — Investigating Electromagnetic Shower Development and Inverse Material Characterisation Using GEANT4

A GEANT4-based Monte Carlo simulation of high-energy electron showers through various absorber materials, developed as part of a CERN Beams for Schools experiment proposal. The simulation models a 4×4 lead-glass calorimeter array and implements both forward shower prediction and inverse material characterisation via observable extraction.

---

## Physics Overview

When a high-energy electron enters matter it undergoes **bremsstrahlung**, emitting a photon. That photon undergoes **pair production**, creating a new electron–positron pair. The cycle repeats, producing a cascade of lower-energy particles — an **electromagnetic shower**.

The shower's shape encodes information about the material it passed through, specifically its **radiation length** X₀ — the characteristic distance over which an electron loses ~63% of its energy. This simulation investigates:

- **Forward problem** — given a known material, predict shower behaviour
- **Inverse problem** — given shower observables, estimate the material's radiation length

---

## Detector Layout

```
Beam (e⁻, 1–4 GeV)
     │
     ▼
┌─────────────┐
│ Absorber    │  ← material under test (Al / Fe / Cu / Pb / composite)
│ slab        │    thickness varied per run
└─────────────┘
      ·
      · (gap)
      ·
┌────┬────┬────┬────┐
│    │    │    │    │
├────┼────┼────┼────┤
│    │ ██ │ ██ │    │  ← 4×4 lead-glass calorimeter blocks
├────┼────┼────┼────┤     each 5×5×30 cm
│    │ ██ │ ██ │    │     central 2×2 = "core"
├────┼────┼────┼────┤
│    │    │    │    │
└────┴────┴────┴────┘
```

Each calorimeter block records the total energy deposited per event. From the 16-block readout, three shower observables are extracted:

| Observable | Description |
|---|---|
| `Etotal` | Total energy deposited across all 16 blocks |
| `CoreFraction` | Energy in central 2×2 blocks / total — measures compactness |
| `ShowerWidth` | Energy-weighted 2D radial RMS — measures lateral spread (proxy for Molière radius) |

---

## Requirements

- **GEANT4 11.x** with data packages
- **CMake 3.16+**
- **Python 3.8+** with `pandas`, `numpy`, `matplotlib`, `seaborn` for analysis



---

## Building

```bash
# Source the GEANT4 environment first (required every session)
source /usr/bin/geant4.sh

mkdir build
cd build
cmake ..
make -j4
```

> **Nushell users:** GEANT4's setup script is bash-only. Drop into bash for the build and run steps:
> ```
> bash
> source /usr/bin/geant4.sh
> cd build && make -j4
> ```

---

## File Structure

```
BL4S_Sim/
├── BL4S_sim.cc                   ← main entry point
├── CMakeLists.txt
│
├── include/
│   ├── ActionInitialisation.hh   ← registers all user actions
│   ├── CalorimeterSD.hh          ← sensitive detector + hit class
│   ├── DetectorConstruction.hh   ← geometry
│   └── DetectorMessenger.hh      ← macro command interface
│
├── src/
│   ├── ActionInitialisation.cc   ← primary generator, run action, event action
│   ├── CalorimeterSD.cc          ← energy scoring per block
│   ├── DetectorConstruction.cc   ← world + absorber + calorimeter volumes
│   └── DetectorMessenger.cc      ← /det/ macro commands
│
├── macros/
│   ├── run_iron.mac              ← iron thickness scan (5 thicknesses × 1000 events)
│   ├── run_material_scan.mac     ← Al / Fe / Cu / Pb comparison
│   └── vis.mac                   ← interactive visualisation (5 events)
│
└── analyse_shower.py             ← CSV analysis + plots
```

---

## Running

```bash
cd build

# Iron absorber, 5 thickness steps
./BL4S_sim macros/run_iron.mac

# Material comparison: Al, Fe, Cu, Pb
./BL4S_sim macros/run_material_scan.mac

# Interactive visualisation
./BL4S_sim
```

### Changing parameters via macro

```
# Change material only (no /det/update needed — LV material updated directly)
/det/setAbsorberMaterial G4_Pb

# Change thickness (triggers automatic geometry rebuild internally)
/det/setAbsorberThickness 20 mm

# /det/update explicitly forces a full geometry rebuild — only needed if you
# want to manually trigger a rebuild for another reason.
/gun/energy 4 GeV
/run/beamOn 2000
```

### Available materials (NIST database)

| Macro name | Material | Z | X₀ (cm) |
|---|---|---|---|
| `G4_Al` | Aluminium | 13 | 8.90 |
| `G4_Fe` | Iron | 26 | 1.76 |
| `G4_Cu` | Copper | 29 | 1.44 |
| `G4_Pb` | Lead | 82 | 0.56 |

---

## Output & Analysis

Each run produces `shower_output.csv` in the build directory. Every row is one simulated event:

```
event_id, material, absorber_thickness_mm, absorber_x0,
beam_energy_GeV, beam_px_MeVc, beam_py_MeVc, beam_pz_MeVc,
vertex_x_mm, vertex_y_mm,
E_block_0_GeV, ..., E_block_15_GeV,
Etotal_GeV, CoreFraction, ShowerWidth_cm
```

### Analyse results

```bash
pip install pandas numpy matplotlib seaborn
python analyse_shower.py shower_output.csv
```

Produces 7 plots in the `plots/` directory:
- `1_energy_resolution.png` — energy resolution σ/E vs beam energy per material
- `2_shower_attenuation.png` — mean Etotal vs absorber thickness (X₀ units)
- `3_core_fraction.png` — core fraction vs absorber thickness
- `4_shower_width.png` — lateral shower width (2D radial RMS) vs absorber thickness
- `5_block_energy_map.png` — mean energy heatmap for each of the 16 blocks
- `6_lateral_profile.png` — normalised energy fraction per column, per material
- `7_correlation_matrix.png` — Pearson correlation matrix of all observables

---

## Noise Filtering

Three levels of noise rejection are applied automatically before writing to CSV:

| Cut | Location | What it removes |
|---|---|---|
| Sub-threshold deposits | `CalorimeterSD.cc` | Step deposits below 0.5 MeV (electronic noise) |
| Beam-miss events | `EventAction` | Events where `Etotal < 5%` of beam energy |
| Statistical outliers | `EventAction` | Events more than 3.5σ from running mean (Welford algorithm) |

Thresholds are defined as constants and easy to tune:
```cpp
// CalorimeterSD.cc
static constexpr G4double kBlockThreshold = 0.5 * MeV;

// ActionInitialisation.cc (NoiseFilter struct)
static constexpr double kMinFraction  = 0.05;   // beam-miss cut
static constexpr double kOutlierSigma = 3.5;    // outlier cut
```

---

## Connecting to the Experiment

### Forward problem
Expected trends to verify against theory and real data:

| Trend | Observable |
|---|---|
| Higher Z → shorter X₀ → shower starts earlier | `Etotal` saturates at smaller thickness for Pb vs Al |
| Higher Z → more compact shower | `CoreFraction` higher for Pb |
| Higher Z → narrower lateral profile | `ShowerWidth` smaller for Pb |

### Inverse problem (ML)
The CSV output is structured as a training dataset for regression:

```
features: [Etotal, CoreFraction, ShowerWidth, E_block_0..15]
target:   radiation length X₀  (known from GEANT4 material properties)
```

Train a model (e.g. Random Forest, neural network) on simulated data, then apply it to a blind material to estimate its X₀ within a target uncertainty of ±15–20%.

### Heitler model check
The depth of shower maximum scales as:
```
t_max = ln(E₀ / Ec) / ln(2)   [in units of X₀]
```
Scan beam energy from 1–4 GeV with a fixed absorber and verify this logarithmic dependence through the change in `CoreFraction` vs thickness.

---

## Acknowledgements

Developed as part of a [CERN Beams for Schools](https://beams.cern.ch/) experiment proposal on electromagnetic shower development and inverse material characterisation.

Physics references:
- [PDG Review — Passage of Particles Through Matter](https://pdg.lbl.gov/2023/reviews/rpp2023-rev-passage-particles-matter.pdf)
- W. Heitler, *The Quantum Theory of Radiation* (1954)
- [GEANT4 Physics Reference Manual](https://geant4-userdoc.web.cern.ch/UsersGuides/PhysicsReferenceManual/html/index.html)
