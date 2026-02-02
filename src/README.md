# Source Code & Data

This directory contains the MATLAB source code and dataset used for the system identification study of a microgrid Battery Energy Storage System (BESS).

---

## Contents

- `final_project_code.m`  
  End-to-end MATLAB script that performs:
  - Data loading and preprocessing
  - Non-parametric identification (autocorrelation, input spectrum, FRF estimation)
  - FIR model estimation and validation
  - State-space realization and order selection
  - ARX model identification and validation
  - Simulation vs prediction performance comparison
  - Figure generation used in the final report

  This script reproduces **all analysis and figures** presented in the final project report.

- `rand_Pinv_Pbess.mat`  
  Input–output dataset used for system identification.

---

## Data Acknowledgment

The dataset (`rand_Pinv_Pbess.mat`) was **provided by Professor Raymond de Callafon and the System Identification and Control Laboratory** for academic use as part of MAE 283A.

The data represents measured input–output signals from a real microgrid system and is included here **solely for educational and demonstration purposes**.

---

## Notes on Usage

- The code is written to run as a single script for clarity and reproducibility.
- No additional toolboxes beyond standard MATLAB system identification and signal processing functionality are required.
- Paths are assumed to be relative to the `src/` directory.

---

## Reproducibility

Running `final_project_code.m` will:
1. Load the provided dataset
2. Execute the full system identification pipeline
3. Generate all figures referenced in the final report

---

## Author

**Shivharsh Kand**  
Graduate Student — Controls & Mechatronics
