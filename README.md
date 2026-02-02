# Microgrid System Identification (BESS Dynamics)

This repository presents an end-to-end **system identification workflow** applied to a real-world **Battery Energy Storage System (BESS)** in a microgrid. The project focuses on extracting reliable dynamic models from measured data and rigorously validating them using industry-relevant criteria.

The emphasis is on **model credibility, validation, and engineering decision-making**, not just curve fitting.

---

## Project Summary

- **System:** Inverter-based Battery Energy Storage System (BESS)
- **Application:** Microgrid power dynamics
- **Input:** Random excitation of real power demand
- **Output:** Measured real power response
- **Objective:** Identify low-order dynamic models suitable for prediction and control

---

## Engineering Approach

### Signal & Excitation Analysis
- Evaluated input autocorrelation and power spectral density
- Identified **low-frequency-dominated, non-white excitation**
- Assessed frequency ranges where models are reliable

### Model Identification
- **Non-parametric frequency response estimation** (SPA)
- **FIR modeling**
  - Order selection based on impulse response decay and noise sensitivity
  - Confidence bounds and residual testing
- **State-space realization**
  - Model order selection using Hankel singular values
  - Cross-validation against FIR dynamics
- **ARX modeling**
  - Structure selection using residual correlation tests
  - Minimal-order model satisfying validation criteria

### Validation Strategy
- Residual–input cross-correlation tests
- Comparison of simulation vs one-step-ahead prediction
- Explicit rejection of models that fail validation, even when visually accurate

---

## Key Outcomes

- FIR model captured long-memory dynamics and **passed validation**
- State-space realization showed good impulse fit but **failed residual tests**
- Low-order ARX model:
  - Passed all validation criteria
  - Delivered strong one-step-ahead prediction performance
  - Provided the best balance of simplicity and accuracy

**Takeaway:** Model validity matters more than model complexity.

---

## Why This Matters (Industry Context)

- Demonstrates **data-driven modeling under non-ideal excitation**
- Shows disciplined **model validation and rejection**
- Directly relevant to:
  - Energy systems & power electronics
  - Robotics and control systems
  - System identification for real hardware
- Reflects engineering judgment expected in production-level modeling

---

## Repository Structure

