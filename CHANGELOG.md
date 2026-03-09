# Changelog

All notable changes to the **FRF Statistics** MATLAB library will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to Semantic Versioning.

## [Unreleased / Current Update] - 2026-03-09

### Added
- **Input Validation:** Added dimension checking to `FRF_Plotbands.m` to ensure the `time` and `avg` vectors are the exact same length, preventing cryptic MATLAB indexing errors.
- **Didactic Headers:** Completely rewrote the internal MATLAB help blocks (the comments under the `function` declarations) for all core functions. Typing `help FRF_Supervector` (or any other function) now provides detailed, formatted documentation on tensor dimensions, inputs, and outputs.
- **Out-of-Bounds Extraction:** Added a code snippet and methodology in the Reference Manual demonstrating how to dynamically extract and plot boundary violations using logical indexing.
- **Industrial Application Examples:** Expanded the Reference Manual to include a new subsection detailing Process Control and MIMO anomaly detection, bridging the gap between clinical neurology and chemical/industrial engineering.
- **Parameter Tuning Guidance:** Added formal recommendations in the manual for tuning the `alpha` parameter specifically for predictive maintenance and early-warning alarm sensitivity.

### Changed
- **Band Matrix Standardization:** Standardized the `band` output matrix order across `FRF_PredictionBand`, `FRF_ConfidenceBand`, `FRF_ConfidenceBandDifference`, and `FRF_MinimalPredictionBand`. 
  - *Row 1* is now strictly the lower bound (floor).
  - *Row 2* is now strictly the upper bound (ceiling).
- **Vector Orientation:** Forced all time-series outputs (such as `t` and `x` in `FRF_pseudoimpulse` and `FRF_Supervector`) to consistently output as column vectors ($L \times 1$). This prevents implicit expansion crashes during logical comparisons.
- **Manual Clarifications:** - Added a dedicated warning regarding 0 Hz (DC component) offsets. The library preserves the user's data structure but strongly recommends manual filtering of the DC component to prevent artificial static offsets in time-domain Pseudo Impulse Responses (PIRs).
  - Clarified the `values` output array in prediction band functions. Explicitly debunked the misconception that it holds the raw PIRs for every run, and advised users to utilize the tilde (`~`) operator to keep workspaces clean and save RAM.

### Fixed
- Fixed potential dimensional mismatch issues in `FRF_Plotbands` by enforcing strict validation prior to plotting.

---

### Links
* **GitHub Repository:** [https://github.com/VittorioFreiburg/FRF-statistics](https://github.com/VittorioFreiburg/FRF-statistics)
* **Reference Manual (Zenodo):** [https://zenodo.org/records/18924563](https://zenodo.org/records/18924563)