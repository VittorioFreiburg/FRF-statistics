# Statistical Analysis of Frequency Response Functions (FRF) - MATLAB Library

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) 


## Overview
This repository contains a comprehensive MATLAB library for the statistical analysis of Frequency Response Functions (FRFs). Performing statistics directly on complex numbers in the frequency domain is mathematically cumbersome and often non-intuitive. 

This library solves that problem by transforming complex frequency-domain data into time-domain Pseudo Impulse Responses (PIRs). Once in the time domain, it leverages non-parametric bootstrapping techniques to facilitate robust statistical evaluation without assuming a strict normal distribution of the data.

*Note: This repository is the official, actively maintained continuation of the FRF Statistics Library originally developed at Uniklinik Freiburg ( https://github.com/mcufidim/FRF-statistics )

## Key Features
* **Confidence & Prediction Bands:** Compute statistically rigorous boundaries for expected average behavior and individual future measurements.
* **Unpaired Group Comparisons:** Test for significant differences between independent groups of FRFs.
* **MIMO System Analysis:** Concatenate Multi-Input Multi-Output arrays into combined "Supervectors".
* **PERMANOVA Testing:** Perform Permutational Multivariate Analysis of Variance on complex FRF datasets.
* **Sample Evaluation:** Calculate exact empirical probabilities (PDF/CDF) and minimal encompassing prediction bands for individual test samples against a historical baseline.

## Installation
This library relies entirely on standard MATLAB functions and requires no external toolboxes. To install:
1. Clone or download this repository.
2. Add the folder to your MATLAB path:
```matlab
addpath('C:\Path\To\Your\Download\FRF_Statistics_Library');
savepath; % Optional: saves the path for future sessions
Quick Start
Here is a minimal working example to generate mock FRF data and calculate a 95% confidence band:

Matlab
% 1. Setup mock frequency vector and parameters
sf = 100;                 % Sampling frequency (Hz)
sample_time = 1/sf;       % Sample time
phi = 0.5:0.5:10;         % Frequencies from 0.5 to 10 Hz
N_samples = 15;           % Number of subjects/trials

% 2. Generate mock complex FRF data (N_samples x Length of phi)
mock_FRFs = complex(rand(N_samples, length(phi)) + 1, ...
                    rand(N_samples, length(phi)) - 0.5);

% 3. Calculate the 95% Confidence Band (alpha = 0.95, 1000 bootstraps)
alpha = 0.95;
B = 1000;
[avg, sigma, band, Cc] = FRF_ConfidenceBand(mock_FRFs, phi, sample_time, alpha, B);

% 4. Generate a time vector and plot
[~, time] = FRF_Supervector(mock_FRFs(1,:), phi, sample_time);
figure;
FRF_Plotbands(time, avg, band, [0 0.4 0.7]);
title('95% Confidence Band for Average PIR');
xlabel('Time (s)');
ylabel('Amplitude');
Example Scripts
To help you get started with real analysis pipelines, this repository includes three detailed example scripts:

SCRIPT_Example.m: Demonstrates standard confidence bands and unpaired tests for significant differences between two independent groups (e.g., Healthy vs. Patients).

SCRIPT_Example_MIMO.m: Expands the workflow to multidimensional MIMO arrays and executes a PERMANOVA test.

SCRIPT_MinBandExample.m: Focuses on evaluating a single test sample against a baseline, calculating encompassing bands and exact probabilities.

Documentation
A comprehensive Reference Manual (PDF) is included in this repository. It provides a didactic, step-by-step explanation of every function, its parameters, and the underlying statistical theory.

How to Cite
If you use this library in your research, please cite the manual/software repository directly:

Snippet di codice
@techreport{lippi2026frfmanual,
  title       = {Statistical Analysis of Frequency Response Functions: A MATLAB Library Reference Manual},
  author      = {Lippi, Vittorio},
  year        = {2026},
  month       = {March},
  type        = {Technical Report},
  institution = {Calejo},
  doi         = {10.5281/zenodo.18888167},
  url         = {https://doi.org/10.5281/zenodo.18888167}
}
Depending on the specific methodologies utilized from this library, please also consider citing the primary research papers where these algorithms were introduced:

Bootstrapping & PIR Transformation: Lippi, V. (2025). Bootstrap Prediction and Confidence Bands for Frequency Response Functions in Posturography. Experimental Techniques.

Unpaired Testing: Lippi, V. (2025). Unpaired Test for the Comparison of Frequency Response Functions Groups. arXiv preprint arXiv:2505.23778.

MIMO & PERMANOVA: Lippi, V., Johard, L., & Landucci, G. (2026). MIMO and Multi-Group Comparison Problem in Process Control... ICINCO (Under Review).

License
GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007