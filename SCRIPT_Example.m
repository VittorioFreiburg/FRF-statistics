%% SCRIPT_Example.m
% This script demonstrates how to use the FRF Statistics Library to compute
% and visualize Bootstrap-based Prediction Bands for Frequency Response Functions.
% It covers data loading, phase unwrapping, bootstrap thresholding, and 
% time-domain visualization (Phase Impulse Responses).

clear; clc; close all;

%% 1. Parameters & Setup
% Moving these to the top makes the script much easier to use as a template.
f1 = [0.05 0.15 0.30 0.40 0.55 0.70 0.90 1.10 1.35 1.75 2.20]; % Frequencies (Hz)
sf = 22;               % Sampling frequency (Hz)
sample_time = 1 / sf;  % Sample time (s)
alpha = 0.95;          % Confidence level (95%)
B = 1000;              % Number of bootstrap repetitions

%% 2. Load Data
% SET1: The mean/reference group used to build the prediction band
% FRF_All: The test samples we want to evaluate against the band
load('SET1.mat');
load('FRF_All.mat');

% Combine test FRFs into a single matrix for plotting
FRF_test = [FRF_DEC; FRF_IC; FRF_EM]; 

% We will use SET1 to build our statistical model
SET = SET1; 

%% 3. Frequency Domain Visualization (Bode Plot)
figure(1);
set(gcf, 'Name', 'Frequency Domain: Amplitude and Phase', 'Position', [100 100 600 800]);

% --- SUBPLOT 1: AMPLITUDE ---
subplot(2,1,1);
% Plot mean (SET) in grey
plot(f1, abs(SET), 'Color', [0.5 0.5 0.5]); hold on;
% Plot test samples in color
plot(f1, abs(FRF_test), 'LineWidth', 1.5);
xlim([0.05 2.2]);
ylabel('Gain (°/°)');
title('EC FRFs (Body sway)');

% Apply custom colors for the test sets
K = [1 0 1; 0 1 1; 1 0.5 0.5]; 
colororder(K);

% --- SUBPLOT 2: PHASE ---
subplot(2,1,2);
% Phase Unwrapping: angle() maps to [-pi, pi]. recompactUp fixes 2*pi jumps.
PH_set = recompactUp(angle(SET')');
plot(f1, PH_set, 'Color', [0.5 0.5 0.5]); hold on;

PH_test = recompactUp(angle(FRF_test')');
plot(f1, PH_test, 'LineWidth', 1.5);
xlim([0.05 2.2]);
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
colororder(K);

% Save Figure
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', screenposition(3:4));
print(gcf, 'DataPrediction', '-dpdf', '-r0');


%% 4. Bootstrap Statistical Analysis
disp('Computing Prediction Band via Bootstrap...');
% This function maps the FRFs to the time domain, resamples them B times, 
% and identifies the statistical threshold (Cp) that encompasses alpha% of the variations.
[avg, sigma, band, Cp, chist, values] = FRF_PredictionBand(SET, f1, sample_time, alpha, B);


%% 5. Compute Time-Domain Signals for Plotting
% Convert the mean SET into Phase Impulse Responses (PIRs)
ntrials = size(SET, 1);
% ns will be dynamically determined by the first call to FRF_pseudoimpulse
[x_temp, t] = FRF_pseudoimpulse(SET(1,:), f1, sf);
ns = length(t);
yt = zeros(ntrials, ns);

for i = 1:ntrials 
    [x, ~] = FRF_pseudoimpulse(SET(i,:), f1, sf);
    yt(i,:) = x;
end


%% 6. Plot the Statistical Histogram
figure(2);
set(gcf, 'Name', 'Bootstrap CDF and Threshold');

% Plot the Cumulative Distribution Function (CDF) of the Bootstrap statistics
plot(values(1:end-1), chist, 'LineWidth', 2); hold on;

% Label the axes using LaTeX for mathematical symbols
xlabel('$C_c$', 'Interpreter', 'latex');
ylabel('$\max \left( \left| \hat{x}_b - \hat{x} \right| / \hat{\sigma}_b \right)$', 'Interpreter', 'latex', 'Fontsize', 11);

% Add threshold lines
yline(alpha, 'LineWidth', 1.5, 'LineStyle', '-.', 'Color', [1 0 0], 'Label', '\alpha = 95%', 'Fontsize', 10);
xline(Cp, 'LineWidth', 1.5, 'LineStyle', '-.', 'Color', [1 0 1], 'Label', 'C_p', 'Fontsize', 10, 'LabelVerticalAlignment', 'middle');

% Save Figure
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', screenposition(3:4));
print(gcf, 'HistogramPrediction', '-dpdf', '-r0');


%% 7. Time Domain Prediction Band Visualization
figure(3);
set(gcf, 'Name', 'Time Domain Prediction Band');

UPPER = band(1,:);
LOWER = band(2,:);

% Plot signals (grey)
plot(t, yt, 'Color', [0.8 0.8 0.8]); hold on;

% Plot mean (black) and 95% bands (blue dashed)
plot(t, avg, 'Color', [0 0 0], 'LineWidth', 1.5);
plot(t, UPPER, 'Color', [0.1 0.1 1], 'LineWidth', 1.5, 'LineStyle', '-.');
plot(t, LOWER, 'Color', [0.1 0.1 1], 'LineWidth', 1.5, 'LineStyle', '-.');
yline(0, 'k-.');

% Convert and plot the test samples to see if they fall outside the band
ntrials_test = size(FRF_test, 1);
yx = zeros(ntrials_test, ns);
for i = 1:ntrials_test
    [x, ~] = FRF_pseudoimpulse(FRF_test(i,:), f1, sf);
    yx(i,:) = x;
end
plot(t, yx, 'LineWidth', 1.5);
colororder([0 1 1; 1 0.5 0.5; 1 0 1]); % Custom colors for test samples

% Clean up the plot
legend({'Mean $x_i(t)$', '$\hat{x}(t)$', 'Upper $95\%$', 'Lower $95\%$'}, 'Interpreter', 'latex', 'Location', 'best');
xlabel('Time (s)');
ylabel('Pseudo Impulse Response (°)');
title('95% Prediction Band in Time Domain');

% Save Figure
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], 'PaperSize', screenposition(3:4));
print(gcf, 'TimeDomainPrediction', '-dpdf', '-r0');

disp('Analysis Complete.');