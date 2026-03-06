%% SCRIPT_Example_MIMO.m
% A didactic example of analyzing Multi-Input Multi-Output (MIMO) systems
% using Pseudo-Impulse Responses (PIRs) and PERMANOVA.
%
% THE SCENARIO:
% Imagine we are testing human balance. A platform tilts, and we measure
% two things simultaneously for each subject:
%   Channel 1: Body Sway (how much they wobble)
%   Channel 2: Muscle Activity (EMG of the calf muscles)
%
% Because we are measuring two things at once, this is a MIMO system!
% We want to compare three groups of 15 subjects to see if their balance
% strategies are statistically different.

clear; clc; close all;

%% 1. Setup and Synthetic Data Generation
disp('1. Generating synthetic MIMO data for 3 groups...');

f1 = [0.05 0.15 0.30 0.40 0.55 0.70 0.90 1.10 1.35 1.75 2.20]; % Frequencies (Hz)
sf = 22;               % Sampling frequency (Hz)
sample_time = 1 / sf;  % Sample time (s)
N = 15;                % Number of subjects per group
M = 2;                 % Number of measurement channels (Sway and EMG)

% Pre-allocate 3D arrays: (Subjects x Frequencies x Channels)
FRF_G1 = zeros(N, length(f1), M); % Group 1: Healthy A
FRF_G2 = zeros(N, length(f1), M); % Group 2: Pathological (Different dynamics)
FRF_G3 = zeros(N, length(f1), M); % Group 3: Healthy B (Same as G1, just different noise)

% Define the "True" underlying dynamics for our channels
% (Using simple complex transfer functions for demonstration)
Sway_Healthy   = 1.0 ./ (1 + 1i * (f1 / 0.5));
EMG_Healthy    = 0.5 ./ (1 + 1i * (f1 / 1.0));

Sway_Pathology = 1.8 ./ (1 + 1i * (f1 / 0.3)); % Larger sway, different phase
EMG_Pathology  = 0.5 ./ (1 + 1i * (f1 / 1.0)); % Muscle response stays the same

% Generate the subjects by adding random biological noise
noise_level = 0.15;
for i = 1:N
    % Group 1: Healthy A
    FRF_G1(i, :, 1) = Sway_Healthy   + noise_level * (randn(1,11) + 1i*randn(1,11));
    FRF_G1(i, :, 2) = EMG_Healthy    + noise_level * (randn(1,11) + 1i*randn(1,11));

    % Group 2: Pathological
    FRF_G2(i, :, 1) = Sway_Pathology + noise_level * (randn(1,11) + 1i*randn(1,11));
    FRF_G2(i, :, 2) = EMG_Pathology  + noise_level * (randn(1,11) + 1i*randn(1,11));

    % Group 3: Healthy B (Same underlying math as Group 1, just new random noise)
    FRF_G3(i, :, 1) = Sway_Healthy   + noise_level * (randn(1,11) + 1i*randn(1,11));
    FRF_G3(i, :, 2) = EMG_Healthy    + noise_level * (randn(1,11) + 1i*randn(1,11));
end

%% 1.5 Visualize the generated MIMO FRFs in Frequency Domain
% Before running complex statistics, let's visually inspect the data.
% We plot the Amplitude (Modulus) and Phase for both channels.
% This shows us exactly how the Pathological group differs physically.
disp('1.5 Visualizing the Frequency Response Functions (Amplitude and Phase)...');

figure('Name', 'MIMO System: Frequency Domain', 'Position', [150 150 1000 700]);

% Define distinct colors for the 3 groups
color_G1 = [0, 0.4470, 0.7410]; % Blue (Healthy A)
color_G2 = [0.8500, 0.3250, 0.0980]; % Red (Pathological)
color_G3 = [0.4660, 0.6740, 0.1880]; % Green (Healthy B)

% -------------------------------------------------------------------------
% COLUMN 1: CHANNEL 1 (BODY SWAY)
% -------------------------------------------------------------------------

% Subplot 1: Amplitude of Body Sway
subplot(2, 2, 1); hold on;
plot(f1, abs(FRF_G3(:,:,1)), 'Color', color_G3);
plot(f1, abs(FRF_G1(:,:,1)), 'Color', color_G1);
plot(f1, abs(FRF_G2(:,:,1)), 'Color', color_G2);
title('Channel 1: Body Sway (Amplitude)');
ylabel('Gain (Modulus)');
xlim([min(f1) max(f1)]);
grid on;

% Subplot 3: Phase of Body Sway
% Note: We use recompactUp to unwrap the 2*pi phase jumps for a clean plot
subplot(2, 2, 3); hold on;
plot(f1, recompactUp(angle(FRF_G3(:,:,1))), 'Color', color_G3);
plot(f1, recompactUp(angle(FRF_G1(:,:,1))), 'Color', color_G1);
plot(f1, recompactUp(angle(FRF_G2(:,:,1))), 'Color', color_G2);
title('Channel 1: Body Sway (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
xlim([min(f1) max(f1)]);
grid on;

% -------------------------------------------------------------------------
% COLUMN 2: CHANNEL 2 (MUSCLE EMG)
% -------------------------------------------------------------------------

% Subplot 2: Amplitude of Muscle EMG
subplot(2, 2, 2); hold on;
plot(f1, abs(FRF_G3(:,:,2)), 'Color', color_G3);
plot(f1, abs(FRF_G1(:,:,2)), 'Color', color_G1);
plot(f1, abs(FRF_G2(:,:,2)), 'Color', color_G2);
title('Channel 2: Muscle EMG (Amplitude)');
ylabel('Gain (Modulus)');
xlim([min(f1) max(f1)]);
grid on;

% Subplot 4: Phase of Muscle EMG
subplot(2, 2, 4); hold on;
plot(f1, recompactUp(angle(FRF_G3(:,:,2))), 'Color', color_G3);
plot(f1, recompactUp(angle(FRF_G1(:,:,2))), 'Color', color_G1);
plot(f1, recompactUp(angle(FRF_G2(:,:,2))), 'Color', color_G2);
title('Channel 2: Muscle EMG (Phase)');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
xlim([min(f1) max(f1)]);
grid on;

% --- Add a Clean Legend ---
% We create "dummy" lines with thick widths just for the legend so it looks clean
subplot(2, 2, 2);
h1 = plot(NaN, NaN, 'Color', color_G1, 'LineWidth', 2);
h2 = plot(NaN, NaN, 'Color', color_G2, 'LineWidth', 2);
h3 = plot(NaN, NaN, 'Color', color_G3, 'LineWidth', 2);
legend([h1, h2, h3], {'Group 1 (Healthy A)', 'Group 2 (Pathological)', 'Group 3 (Healthy B)'}, 'Location', 'northeast');

%% 2. Pairwise Comparison: Healthy vs Pathological
% Let's ask: "Are Group 1 and Group 2 statistically different?"
% We calculate the Confidence Band of the DIFFERENCE between their means.
% If the difference band excludes zero, the groups are physically different!
disp('2. Running Pairwise Test: Group 1 vs Group 2 (Should find a difference)...');

alpha = 0.95; % We want to be 95% confident
B = 500;      % Outer bootstrap iterations
Bs = 50;      % Inner bootstrap iterations for standard deviation

% Run the test
[avg_1v2, ~, band_1v2, ~, ~, ~] = FRF_ConfidenceBandDifference(FRF_G1, FRF_G2, f1, sample_time, alpha, B, Bs);

% We need the time vector for plotting. We can get it quickly using FRF_Supervector.
% Remember, the supervector stitches Channel 1 and Channel 2 end-to-end!
% 1. Get the time vector for a SINGLE channel
[~, t_single] = FRF_Supervector(FRF_G1(1,:,:), f1, sample_time);

% 2. Create an artificial "stitched" time vector for the supervector plot
dt = t_single(2) - t_single(1); % Time step
% Shift the second channel so it starts right after the first one ends
t_super = [t_single, t_single + t_single(end) + dt];

figure(); set(gcf, 'Name', 'Pairwise Comparisons', 'Position', [100 100 800 600]);

subplot(2,1,1);
plot(t_super, avg_1v2, 'k', 'LineWidth', 2); hold on;
plot(t_super, band_1v2(1,:), 'b-.', 'LineWidth', 1.5);
plot(t_super, band_1v2(2,:), 'b-.', 'LineWidth', 1.5);
yline(0, 'r--', 'LineWidth', 1.5); % The ZERO line

% Find the split point between Channel 1 and Channel 2 for visualization
midpoint = t_super(end) / 2;
xline(midpoint, 'k:', 'LineWidth', 2);
text(midpoint/2, max(band_1v2(1,:))*0.8, 'Body Sway', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(midpoint + midpoint/2, max(band_1v2(1,:))*0.8, 'Muscle EMG', 'HorizontalAlignment', 'center', 'FontSize', 12);

title('Group 1 (Healthy) vs Group 2 (Pathological)');
ylabel('Difference (Units)');
legend('Mean Difference', '95% Confidence Band', '', 'Zero Line');
% Notice how the blue band does NOT touch the red zero line in the "Body Sway" section!


%% 3. Pairwise Comparison: Healthy A vs Healthy B
% Now let's ask: "Are Group 1 and Group 3 statistically different?"
% Since they are both healthy, the math should confidently tell us they are the same.
disp('3. Running Pairwise Test: Group 1 vs Group 3 (Should NOT find a difference)...');

[avg_1v3, ~, band_1v3, ~, ~, ~] = FRF_ConfidenceBandDifference(FRF_G1, FRF_G3, f1, sample_time, alpha, B, Bs);

subplot(2,1,2);
plot(t_super, avg_1v3, 'k', 'LineWidth', 2); hold on;
plot(t_super, band_1v3(1,:), 'b-.', 'LineWidth', 1.5);
plot(t_super, band_1v3(2,:), 'b-.', 'LineWidth', 1.5);
yline(0, 'r--', 'LineWidth', 1.5); % The ZERO line
xline(midpoint, 'k:', 'LineWidth', 2);

title('Group 1 (Healthy A) vs Group 3 (Healthy B)');
xlabel('Concatenated Time (s)');
ylabel('Difference (Units)');
% Notice how the red zero line sits comfortably INSIDE the blue band everywhere!
% This means the difference between them is statistically zero.


%% 4. Global Multi-Group Comparison (PERMANOVA)
% If we had 10 groups, doing pairwise tests for all of them would inflate our
% risk of a false positive (finding a difference by pure random luck).
% PERMANOVA solves this by doing ONE global test first[cite: 83].
disp('4. Running Global PERMANOVA across all 3 groups...');

% Combine all data into one giant dataset (45 subjects x 11 frequencies x 2 channels)
FRF_ALL = cat(1, FRF_G1, FRF_G2, FRF_G3);

% Create labels so the math knows who is in which group
% 1 = Healthy A, 2 = Pathological, 3 = Healthy B
group_labels = [ones(N,1); 2*ones(N,1); 3*ones(N,1)];

% Run PERMANOVA with 999 permutations
% (PERMANOVA uses Euclidean distance, so the supervectors are automatically
% Z-score standardized internally so Body Sway doesn't drown out the Muscle EMG [cite: 84, 87]).
B_perm = 999;
[p_value, F_obs, R2] = FRF_permanova(FRF_ALL, group_labels, f1, sample_time, B_perm);

disp('--------------------------------------------------');
disp('                PERMANOVA RESULTS                 ');
disp('--------------------------------------------------');
fprintf('Pseudo-F Statistic : %.2f\n', F_obs);
fprintf('Effect Size (R^2)  : %.2f%% (Variance explained by the groups)\n', R2 * 100);
fprintf('p-value            : %.4f\n', p_value);

if p_value < 0.05
    disp('CONCLUSION: The global test is SIGNIFICANT (p < 0.05).');
    disp('At least one group has a completely different balance strategy!');
else
    disp('CONCLUSION: The global test is NOT SIGNIFICANT.');
    disp('All groups share the same underlying balance strategy.');
end
disp('--------------------------------------------------');