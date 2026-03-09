function [avg,sigma,band,Cc,chist,values] = FRF_ConfidenceBandDifference(FRF1,FRF2,phi,sample_time,alpha,B,Bs)
% FRF_CONFIDENCEBANDDIFFERENCE Computes confidence bands on the difference 
% between the means of two independent groups (Unpaired Test).
%
%   [avg, sigma, band, Cc, chist, values] = FRF_ConfidenceBandDifference(FRF1, FRF2, phi, dt, alpha, B, Bs)
%   evaluates whether two distinct sets of FRFs are statistically different 
%   from one another. It calculates the average difference between the groups' 
%   PIRs and generates a confidence band around that difference. If the resulting 
%   band does not cross the zero line, the difference is statistically significant.
%
%   INPUTS:
%       FRF1        : Complex frequency response data for Group 1 (N1 x F).
%       FRF2        : Complex frequency response data for Group 2 (N2 x F).
%                     Note: N1 must equal N2.
%       phi         : Vector of evaluated frequencies (in Hz).
%       sample_time : Time step resolution (dt) for internal PIR conversion.
%       alpha       : The desired statistical significance level (e.g., 0.95).
%       B           : The number of outer bootstrap repetitions used to calculate 
%                     the main statistic distribution.
%       Bs          : The number of inner bootstrap repetitions used to estimate 
%                     the standard deviation for studentization.
%
%   OUTPUTS:
%       avg         : The difference between the average PIRs of the two groups.
%       sigma       : The measure of variation (standard deviation of the difference).
%       band        : A 2-row matrix defining the statistical boundaries:
%                     - Row 1: The Lower Bound (Floor)
%                     - Row 2: The Upper Bound (Ceiling)
%       Cc          : The threshold constant obtained by the bootstrap.
%       chist       : The cumulative histogram of the bootstrap statistics.
%       values      : The corresponding boundary values for the histogram.
%
%   NOTES:
%       Ensure your data does not contain unwanted 0 Hz DC offsets before processing.

N1 = size(FRF1, 1);
N2 = size(FRF2, 1);

% --- INPUT VALIDATION ---
if N1 ~= N2
    error('FRF_ConfidenceBandDifference:UnequalSampleSizes', ...
        'FRF1 and FRF2 must contain the same number of samples (N1 = %d, N2 = %d).', N1, N2);
end
N = N1; % Use a single N for both sets moving forward

% --- UNIVERSAL DATA INGESTION ---
% Replace the old for-loops with the universal engine.
% y1 and y2 will be N x L matrices (where L is time steps * channels).
[y1, ~] = FRF_Supervector(FRF1, phi, sample_time);
[y2, ~] = FRF_Supervector(FRF2, phi, sample_time);

ns = size(y1, 2); % Total number of concatenated features

xm = mean(y1, 1) - mean(y2, 1);

%% --- VECTORIZED INITIAL STD ESTIMATE ---
idx_s1 = randi(N, N * Bs, 1);
idx_s2 = randi(N, N * Bs, 1);

mean_s1 = reshape(mean(reshape(y1(idx_s1, :), N, Bs, ns), 1), Bs, ns);
mean_s2 = reshape(mean(reshape(y2(idx_s2, :), N, Bs, ns), 1), Bs, ns);

sx = std(mean_s1 - mean_s2, 0, 1);

%% --- MAIN BOOTSTRAP LOOP ---
STAT = zeros(1, B);

for b = 1:B
    % 1. Resample outer groups
    yb1 = y1(randi(N, 1, N), :);
    yb2 = y2(randi(N, 1, N), :);

    xb = mean(yb1, 1) - mean(yb2, 1);

    % 2. Vectorized Inner Bootstrap
    idx_b1 = randi(N, N * Bs, 1);
    idx_b2 = randi(N, N * Bs, 1);

    mean_ybb1 = reshape(mean(reshape(yb1(idx_b1, :), N, Bs, ns), 1), Bs, ns);
    mean_ybb2 = reshape(mean(reshape(yb2(idx_b2, :), N, Bs, ns), 1), Bs, ns);

    % 3. Calculate inner standard deviation
    sb = std(mean_ybb1 - mean_ybb2, 0, 1);

    % 4. Calculate studentized max deviation
    STAT(b) = max(abs(xm - xb) ./ sb);
end

%% Histogram
STAT = sort(STAT);
[chist, values] = histcounts(STAT, 1000, 'Normalization', 'cdf');
Cc = values(find(chist > alpha, 1, 'first'));

avg = xm;
sigma = sx;

band = [avg - Cp*sigma; avg + Cp*sigma]; % Row 1: Lower bound, Row 2: Upper bound

end