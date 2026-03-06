function [avg,sigma,band,Cc,chist,values] = FRF_ConfidenceBandDifference(FRF1,FRF2,phi,sample_time,alpha,B,Bs)
%[AVG,SIGMA,band,CC,CHIST,VALUES] =
%FRF_CONFIDENCEBANDDIFFERENCE(FRFS1,FRFS2,PHI,SAMPLE_TIME,B)
% Confidence bands on the difference between the means of two groups FRF1
% andd FRF2.
%
% B is the number of bootstrap repetitions
% Bs is the number of bootstrap repetitions used to estimate STD
%
% avg is the difference between average PIRs of the groups, sigma is the
% measure of variation ?ˆx(t), band a two row matrix with the boundaries of
% the band. Cp is the threshold constant obtained by the bootstrap.  FRFS
% is a matrix where each row represents a FRF of the set, phi is the vector
% of frequencies  and SAMPLE_TIME is the sample time of the PIRs. Chist is
% a vector representing the cumulative histogram for the values  returned
% in VALUES.

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

band = [avg + Cc * sigma; avg - Cc * sigma];

end