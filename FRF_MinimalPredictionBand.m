function [avg,sigma,band,Cp,chist,values,alpha] = FRF_MinimalPredictionBand(X,FRFs,phi,sample_time,B)
% FRF_MINIMALPREDICTIONBAND Computes the minimum threshold to encompass a tested FRF.
%
%   [avg, sigma, band, Cp, chist, values, alpha] = FRF_MinimalPredictionBand(X, FRFs, phi, dt, B)
%   evaluates how abnormal a specific tested sample (X) is compared to a baseline 
%   population (FRFs). It calculates the tightest possible boundary level (alpha) 
%   needed so that the tested sample just barely fits inside the band.
%
%   INPUTS:
%       X           : The single complex FRF test sample to be evaluated (1 x F).
%       FRFs        : The baseline/historical complex FRF data (N x F).
%       phi         : Vector of evaluated frequencies (in Hz).
%       sample_time : Time step resolution (dt) for internal PIR conversion.
%       B           : Number of bootstrap iterations.
%
%   OUTPUTS:
%       avg         : The baseline average PIR.
%       sigma       : The measure of baseline variation (standard deviation).
%       band        : A 2-row matrix containing the boundaries of the band 
%                     (Row 1: Lower bound, Row 2: Upper bound) scaled precisely to X.
%       Cp          : The minimal threshold constant required to include X.
%       chist       : The cumulative histogram for the values returned in 'values'.
%       values      : The threshold values corresponding to the histogram.
%       alpha       : The empirical probability level associated with encompassing X. 
%                     (e.g., an alpha of 0.99 means X is highly abnormal).

    % --- UNIVERSAL DATA INGESTION ---
    % xt will be 1 x L (where L is time steps * channels)
    [xt, ~] = FRF_Supervector(X, phi, sample_time);
    
    % yt will be N x L
    [yt, t] = FRF_Supervector(FRFs, phi, sample_time);

    N = size(yt, 1);

    % Baseline statistics
    xm = mean(yt, 1);
    sx = std(yt, 0, 1);
    
    % Safety: Prevent division by zero if a channel has zero variance
    sx_safe = sx;
    sx_safe(sx_safe == 0) = 1;

    %% --- VECTORIZED BOOTSTRAP LOOP ---
    STAT = zeros(N, B);

    for b = 1:B 
        % Resample with replacement
        resamp = randi(N, 1, N);
        yb = yt(resamp, :);
        
        xb = mean(yb, 1);
        sb = std(yb, 0, 1);
        sb(sb == 0) = 1; % Safety against zero variance in bootstrap sample
        
        % Vectorized calculation of studentized max deviation 
        % Computes all N deviations simultaneously across all L features
        STAT(:, b) = max(abs(yt - xb) ./ sb, [], 2);
    end

    %% --- HISTOGRAM & THRESHOLD ---
    % Flatten and sort the statistic
    STAT = sort(STAT(:)');
    [chist, values] = histcounts(STAT, 1000, 'Normalization', 'cdf'); 

    % Find the minimal Cp required to encompass the specific test sample X
    Cp = max(abs(xt - xm) ./ sx_safe);

    % Find the corresponding empirical alpha level
    alpha = chist(find(values > Cp, 1, 'first'));
    
    % Edge case: if the test sample is so extreme it exceeds all bootstrap values
    if isempty(alpha)
        alpha = 1; 
    end

    avg = xm;
    sigma = sx;

    % The band automatically scales to physical units (SISO) or raw concatenated units (MIMO)
    band = [avg - Cp*sigma; avg + Cp*sigma]; % Row 1: Lower bound, Row 2: Upper bound
end