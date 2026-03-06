function [avg,sigma,band,Cp,chist,values,alpha] = FRF_MinimalPredictionBand(X,FRFs,phi,sample_time,B)
%[AVG,SIGMA,VOL,CP,CHIST,VALUES,ALPHA] =
%FRF_MINIMALPREDICTIONBAND(X,FRFS,PHI,SAMPLE_TIME,B)
% Computes the minimum Cp that includes the tested FRF X and the empirical
% ALPHA associated with it.
% where avg is average PIR, sigma is the measure of variation ?ˆx(t), band
% a two row matrix with the boundaries of the band. Cp is the threshold
% constant obtained by the bootstrap.  FRFS is a matrix where each row
% represents a FRF of the set, phi is the vector of frequencies  and
% SAMPLE_TIME is the sample time of the PIRs. Chist is a vector
% representing the cumulative histogram for the values  returned in VALUES.
% 
% Fully supports both SISO and MIMO arrays via FRF_Supervector.

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
    band = [avg + Cp * sigma; avg - Cp * sigma];
end