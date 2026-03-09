function [p_value, F_obs, R2, perm_F] = FRF_permanova(FRFs, groups, phi, sample_time, B)
% FRF_PERMANOVA Performs Permutational Multivariate Analysis of Variance on FRFs.
%
%   [p_value, F_obs, R2, perm_F] = FRF_permanova(FRFs, groups, phi, sample_time, B)
%   tests whether the variance between multiple groups of FRFs significantly 
%   exceeds the variance within the groups. It operates on the time-domain 
%   Supervectors and automatically applies Z-Score Standardization to ensure 
%   variables with large magnitudes do not unfairly dominate the analysis.
%
%   INPUTS:
%       FRFs        : Complex frequency response data (supports MIMO arrays).
%       groups      : A vector of group labels corresponding to each sample in FRFs.
%       phi         : Vector of evaluated frequencies (in Hz).
%       sample_time : Time step resolution (dt).
%       B           : Number of permutations used to build the null distribution.
%
%   OUTPUTS:
%       p_value     : The statistical significance (probability of observing the 
%                     variance by random chance). Typically, p < 0.05 rejects 
%                     the null hypothesis.
%       F_obs       : The observed pseudo-F statistic for the true groupings.
%       R2          : The R-squared value, indicating the proportion of variance 
%                     explained by the grouping.
%       perm_F      : The vector of pseudo-F statistics generated across all 
%                     B permutations (useful for plotting the null distribution).

    % 1. Get the raw stacked PIR Supervectors
    [S_raw, ~] = FRF_Supervector(FRFs, phi, sample_time);
    
    % 2. Z-Score Standardization (CRITICAL FOR PERMANOVA)
    % As per the methodology, this ensures variables with large magnitudes 
    % (e.g., pressure) don't overpower small ones (e.g., composition).
    S_mean = mean(S_raw, 1);
    S_std  = std(S_raw, 0, 1);
    S_std(S_std == 0) = 1; % Prevent division by zero
    
    S = (S_raw - S_mean) ./ S_std;
    
    % 3. Group Setup
    [N, L] = size(S);
    [G_labels, ~, G_idx] = unique(groups);
    G = length(G_labels); 
    
    df_A = G - 1;
    df_W = N - G;
    
    % 4. Calculate Total Sum of Squares (SS_T)
    grand_centroid = mean(S, 1);
    SS_T = sum(sum((S - grand_centroid).^2));
    
    % 5. Helper Function: Within-Group SS (SS_W)
    function ssw = calc_SSW(labels)
        ssw = 0;
        for g = 1:G
            group_samples = S(labels == g, :);
            group_centroid = mean(group_samples, 1);
            ssw = ssw + sum(sum((group_samples - group_centroid).^2));
        end
    end
    
    % 6. Calculate Observed Statistics
    SS_W_obs = calc_SSW(G_idx);
    SS_A_obs = SS_T - SS_W_obs;
    
    F_obs = (SS_A_obs / df_A) / (SS_W_obs / df_W);
    R2 = SS_A_obs / SS_T;
    
    % 7. Permutation Testing
    perm_F = zeros(B, 1);
    for b = 1:B
        shuffled_idx = G_idx(randperm(N));
        SS_W_perm = calc_SSW(shuffled_idx);
        SS_A_perm = SS_T - SS_W_perm;
        perm_F(b) = (SS_A_perm / df_A) / (SS_W_perm / df_W);
    end
    
    p_value = sum(perm_F >= F_obs) / B;
end