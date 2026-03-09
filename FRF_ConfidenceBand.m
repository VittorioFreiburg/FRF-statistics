function [avg,sigma,band,Cc,chist,values] = FRF_ConfidenceBand(FRFs,phi,sample_time,alpha,B)
% FRF_CONFIDENCEBAND Calculates a non-parametric bootstrap confidence band.
%
%   [avg, sigma, band, Cc, chist, values] = FRF_ConfidenceBand(FRFs, phi, dt, alpha, B)
%   estimates the true average behavior of a system and generates statistical 
%   boundaries (confidence bands) indicating where the population mean lies. 
%   It utilizes a resampling bootstrap technique on the time-domain Pseudo 
%   Impulse Responses (PIRs).
%
%   INPUTS:
%       FRFs        : Complex frequency response data (N samples x F frequencies).
%                     Fully supports 3D matrices for MIMO arrays.
%       phi         : Vector of evaluated frequencies (in Hz).
%       sample_time : Time step resolution (dt) for the internal PIR conversion.
%       alpha       : The desired statistical significance level (e.g., 0.95 
%                     for a 95% confidence band).
%       B           : Number of bootstrap iterations (e.g., 1000). Higher 
%                     values increase boundary precision but take longer to compute.
%
%   OUTPUTS:
%       avg         : The mean Pseudo Impulse Response across all samples.
%       sigma       : The standard deviation of the samples.
%       band        : A 2-row matrix defining the statistical boundaries:
%                     - Row 1: The Lower Bound (Floor)
%                     - Row 2: The Upper Bound (Ceiling)
%       Cc          : The threshold multiplier constant derived from the bootstrap.
%       chist       : The cumulative histogram of the bootstrap statistics.
%       values      : The corresponding boundary values for the histogram.
%
%   EXAMPLE:
%       % Compute a 95% confidence band using 1000 bootstraps
%       [avg, ~, band, ~] = FRF_ConfidenceBand(my_data, freqs, 0.01, 0.95, 1000);
%
%   NOTES:
%       This function relies on FRF_Supervector to convert frequency data 
%       into the time domain prior to computing the statistics. Ensure your 
%       data does not contain unwanted 0 Hz DC offsets.

N=size(FRFs,1); %number of FRFs

[yt, t] = FRF_Supervector(FRFs, phi, sample_time);

    xm=mean(yt);
    sx=std(yt);


% Pre-allocate
STAT = zeros(1, B);

for b = 1:B 
    % 1. Resample with replacement
    resamp = randi(N, 1, N);
    yb = yt(resamp, :);
    
    % 2. Calculate Bootstrap statistics
    xb = mean(yb, 1);
    sb = std(yb, 0, 1);
    
    % 3. Calculate maximum deviation of the means
    STAT(b) = max(abs(xm - xb) ./ sb);
end


%% Histogram

STAT=sort(STAT);
%H=histogram(STAT,1000,'Normalization','cdf','Edgecolor','none')
[chist, values] = histcounts(STAT,1000,'Normalization','cdf'); 
Cc=values(find(chist>alpha,1,'first'));

avg=xm;
sigma=sx;

band = [avg - Cp*sigma; avg + Cp*sigma]; % Row 1: Lower bound, Row 2: Upper bound

end

