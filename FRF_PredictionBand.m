function [avg,sigma,band,Cp,chist,values] = FRF_PredictionBand(FRFs,phi,sample_time,alpha,B)
% FRF_PREDICTIONBAND Calculates a non-parametric prediction band for future samples.
%
%   [avg, sigma, band, Cp, chist, values] = FRF_PredictionBand(FRFs, phi, dt, alpha, B)
%   computes a statistical boundary indicating where the next *single* observation 
%   of a system is expected to fall. Unlike a confidence band (which bounds the 
%   average), a prediction band accounts for the natural spread of individual 
%   data points and is therefore inherently wider.
%
%   INPUTS:
%       FRFs        : Complex frequency response baseline data (N x F).
%       phi         : Vector of evaluated frequencies (in Hz).
%       sample_time : Time step resolution (dt) for internal PIR conversion.
%       alpha       : The desired probability level for the band (e.g., 0.95 
%                     means a 95% chance the next single sample falls inside).
%       B           : Number of bootstrap iterations (e.g., 1000).
%
%   OUTPUTS:
%       avg         : The average Pseudo Impulse Response (PIR).
%       sigma       : The standard deviation of the samples.
%       band        : A 2-row matrix defining the statistical boundaries:
%                     - Row 1: The Lower Bound (Floor)
%                     - Row 2: The Upper Bound (Ceiling)
%       Cp          : The threshold constant obtained by the bootstrap.
%       chist       : The cumulative histogram of the bootstrap statistics.
%       values      : The corresponding boundary values for the histogram.

N=size(FRFs,1); %number of FRFs
[yt, t] = FRF_Supervector(FRFs, phi, sample_time);

    xm=mean(yt);
    sx=std(yt);


STAT=zeros(1,B*N);


STAT = zeros(N, B);

for b = 1:B
    % 1. Resample with replacement
    resamp = randi(N, 1, N);
    yb = yt(resamp, :);
    
    % 2. Calculate Bootstrap statistics (specifying dimension 1 for safety)
    xb = mean(yb, 1);
    sb = std(yb, 0, 1);
    
    % 3. VECTORIZED DEVIATION CALCULATION
    % MATLAB automatically expands xb and sb to match the dimensions of yt.
    % max(..., [], 2) finds the max value across the time samples (columns) 
    % for each of the N individual rows simultaneously.
    STAT(:, b) = max(abs(yt - xb) ./ sb, [], 2); 
end

% Flatten the matrix into a 1x(B*N) vector for the histogram
STAT = STAT(:)';

%% Histogram

STAT=sort(STAT);
%H=histogram(STAT,1000,'Normalization','cdf','Edgecolor','none')
[chist, values] = histcounts(STAT,1000,'Normalization','cdf'); 
Cp=values(find(chist>alpha,1,'first'));

avg=xm;
sigma=sx;

band = [avg - Cp*sigma; avg + Cp*sigma]; % Row 1: Lower bound, Row 2: Upper bound
end

