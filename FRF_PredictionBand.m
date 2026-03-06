function [avg,sigma,band,Cp,chist,values] = FRF_PredictionBand(FRFs,phi,sample_time,alpha,B)
%[AVG,SIGMA,VOL,CP,CHIST,VALUES] =
%FRF_PREDICTIONBAND(FRFS,PHI,SAMPLE_TIME,B)
% where avg is average PIR, sigma is the measure of variation ?ˆx(t), band
% a two row matrix with the boundaries of the band. Cp is the threshold
% constant obtained by the bootstrap.  FRFS is a matrix where each row
% represents a FRF of the set, phi is the vector of frequencies  and
% SAMPLE_TIME is the sample time of the PIRs. Chist is a vector
% representing the cumulative histogram for the values  returned in VALUES.

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

band=[avg+Cp*sigma;avg-Cp*sigma];
end

