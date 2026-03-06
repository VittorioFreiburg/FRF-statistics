function [avg,sigma,band,Cc,chist,values] = FRF_ConfidenceBand(FRFs,phi,sample_time,alpha,B)
%[AVG,SIGMA,band,CC,CHIST,VALUES] =
%FRF_CONFIDENCEBAND(FRFS,PHI,SAMPLE_TIME,B)
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

band=[avg+Cc*sigma;avg-Cc*sigma];

end

