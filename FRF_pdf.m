function [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(varargin)
% [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(X,FRFs,phi,sample_time,B,metric)
%
% estimates the cumulative densitiy function and the density function
% associated to the sample X and a STD on their estimation. 
% the input METRIC defines the measure used to define the distance
% between X and the mean of the sample FRFs. By default it is the sum of
% squared residuals
% distance, if METRIC is a function handle it is applied directly. The
% following strings can be specified:
% - 'squared' sum of squared residuals
% - 'max' maximum difference between two samples
% 
    function [cdf, pdf, sigma_cdf, sigma_pdf] = FRF_pdf(varargin)
% [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(X,FRFs,phi,sample_time,B,metric)
%
% Estimates the cumulative density function and the probability density 
% function associated with the sample X, and a STD on their estimation.
%  
%
% Fully supports both SISO and MIMO arrays via FRF_Supervector.

if nargin == 5
   metric = @(x,y) sum((x-y).^2, 2);
elseif nargin == 6
    metric = varargin{6};
    if isa(metric, 'string') || isa(metric, 'char')
        if metric == "squared"
            metric = @(x,y) sum((x-y).^2, 2);
        elseif metric == "max"
            metric = @(x,y) max(abs(x-y), [], 2);
        else
            error([char(metric), ' is not a valid metric']);
        end
    elseif ~isa(metric, 'function_handle')
        error('METRIC must be a string or a function handle');
    end
else
    error('Input arguments must be 5 or 6');
end

X = varargin{1};
FRFs = varargin{2};
phi = varargin{3};
sample_time = varargin{4};
B = varargin{5};

% --- UNIVERSAL DATA INGESTION ---
% xt will be 1 x L (where L is time steps * channels)
[xt, ~] = FRF_Supervector(X, phi, sample_time);

% yt will be N x L
[yt, ~] = FRF_Supervector(FRFs, phi, sample_time);

N = size(yt, 1); % number of samples
Ds = max(1, fix(N/20)); % Smoothing window for derivative

% --- BOOTSTRAP LOOP ---
STAT = zeros(1, B);
dSTAT = zeros(1, B);

for b = 1:B 
    % Resample with replacement
    resamp = randi(N, 1, N);
    yb = yt(resamp, :);
   
    xb = mean(yb, 1);
    
    % Vectorized metric calculation
    % implicitly handles L features (SISO or concatenated MIMO)
    es = sort(metric(yb, xb)); 
    et = metric(xt, xb);
    
    idx = find(es > et, 1, 'first');
    if isempty(idx)
        idx = N;
    end
    STAT(b) = idx / N;
        
    % Smoothing window for PDF estimation
    i1 = idx - Ds;
    i2 = idx + Ds;
    
    if i1 < 1
        i1 = 1;
        i2 = 1 + Ds;
    end
    
    if i2 > N
        i2 = N;
        i1 = N - Ds;
    end
    
    % Protect against division by zero if es(i2) == es(i1)
    diff_es = es(i2) - es(i1);
    if diff_es == 0
        diff_es = 1e-10; % Small epsilon
    end
    
    dSTAT(b) = Ds / (N * diff_es);
end

cdf = mean(STAT);
sigma_cdf = std(STAT);

pdf = mean(dSTAT);
sigma_pdf = std(dSTAT);

end
