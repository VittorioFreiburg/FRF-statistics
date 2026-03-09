function [cdf,pdf,sigma_cdf,sigma_pdf] = FRF_pdf(varargin)
% FRF_PDF Estimates the empirical density functions associated with a sample.
%
%   [cdf, pdf, sigma_cdf, sigma_pdf] = FRF_pdf(X, FRFs, phi, sample_time, B, metric)
%   calculates the exact statistical probability of observing a specific sample (X) 
%   given a baseline dataset (FRFs). It estimates both the Cumulative Density 
%   Function (CDF) and Probability Density Function (PDF) alongside their 
%   respective standard deviations.
%
%   INPUTS:
%       X           : The single complex FRF test sample (1 x F).
%       FRFs        : The baseline complex FRF data (N x F).
%       phi         : Vector of evaluated frequencies (in Hz).
%       sample_time : Time step resolution (dt).
%       B           : Number of bootstrap iterations.
%       metric      : (Optional) Defines the distance measure between X and the 
%                     sample mean. Can be a custom function handle or a string:
%                     - 'squared' : Sum of squared residuals (Default).
%                     - 'max'     : Maximum difference between two samples.
%
%   OUTPUTS:
%       cdf         : The estimated cumulative density probability of X.
%       pdf         : The estimated probability density of X.
%       sigma_cdf   : Standard deviation on the CDF estimation.
%       sigma_pdf   : Standard deviation on the PDF estimation.

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
