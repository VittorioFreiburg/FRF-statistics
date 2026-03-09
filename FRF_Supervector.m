function [S, t] = FRF_Supervector(FRFs, phi, sample_time)
% FRF_SUPERVECTOR Converts MIMO Frequency Response Functions into raw PIRs.
%
%   [S, t] = FRF_Supervector(FRFs, phi, sample_time) translates complex
%   frequency-domain data into time-domain Pseudo Impulse Responses (PIRs). 
%   For Multi-Input Multi-Output (MIMO) systems, it concatenates the channels 
%   horizontally into a single, flattened "Supervector" to allow for 
%   simultaneous multivariate statistical analysis.
%
%   INPUTS:
%       FRFs        : The complex frequency response data. This can be:
%                     1) A 2D matrix (N x F) for standard SISO data.
%                     2) A 3D matrix (N x F x M) for MIMO data with shared parameters.
%                     3) A 1 x M cell array of 2D matrices for heterogeneous MIMO data.
%                     (N = samples/subjects, F = frequencies, M = channels)
%       phi         : Vector of frequencies (in Hz) corresponding to the columns 
%                     of the FRFs. Can be a cell array for heterogeneous data.
%       sample_time : The desired time step (dt) for the output time signal 
%                     (i.e., 1 / sampling_frequency).
%
%   OUTPUTS:
%       S           : The resulting Supervector matrix (N x L), where N is the 
%                     number of samples and L is the total concatenated time 
%                     steps across all channels.
%       t           : The generated time vector (column vector) corresponding 
%                     to the time steps in S. If inputs are heterogeneous, 
%                     this will be a cell array of time vectors.
%
%   IMPORTANT NOTE ON DC OFFSETS (0 Hz):
%       Time-domain transformations are highly sensitive to 0 Hz components. 
%       If your 'phi' vector starts at 0 Hz, it may introduce an artificial 
%       static offset (DC bias) in your resulting PIRs. If this is undesired, 
%       it is highly recommended that you manually remove the 0 Hz column 
%       from both 'FRFs' and 'phi' before calling this function.

% 1. Input Normalization (same as before)
if ~iscell(FRFs)
    [N, ~, M] = size(FRFs);
    FRF_cell = cell(1, M);
    for m = 1:M, FRF_cell{m} = FRFs(:, :, m); end
    return_single_t = true;
else
    M = length(FRFs);
    N = size(FRFs{1}, 1);
    FRF_cell = FRFs;
    return_single_t = false;
end

if ~iscell(phi), phi_cell = repmat({phi}, 1, M); else, phi_cell = phi; return_single_t = false; end
if isscalar(sample_time), st_vec = repmat(sample_time, 1, M); else, st_vec = sample_time; return_single_t = false; end

% 2. Generate and Stack PIRs
PIR_cell = cell(1, M);
t_cell = cell(1, M);

for m = 1:M
    FRF_m = FRF_cell{m};
    sf_m = 1 / st_vec(m);

    [x1, t_m] = FRF_pseudoimpulse(FRF_m(1, :), phi_cell{m}, sf_m);
    T_m = length(t_m);
    t_cell{m} = t_m;

    PIR_m = zeros(N, T_m);
    PIR_m(1, :) = x1;

    for i = 2:N
        [x, ~] = FRF_pseudoimpulse(FRF_m(i, :), phi_cell{m}, sf_m);
        PIR_m(i, :) = x;
    end

    % Just store the raw physical PIRs!
    PIR_cell{m} = PIR_m;
end

% 3. Construct Supervector by simple horizontal concatenation
S = cell2mat(PIR_cell);

if return_single_t
    t = t_cell{1}(:); % Force column vector
else
    % Force all individual time vectors in the cell array to be columns
    t = cellfun(@(v) v(:), t_cell, 'UniformOutput', false);
end
end