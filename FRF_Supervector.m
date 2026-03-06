function [S, t] = FRF_Supervector(FRFs, phi, sample_time)
% FRF_SUPERVECTOR Converts MIMO FRFs into concatenated raw PIRs.
% 
% S : N x L matrix of raw PIRs (concatenated horizontally for MIMO)
% t : Time vector (or cell array of time vectors for heterogeneous inputs)

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
    
    if return_single_t, t = t_cell{1}; else, t = t_cell; end
end