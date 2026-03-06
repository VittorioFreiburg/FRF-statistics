function [x,t] = FRF_pseudoimpulse(y,F,sf)
% FRF_PSEUDOIMPULSE Reconstructs time-domain signal from FRF components
% y  : Complex FRF values
% F  : Vector of frequencies (Hz)
% sf : Sampling frequency (Hz)

% Define tolerance for floating-point GCD calculation
tol = 1e-6; 

% Find the fundamental frequency (bf) using the local helper function
bf = numeric_gcd(F, tol);

% Fallback: if frequencies are truly irrational/incommensurate, bf approaches 0.
% We cap it to prevent infinite time vectors.
if bf < tol
    bf = min(F); 
end

% Define the time vector based on the period (1/bf) and sampling rate (1/sf)
step = 1/sf;
t = 0:step:(1/bf);

n = length(F);
x = zeros(size(t)); % Pre-allocate x for speed

% Manual Inverse Fourier Transform (Time-domain reconstruction)
for i = 1:n
    x = x + real(y(i))*cos(2*pi*F(i)*t) - imag(y(i))*sin(2*pi*F(i)*t);
end

end % <-- End of the main function


%% ========================================================================
%  LOCAL FUNCTIONS
%  ========================================================================

function g = numeric_gcd(numbers, tol)
% NUMERIC_GCD Finds the greatest common divisor of an array of real numbers
% within a specified floating-point tolerance.

    g = numbers(1);
    for i = 2:length(numbers)
        a = g;
        b = numbers(i);
        
        % Euclidean algorithm adapted for floating-point numbers
        while b > tol
            temp = mod(a, b);
            a = b;
            b = temp;
        end
        g = a;
    end
end