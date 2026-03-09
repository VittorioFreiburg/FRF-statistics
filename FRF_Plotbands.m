function h = FRF_Plotbands(time, avg, band, color)
%
%   h = FRF_Plotbands(time, avg, band, color) plots the average Pseudo 
%   Impulse Response (PIR) as a solid line and its corresponding statistical 
%   boundaries (confidence or prediction bands) as dash-dotted lines.
%
%   INPUTS:
%       time  : A vector containing the time steps for the x-axis.
%       avg   : A vector containing the average PIR values. Must be the 
%               exact same length as 'time'.
%       band  : A 2xL matrix containing the boundary values, where L is 
%               the length of 'time'. Row 1 must be the lower bound and 
%               Row 2 must be the upper bound.
%       color : A 1x3 RGB array (e.g., [0 0.4 0.7] for blue) defining the 
%               base color of the plot. The average line will automatically 
%               be plotted in a slightly darker/muted shade of this color.
%
%   OUTPUTS:
%       h     : The handle to the current figure (gcf).
%
%   EXAMPLE:
%       % After computing a 95% confidence band:
%       [avg, sigma, band, Cc] = FRF_ConfidenceBand(FRFs, phi, dt, 0.95, 1000);
%       
%       figure;
%       FRF_Plotbands(t, avg, band, [0 0.4 0.7]); % Plot in blue
%       xlabel('Time (s)'); ylabel('Amplitude');
%       title('System Response with 95% Confidence Band');
%
%   NOTES:
%       - This function automatically applies 'hold on' to the current figure.
%       - It includes an internal safety check to ensure 'time' and 'avg' 
%         dimensions match, preventing cryptic MATLAB indexing errors.

% --- INPUT VALIDATION ---
if length(time) ~= length(avg)
    error('FRF Error: The time vector and average vector must be the exact same length.');
end

h = gcf;
plot(time, avg, 'color', color/2, 'linewidth', 2); hold on;
plot(time,band(1,:),'color',color,'linewidth',1.5,'linestyle','-.');
plot(time,band(2,:),'color',color,'linewidth',1.5,'linestyle','-.');