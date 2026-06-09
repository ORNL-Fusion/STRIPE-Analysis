function [IntensityMean,rngPlasma,t0] = GetDataDuringPlasma(f,t_f,options)
% To determine the start and end of the RF pulse we use the following
% method:
% For each frame, we calculate the mean pixel intensity for the entire image
% This operation results in a 1D vector called "intensityRaw_mean"
% We use a diff operation to find the starting and ending edges of
% "intensityRaw_mean", these edges are called n1 and n2
% We apply an offset to n1 and n2 to widen the time window
% From these edges we define a range which we use to extract a subset of
% data from the main raw data

% Offsets to widen the RF pulse window:
n_Before = options.n_Before;
n_After = options.n_After;

% Define start and end of pulse and define data range:
[Nx,Ny,Nz] = size(f); 

% For each frame, compute the mean intensity over the entire frame:
for ii = 1:Nz
   IntensityMean(ii) = mean(mean(f(1:20:Nx,1:20:Ny,ii))); 
end
    
% Find the start and end of the RF pulse:
[~,n1] = max(diff(IntensityMean,1));
[~,n2] = min(diff(IntensityMean,1));
    
% Start of the RF pulse relative to the raw data's time trace
t0 = t_f(n1);
    
% Define the time window over which the RF is on:
% Include some points before and after the RF
rngPlasma = [n1-n_Before:n2+n_After];

% Get data only during the plasma shot:
g   = f(:,:,rngPlasma);
t_g = t_f(rngPlasma) - t0;

end

