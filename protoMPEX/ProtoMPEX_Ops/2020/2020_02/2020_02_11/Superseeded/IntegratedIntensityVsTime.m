function [flag] = IntegratedIntensityVsTime(f,options)
% To determine the start and end of the RF pulse we use the following
% method:
% For each frame, we calculate the mean pixel intensity for the entire image
% This operation results in a 1D vector called "intensityRaw_mean"
% We use a diff operation to find the starting and ending edges of
% "intensityRaw_mean", these edges are called n1 and n2
% We apply an offset to n1 and n2 to widen the time window
% From these edges we define a range which we use to extract a subset of
% data from the main raw data
% 
% Offsets to reject data from the entire time trace
n1_offset = options.n1_offset;
n2_offset = options.n2_offset;
% Offsets to widen the RF pulse window
n_Before = options.n_Before;
n_After = options.n_After;
frameRate = options.frameRate;
dt = 1/frameRate;

% Define start and end of pulse and define data range
[Nx,Ny,Nz] = size(f); 

for ii = 1:Nz
   f_mean(ii) = mean(mean(f(1:20:Nx,1:20:Ny,ii))); 
end
t_f_mean = 0:dt:(Nz-1)*dt;
    
% Find the start and end of the RF pulse:
[~,n1] = max(diff(f_mean(n1_offset:end),1));
n1 = n1 + n1_offset - 1;
[~,n2] = min(diff(f_mean(n2_offset:end),1));
n2 = n2 + n2_offset - 1;
    
% Start of the RF pulse relative to the raw data's time trace
t0 = t_f_mean(n1);
    
% Define the time window over which the RF is on:
% Include some points before and after the RF
rng = [n1-n_Before:n2+n_After];

% Plot data to confirm correctness of calculation
if 1
figure; 
set(gcf,'Tag','DefineDataRng')
hold on
f_mean_noOffset = f_mean-min(f_mean);
plot(t_f_mean-t0,f_mean_noOffset,'LineWidth',0.5);
hIR = plot(t_f_mean(rng)-t0,f_mean_noOffset(rng),'k.','LineWidth',2);
title('intensityRaw_mean_noOffset','Interpreter','none')
ylabel('{\Delta}intensity')
xlabel('t [s]')
xlim([-0.5,2])
box on; set(gcf,'color','w')
end

return
% Extract the relevent data from "intensityRaw":
intensity = f(:,:,rng);
t_intensity = 0:dt:(length(rng)-1)*dt;


end

