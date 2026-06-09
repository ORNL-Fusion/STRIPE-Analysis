% Cap probe data analysis from 2020_06_17:
% =========================================================================

clear all
close all

%% 1 - LOAD DATA:

% Load data:
% =========================================================================
c = load('Step_1_GetData_CapProbe_2020_06_17.mat');

%% 2 - CROP DATA:

% Crop data:
% =========================================================================
tStart = 4.15;
tEnd   = 4.7;

% cap probe data:
for ii = 1:numel(c.vrms)
    tdum = c.t_vrms{ii};
    rng = find(tdum>tStart & tdum<tEnd);

    vrms{ii}   = c.vrms{ii}(rng); 
    t_vrms{ii} = c.t_vrms{ii}(rng);

    tdum = c.t_FWD{ii};
    rng = find(tdum>tStart & tdum<tEnd);

    c.REF{ii}    = c.REF{ii}(rng); 
    c.t_REF{ii}  = c.t_REF{ii}(rng); 
    c.FWD{ii}    = c.FWD{ii}(rng); 
    c.t_FWD{ii}  = c.t_FWD{ii}(rng); 
    
    tdum = c.t_isat{ii};
    rng = find(tdum>tStart & tdum<tEnd);
    
    c.isat{ii}   = c.isat{ii}(rng); 
    c.t_isat{ii} = c.t_isat{ii}(rng);  
end

%% 3 - DEFINE RADIAL POSITION OF PROBES:
% =========================================================================
% Based on radius of chamber and then using the 1.5 cm that the probe tip
% extends into the chamber:
c.rprobe = c.r;

%% 4 - ASSEMBLE TIME-RESOLVED RADIAL SCANS:

% Capacitive probe data:
% =========================================================================
% Organize shots in ascending order in "r"
% N = numel(c.shotlist);
% [~,nn] = sort(c.rprobe);
[~,nn] = unique(c.rprobe);
c.shotSeries   = c.shotlist(nn);
c.r_shotSeries = c.rprobe(nn);

% Time resolution of radial scan:
inc_t = 5e-3;
dt = mean(diff(t_vrms{1}));
dii = round(inc_t/dt);

% Assemble 2D arrays:
% -------------------------------------------------------------------------
% For all time:
for ii = 1:numel(t_vrms{1})-dii
    % For all radial positions:
    for ss = 1:numel(c.rprobe(nn))     
        rr = nn(ss);
        
        % Magnitude:
        vrms_m(ii,ss) = mean(vrms{rr}(ii:ii+dii));
        t_Vrms_m(ii)  = mean(t_vrms{rr}(ii:ii+dii));
        
    end
end
% Radial coordinate:
r_Vrms_m = c.rprobe(nn);

% Shape of the 2D arrays is the following:
% vrms_m(time,radius)

% Interpolate 2D arrays:
% -------------------------------------------------------------------------
% Initial grid:
[rr,tt] = meshgrid(r_Vrms_m,t_Vrms_m);
% Final grid:
[RR,TT] = meshgrid(r_Vrms_m,linspace(tStart,tEnd,500));
% Interpolate magnitude:
Vrms.mag = interp2(rr,tt,vrms_m,RR,TT);
% Assign coordinates to data:
Vrms.RR = RR;
Vrms.TT = TT;

% Plot data:
% -------------------------------------------------------------------------
figure('color','w')
surf(Vrms.RR,Vrms.TT,Vrms.mag,'LineStyle','none')
ylim([tStart,tEnd])
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
zlabel('$\widetilde{V}_{RMS}$ [V]','interpreter','latex','Fontsize',12)
view([0,90])

%% 5 - PLOT TIME-RESOLVED RADIAL SCAN DATA:
% Example to compare Bdot probe and Capacitive probe 2D arrays:
% -------------------------------------------------------------------------

saveFig = 1;
figureName = 'Step_2_Cap_2D_TimeResolved';

tRfStart = 4.166;

figure('color','w')
hold on
contourf(Vrms.RR,Vrms.TT,Vrms.mag,50,'LineStyle','none')
line([0,13],[1,1]*tRfStart,[40,40],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([0.0,13])
title('$\widetilde{V}_{RMS}$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([0,90])
colorbar

if saveFig
    saveas(gcf,figureName,'tiffn')
end

% Plot as a function of flux coordinate:
% -------------------------------------------------------------------------

saveFig = 1;
figureName = 'Step_2_Cap_2D_TimeResolved_Xi';

% Reference Magnetic flux at limiter:
B0 = 0.0544;  % [T] 
r0 = 6.256;  %[cm]
phi0 = B0*r0^2;

% Magnetic field at probe locations:
c.B = 0.183; % [T]

% Flux coordinates:
Vrms.Xi  = c.B*(Vrms.RR.^2)/phi0; 

figure('color','w')
hold on
contourf(sqrt(Vrms.Xi),Vrms.TT,Vrms.mag,50,'LineStyle','none')
hL = line([0,2],[1,1]*tRfStart,[40,40],'Color','r','LineWidth',3)
hXi = line([1,1],[0,5],[0.1,0.1],'Color','r','LineWidth',2,'LineStyle','--');
ylim([4.16,4.25])
xlim([0,2])
title('$\widetilde{V}_{RMS}$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([0,90])
colorbar

hLeg = legend([hL,hXi],'Start of RF','$\chi$ = 1');
hLeg.Location = 'north';
hLeg.Interpreter = 'latex';

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% 9 - SAVE DATA:
% Assemble data package:
% =========================================================================

% Capacitive probe data:
Cap.Vrms = Vrms;
Cap.shotlist = c.shotSeries;
Cap.rprobe   = c.r_shotSeries;
Cap.dateOfExperiment  = '2020_06_17';
Cap.comment = 'Capacitive probe data, Window-limiter magnetic configuration, see shot summaries for 2020-06-17';
Cap.t_rfStart = tRfStart;
Cap.zLoc = 'Spool 6.5';

% Assemble list of variables to save:
varlistSave = {'Cap'};
save('Cap_ProbeData_2020_06_17.mat',varlistSave{:})
