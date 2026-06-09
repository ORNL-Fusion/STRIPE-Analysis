% Bdot probe data analysis from 2020_06_09:
% =========================================================================

clear all
close all

% B_TANG -> Yellow -> S1

%% 1 - LOAD DATA:

% Load Bdot probe data:
% =========================================================================
b = load('Step_1_GetData_BdotProbe_2020_06_09.mat');

% Read spreadsheet:
% =========================================================================
c = load('Step_1_GetData_CapProbe_2020_06_09.mat');

%% 2 - CROP DATA:

% Crop data:
% =========================================================================
tStart = 4.15;
tEnd   = 4.7;

% bdot probe data:
for ii = 1:numel(b.VMAG_1)
    tdum = b.t_V{ii};
    rng = find(tdum>tStart & tdum<tEnd);

    VMAG_1{ii}  = b.VMAG_1{ii}(rng); 
    VMAG_2{ii}  = b.VMAG_2{ii}(rng); 
    VPH00_1{ii} = b.VPH00_1{ii}(rng); 
    VPH00_2{ii} = b.VPH00_2{ii}(rng); 
    VPH90_1{ii} = b.VPH90_1{ii}(rng);     
    VPH90_2{ii} = b.VPH90_2{ii}(rng);
    VREF_0{ii}  = b.VREF_0{ii}(rng);
    t_V{ii}     = b.t_V{ii}(rng);

    tdum = b.t_FWD{ii};
    rng = find(tdum>tStart & tdum<tEnd);

    b.REF{ii}    = b.REF{ii}(rng); 
    b.t_REF{ii}  = b.t_REF{ii}(rng); 
    b.FWD{ii}    = b.FWD{ii}(rng); 
    b.t_FWD{ii}  = b.t_FWD{ii}(rng); 
    
    tdum = b.t_isat{ii};
    rng = find(tdum>tStart & tdum<tEnd);
    
    b.isat{ii}   = b.isat{ii}(rng); 
    b.t_isat{ii} = b.t_isat{ii}(rng); 
end

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
b.rprobe = 7.5 - (b.r+1.5);
c.rprobe = c.r;

%% 4 - ABSOLUTE VALUES OF BDOT PROBE DATA:
% Extract absolute value of Reference signal S0 (VREF_0):
% =========================================================================
% Attenuation ratios:
att1 = 10^(-b.REF_IN_ATT/20);
att2 = 10^(-b.REF_SPLT_ATT/20);
att3 = 10^(-55/20);

% Calibration factors:
a0 = 0.0055;
b0 = 0.7645;

for ii = 1:numel(b.shotlist)
    S0{ii} = ((VREF_0{ii}-a0)/b0)/att2;
end

% Compute amplitude ratio R1 and R2:
% =========================================================================
for ii = 1:numel(b.shotlist)
    % Calibration factors:
    aa = 1.23;
    bb = 0.61;
    
    % S1:
    A = (VMAG_1{ii} - aa)/bb;
    R1{ii} = (10.^A); % R1 = S1/(S0*att5) where att5 = 10^(-AD8302_ATT/20)
        
    % S2:
    A = (VMAG_2{ii} - aa)/bb;
    R2{ii} = (10.^A); % R2 = S2/(S0*att) where att = 10^(-AD8302_ATT/20)
end

% Extract absolute values of S1 and S2:
% =========================================================================
att5 = 10^(-b.AD8302_ATT/20);
att6 = 10^(-b.BDOT_SPLT_ATT/20);

% Sensitivity of pick up coils:
sense_NORM = 1.58; % [mT/V]
sense_TANG = 2.56; % [mT/V]

for ii = 1:numel(b.shotlist)
    % The reference signal seen by the AD8302 is reduced by a factor att5
    S1{ii} = R1{ii}.*S0{ii}*att5;
    S2{ii} = R2{ii}.*S0{ii}*att5;
    
    % Bdot probe output voltages:
    V_BDOT_NORM{ii} = S1{ii}/att6;
    V_BDOT_TANG{ii} = S2{ii}/att6;
    
    % Absolute RF magnetic field at 13.56 MHz:
    B_NORM{ii} = V_BDOT_NORM{ii}*sense_NORM;
    B_TANG{ii} = V_BDOT_TANG{ii}*sense_TANG;
end

% Compute phase:
% =========================================================================
% calibration factors:
% Upper half:
a_UPP = -0.0248;
b_UPP = +0.0108;

% Lower half:
a_LOW = +3.89;
b_LOW = -0.0109;

for ii = 1:numel(S1) 
    % Loop for each element of V00:
    N = numel(VPH00_1{ii});   
    for jj = 1:N
       % Determine quadrant:
       if VPH90_1{ii}(jj)<1
           % Upper half
           phase_NORM{ii}(jj) = (VPH00_1{ii}(jj) - a_UPP)/b_UPP;
       else
           % Lower half
           phase_NORM{ii}(jj) = (VPH00_1{ii}(jj) - a_LOW)/b_LOW;
       end       
    end
end

for ii = 1:numel(S2) 
    % Loop for each element of V00:
    N = numel(VPH00_2{ii});   
    for jj = 1:N
       % Determine quadrant:
       if VPH90_2{ii}(jj)<1
           % Upper half
           phase_TANG{ii}(jj) = (VPH00_2{ii}(jj) - a_UPP)/b_UPP; % [deg]
       else
           % Lower half
           phase_TANG{ii}(jj) = (VPH00_2{ii}(jj) - a_LOW)/b_LOW; % [deg]
       end       
    end
end

%% 5 - TEST 180-DEGREE PHASE SHIFT:
% =========================================================================
n1 = find(b.shotlist == [30148]);
n2 = find(b.shotlist == [30149]);

% Tangent pick-up coil:
% -------------------------------------------------------------------------
% Arrow is pointing to Target or Dump, so we are sampling the "z" component
% of the RF magnetic field
saveFig = 1;
figureName = 'Step_2_Bz_180-degree Test';

% Magnitude:
figure('color','w')
subplot(2,1,1)
set(gca,'FontName','times','fontsize',11)
hold on
f1 = sgolay_t(B_TANG{n1},3,181);
f2 = sgolay_t(B_TANG{n2},3,181);
hM(1) = plot(t_V{n1},f1,'k');
hM(2) = plot(t_V{n2},f2,'r');
hM(3) = plot(b.t_FWD{n1},b.FWD{n1}*0.5,'g','LineWidth',2);
xlim([tStart,tEnd])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,1])
ylabel('$|\widetilde{B}_z|$ [mT]','interpreter','latex','Fontsize',13)
box on
grid on
L1 = legend(hM,'B-dot, 0 deg','B-dot, 180 deg','Reference: fwd power')
set(L1,'interpreter','latex','Fontsize',10)
textText = ['Shot# ',num2str(b.shotlist(n1)),' ,',num2str(b.shotlist(n2))];
hT = text(4.16,0.9,textText);
set(hT,'FontName','Times','FontSize',10)

% Phase:
subplot(2,1,2)
set(gca,'FontName','times','fontsize',11)
hold on
ph1 = phase_TANG{n1}*pi/180; % [Rad]
ph2 = phase_TANG{n2}*pi/180; % [Rad]
dphase = (ph1-ph2);

hP(1) = plot(t_V{n1},ph1/pi,'k.','markersize',3);
hP(2) = plot(t_V{n2},ph2/pi,'r.','markersize',3);
xlim([tStart,tEnd])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,2])
ylabel('$\theta/\pi$ [Rad]','interpreter','latex','Fontsize',13)
box on
grid on

if saveFig
    saveas(gcf,figureName,'tiffn')
end

% Normal pick-up coil:
% -------------------------------------------------------------------------
% We are sampling the "phi" component of the RF magnetic field
saveFig = 1;
figureName = 'Step_2_Bphi_180-degree Test';

% Magnitude:
figure('color','w')
subplot(2,1,1)
set(gca,'FontName','times','fontsize',11)
hold on
f1 = sgolay_t(B_NORM{n1},3,181);
f2 = sgolay_t(B_NORM{n2},3,181);
hM(1) = plot(t_V{n1},f1,'k');
hM(2) = plot(t_V{n2},f2,'r');
hM(3) = plot(b.t_FWD{n1},b.FWD{n1}*0.5,'g','LineWidth',2);
xlim([tStart,tEnd])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,1])
ylabel('$|\widetilde{B}_{\phi}|$ [mT]','interpreter','latex','Fontsize',13)
box on
grid on
L1 = legend(hM,'B-dot, 0 deg','B-dot, 180 deg','Reference: fwd power')
set(L1,'interpreter','latex','Fontsize',10)
textText = ['Shot# ',num2str(b.shotlist(n1)),' ,',num2str(b.shotlist(n2))];
hT = text(4.16,0.9,textText);
set(hT,'FontName','Times','FontSize',10)

% Phase:
subplot(2,1,2)
set(gca,'FontName','times','fontsize',11)
hold on
ph1 = phase_NORM{n1}*pi/180; % [Rad]
ph2 = phase_NORM{n2}*pi/180; % [Rad]
dphase = (ph1-ph2);
hP(1) = plot(t_V{n1},ph1/pi,'k.','markersize',3);
hP(2) = plot(t_V{n2},ph2/pi,'r.','markersize',3);
xlim([tStart,tEnd])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,2])
ylabel('$\theta/\pi$ [Rad]','interpreter','latex','Fontsize',13)
box on
grid on

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% 6 - COMPARE CAPACITIVE AND BDOT PROBE DATA:

saveFig = 1;
figureName = 'Step_2_CompareMagneticAndCapacitiveSignals';

figure('color','w')
hold on
% Plot RF:
n = 1;
plot(c.t_FWD{n},c.FWD{n},'g','LineWidth',2)
m = 5;
hf(1) = plot(b.t_FWD{m},b.FWD{m},'g','LineWidth',1)

% bdot cap data:
hf(2) = plot(t_V{m},sgolay_t(B_NORM{m},3,121),'r','LineWidth',1)

% Plot cap data:
hf(3) = plot(t_vrms{n},sgolay_t(vrms{n},3,121)/100,'k','LineWidth',2)

legend(hf,{'RF FWD',['B-dot, ',num2str(b.shotlist(m))],['Cap, ',num2str(c.shotlist(n))]})
ylim([0,1.5])
xlim([4.15,4.7])
box on
xlabel('t [sec]','Interpreter','latex','FontSize',12)
set(gca,'PlotBoxAspectRatio',[2.5 1 1])

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% 7 - ASSEMBLE TIME-RESOLVED RADIAL SCANS:

% Bdot probe data:
% =========================================================================
% From shot# 30149 to 30162:
b.shotSeries = b.shotlist(2:end);
[~,nn,~] = intersect(b.shotlist,b.shotSeries);
b.r_shotSeries = b.rprobe(nn);

% Time resolution of radial scan:
inc_t = 5e-3;
dt = mean(diff(t_V{1}));
dii = round(inc_t/dt);

% Assemble 2D arrays:
% -------------------------------------------------------------------------
% For all time:
for ii = 1:numel(t_V{1})-dii
    % For all radial positions:
    for ss = 1:numel(b.rprobe(nn))     
        rr = nn(ss);
        
        % Magnitude:
        Bn(ii,ss) = mean(B_NORM{rr}(ii:ii+dii));
        Bt(ii,ss) = mean(B_TANG{rr}(ii:ii+dii));
        t_B(ii)   = mean(t_V{rr}(ii:ii+dii));
        
        % Phase:
        pn(ii,ss) = mean(phase_NORM{rr}(ii:ii+dii))*pi/180; % [Rad]
        pt(ii,ss) = mean(phase_TANG{rr}(ii:ii+dii))*pi/180; % [Rad]
        t_p(ii)   = mean(t_V{rr}(ii:ii+dii));
        
    end
    pn(ii,:) = unwrap(pn(ii,:));
    pt(ii,:) = unwrap(pt(ii,:));
end
% Radial coordinate:
r_B = b.rprobe(nn);

% Shape of the 2D arrays is the following:
% Bn(time,radius)

% Interpolate 2D arrays:
% -------------------------------------------------------------------------
% Initial grid:
[rr,tt] = meshgrid(r_B,t_B);
% Final grid:
[RR,TT] = meshgrid(r_B,linspace(tStart,tEnd,500));
% Interpolate magnitude:
Bnorm.mag = interp2(rr,tt,Bn,RR,TT);
Btang.mag = interp2(rr,tt,Bt,RR,TT);
% Interpolate phase:
Bnorm.phase = interp2(rr,tt,pn,RR,TT);
Btang.phase = interp2(rr,tt,pt,RR,TT);
% Assign coordinates to data:
Bnorm.RR = RR;
Bnorm.TT = TT;
Btang.RR = RR;
Btang.TT = TT;

% Plot data:
% -------------------------------------------------------------------------
figure('color','w')
mesh(Bnorm.RR,Bnorm.TT,Bnorm.mag)
ylim([tStart,tEnd])
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
zlabel('$\widetilde{B}_r$ [mT]','interpreter','latex','Fontsize',12)
view([0,90])

figure('color','w')
mesh(Btang.RR,Btang.TT,Btang.mag)
ylim([tStart,tEnd])
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
zlabel('$\widetilde{B}_z$ [mT]','interpreter','latex','Fontsize',12)
view([0,90])

figure('color','w')
mesh(Bnorm.RR,Bnorm.TT,Bnorm.phase/pi)
ylim([tStart,tEnd])
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
zlabel('Arg($\widetilde{B}_r$)/$\pi$ [Rad]','interpreter','latex','Fontsize',12)
view([0,90])

figure('color','w')
mesh(Btang.RR,Btang.TT,Btang.phase/pi)
ylim([tStart,tEnd])
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
zlabel('Arg($\widetilde{B}_z$)/$\pi$ [Rad]','interpreter','latex','Fontsize',12)
view([0,90])

% Capacitive probe data:
% =========================================================================
% Organize shots in ascending order in "r"
nn = [10:-1:3,1];
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

%% 8 - PLOT TIME-RESOLVED RADIAL SCAN DATA:
% Example to compare Bdot probe and Capacitive probe 2D arrays:
% -------------------------------------------------------------------------

saveFig = 1;
figureName = 'Step2_MagAndCapa_2D_TimeResolved';

tRfStart = 4.166;

figure('color','w')
subplot(1,2,1)
hold on
contourf(Vrms.RR,Vrms.TT,Vrms.mag,50,'LineStyle','none')
line([0,4],[1,1]*tRfStart,[40,40],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([0.0,4.0])
title('$\widetilde{V}_{RMS}$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([0,90])
colorbar

subplot(1,2,2)
hold on
contourf(Btang.RR,Btang.TT,Btang.mag,50,'LineStyle','none')
line([-1,6],[1,1]*tRfStart,[40,40],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([-0.5,4.5])
title('$\widetilde{B}_z$ [mT]','interpreter','latex','Fontsize',12)
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
figureName = 'Step2_MagAndCapa_2D_TimeResolved_Xi';

% reference Magnetic flux at limiter:
B0 = 0.068;  % [T] 
r0 = 6.095;  %[cm]
phi0 = B0*r0^2;

% Magnetic field at probe locations:
b.B = 0.074;  % [T]
c.B = 0.133; % [T]

% Flux coordinates:
% Bnorm.Xi = b.B*sign(Bnorm.RR).*(Bnorm.RR.^2)/phi0; 
% Btang.Xi = b.B*sign(Btang.RR).*(Btang.RR.^2)/phi0; 
Bnorm.Xi = b.B*(Bnorm.RR.^2)/phi0; 
Btang.Xi = b.B*(Btang.RR.^2)/phi0; 
Vrms.Xi  = c.B*(Vrms.RR.^2)/phi0; 


figure('color','w')
subplot(1,2,1)
hold on
contourf(sqrt(Vrms.Xi),Vrms.TT,Vrms.mag,50,'LineStyle','none')
line([0,4],[1,1]*tRfStart,[40,40],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([0,1])
title('$\widetilde{V}_{RMS}$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([0,90])
colorbar

subplot(1,2,2)
hold on
contourf(sqrt(Btang.Xi),Btang.TT,Btang.mag,50,'LineStyle','none')
line([-1,6],[1,1]*tRfStart,[40,40],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([0,1])
title('$\widetilde{B}_z$ [mT]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([0,90])
colorbar

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% 9 - SAVE DATA:
% Assemble data package:
% =========================================================================

% Bdot probe data:
Bdot.Bnorm = Bnorm;
Bdot.Btang = Btang;
Bdot.shotlist = b.shotSeries;
Bdot.rprobe   = b.r_shotSeries;
Bdot.dateOfExperiment  = '2020_06_09';
Bdot.comment = 'Bdot probe data, MPEX-like limiter magnetic configuration, see shot summaries for 2020-06-09';
Bdot.t_rfStart = tRfStart;
Bdot.zLoc = 'Spool 2.5';

% Capacitive probe data:
Cap.Vrms = Vrms;
Cap.shotlist = c.shotSeries;
Cap.rprobe   = c.r_shotSeries;
Cap.dateOfExperiment  = '2020_06_09';
Cap.comment = 'Capacitive probe data, MPEX-like limiter magnetic configuration, see shot summaries for 2020-06-09';
Cap.t_rfStart = tRfStart;
Cap.zLoc = 'Spool 6.5';

% Assemble list of variables to save:
varlistSave = {'Bdot','Cap'};
save('Bdot_Cap_ProbeData_2020_06_09.mat',varlistSave{:})
