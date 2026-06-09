% Bdot probe data analysis from 2020_06_09:
% =========================================================================

clear all
close all

% Color coding:
% B_NORM -> Yellow -> S1
% B_TANG -> Red    -> S2

% Load data:
% =========================================================================
load('GetData_BdotProbe_2020_06_09.mat');

% Read spreadsheet:
% =========================================================================
% Read data from table where we can find RF power and orientation of probe
metadata = readtable('BdotProbeShotList_Final.xlsx','Sheet',1);

% Plasma radius:
% =========================================================================
% Based on radius of chamber and then using the 1.5 cm that the probe tip
% extends into the chamber:
rp = 7.5 - (r+1.5);

% Extract absolute value of Reference signal S0 (VREF_0):
% =========================================================================
% Attenuation ratios:
att1 = 10^(-REF_IN_ATT/20);
att2 = 10^(-REF_SPLT_ATT/20);
att3 = 10^(-55/20);

% Calibration factors:
a0 = 0.0055;
b0 = 0.7645;

for ii = 1:numel(shotlist)
    S0{ii} = ((VREF_0{ii}-a0)/b0)/att2;
end

% Compute amplitude ratio R1 and R2:
% =========================================================================
for ii = 1:numel(shotlist)
    % Calibration factors:
    a = 1.23;
    b = 0.61;
    
    % S1:
    A = (VMAG_1{ii} - a)/b;
    R1{ii} = (10.^A); % R1 = S1/(S0*att5) where att5 = 10^(-AD8302_ATT/20)
        
    % S2:
    A = (VMAG_2{ii} - a)/b;
    R2{ii} = (10.^A); % R2 = S2/(S0*att) where att = 10^(-AD8302_ATT/20)
    
    % Time traces:
    t_V{ii}    = t_V{ii}(1:end-1);
    t_isat{ii} = t_isat{ii}(1:end-1);
    t_FWD{ii}  = t_FWD{ii}(1:end-1);
    t_REF{ii}  = t_REF{ii}(1:end-1);
end

% Extract absolute values of S1 and S2:
% =========================================================================
att5 = 10^(-AD8302_ATT/20);
att6 = 10^(-BDOT_SPLT_ATT/20);

% Sensitivity of pick up coils:
sense_NORM = 1.58; % [mT/V]
sense_TANG = 2.56; % [mT/V]

for ii = 1:numel(shotlist)
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

%% 180 degree phase shift:
% =========================================================================
n1 = find(shotlist == [30148]);
n2 = find(shotlist == [30149]);

% Tangent pick-up coil:
% -------------------------------------------------------------------------
% Arrow is pointing to Target or Dump, so we are sampling the "z" component
% of the RF magnetic field
saveFig = 1;
figureName = 'Bz_180-degree Test';

% Magnitude:
figure('color','w')
subplot(2,1,1)
set(gca,'FontName','times','fontsize',11)
hold on
f1 = sgolay_t(B_TANG{n1},3,181);
f2 = sgolay_t(B_TANG{n2},3,181);
hM(1) = plot(t_V{n1},f1,'k');
hM(2) = plot(t_V{n2},f2,'r');
hM(3) = plot(t_V{n1},S0{n1}*0.6,'g');
xlim([4.13,4.7])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,1])
ylabel('$|\widetilde{B}_z|$ [mT]','interpreter','latex','Fontsize',13)
box on
grid on
L1 = legend(hM,'B-dot, 0 deg','B-dot, 180 deg','Reference: fwd power')
set(L1,'interpreter','latex','Fontsize',10)
textText = ['Shot# ',num2str(shotlist(n1)),' ,',num2str(shotlist(n2))];
hT = text(4.14,0.9,textText);
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
xlim([4.13,4.7])
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
figureName = 'Bphi_180-degree Test';

% Magnitude:
figure('color','w')
subplot(2,1,1)
set(gca,'FontName','times','fontsize',11)
hold on
f1 = sgolay_t(B_NORM{n1},3,181);
f2 = sgolay_t(B_NORM{n2},3,181);
hM(1) = plot(t_V{n1},f1,'k');
hM(2) = plot(t_V{n2},f2,'r');
hM(3) = plot(t_V{n1},S0{n1}*0.6,'g');
xlim([4.13,4.7])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,1])
ylabel('$|\widetilde{B}_{\phi}|$ [mT]','interpreter','latex','Fontsize',13)
box on
grid on
L1 = legend(hM,'B-dot, 0 deg','B-dot, 180 deg','Reference: fwd power')
set(L1,'interpreter','latex','Fontsize',10)
textText = ['Shot# ',num2str(shotlist(n1)),' ,',num2str(shotlist(n2))];
hT = text(4.14,0.9,textText);
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
xlim([4.13,4.7])
xlabel('time [s]','interpreter','latex','Fontsize',13)
ylim([0,2])
ylabel('$\theta/\pi$ [Rad]','interpreter','latex','Fontsize',13)
box on
grid on

if saveFig
    saveas(gcf,figureName,'tiffn')
end

%% Radial scan
close all

% =========================================================================
% PLOT ALL DATA TO DETERMINE STEADY-STATE SECTION:
% =========================================================================

% From shot# 30149 to 30162:
shotSeries = shotlist(2:end);
[~,nn,~] = intersect(shotlist,shotSeries);


tStart = 4.170;
tEnd   = 4.175;

% tStart = 4.19
% tEnd   = 4.20;

% tStart = 4.275;
% tEnd   = 4.35;
% 
% tStart = 4.4;
% tEnd   = 4.45;
% 
% 
% tStart = 4.5;
% tEnd   = 4.55;
% 
% tStart = 4.58;
% tEnd   = 4.63;

% Magnitude:
% -------------------------------------------------------------------------
figure
for ii = 1:numel(nn)
    kk = nn(ii);
    subplot(3,4,ii)
    hold on
    area([tStart,tEnd],[3,3])
    fN = sgolay_t(B_NORM{kk},3,181);
    fT = sgolay_t(B_TANG{kk},3,181);
    fR = sgolay_t(S0{kk},3,181);
    hN(ii) = plot(t_V{kk},fN,'k');
    hT(ii) = plot(t_V{kk},fT,'r');
    hT(ii) = plot(t_V{kk},fR,'g');
    plot(t_isat{kk},isat{kk}*5,'m')
    xlim([4.13,4.7])
    ylim([0,3])
    title(['r = ',num2str(rp(kk)),' [cm]'])
end

% Phase:
% -------------------------------------------------------------------------
figure
for ii = 1:numel(nn)
    kk = nn(ii);
    subplot(3,4,ii)
    hold on
    area([tStart,tEnd],[2,2])
    fR = sgolay_t(S0{kk},3,181);
    hN(ii) = plot(t_V{kk},phase_NORM{kk}/180,'k.','markersize',3);
    hT(ii) = plot(t_V{kk},phase_TANG{kk}/180,'r.','markersize',3);
%     hT(ii) = plot(t_V{kk},fR,'g');
    plot(t_isat{kk},isat{kk}*5,'m')
    xlim([4.13,4.7])
    ylim([0,2])
    title(['r = ',num2str(rp(kk)),' [cm]'])
end

% Choose steady-state values:
% -------------------------------------------------------------------------
for ii = 1:numel(nn)
    kk = nn(ii);

    rng = find(t_V{kk} > tStart & t_V{kk} < tEnd);
    % Magnitude:
    Bn(ii).mag = B_NORM{kk}(rng);
    Bt(ii).mag = B_TANG{kk}(rng);
    
    % Phase in [Rad]:
    Bn(ii).phase = phase_NORM{kk}(rng)*pi/180;
    Bt(ii).phase = phase_TANG{kk}(rng)*pi/180;
    
    % time:
    Bn(ii).t = t_V{kk}(rng);
    Bt(ii).t = t_V{kk}(rng);
    
    % radial location:
    Bn(ii).rp = rp(kk);
    Bt(ii).rp = rp(kk);
end

% =========================================================================
% PLOT THE SELECTED STEADY-STATE DATA:
% =========================================================================

% Magnitude:
% -------------------------------------------------------------------------
figure
subplot(1,2,1)
hold on
for kk = 1:numel(Bn)
    rr = ones(size(Bn(kk).mag));
    plot3(Bn(kk).t,Bn(kk).rp*rr,Bn(kk).mag)
end
view([30,30])

subplot(1,2,2)
hold on
for kk = 1:numel(Bn)
    rr = ones(size(Bt(kk).mag));
    plot3(Bt(kk).t,Bt(kk).rp*rr,Bt(kk).mag)
end
view([30,30])

% Produce steady-state values:
for kk = 1:numel(Bn)
    rprobe(kk) = Bn(kk).rp;
    Bn_ss.mag(kk)  = mean(Bn(kk).mag);
    Bn_ss.dmag(kk) = std(Bn(kk).mag,1);
    Bt_ss.mag(kk)  = mean(Bt(kk).mag);
    Bt_ss.dmag(kk) = std(Bt(kk).mag,1);
end

% Phase:
% -------------------------------------------------------------------------
figure

offset_1 = zeros(size(Bn,2),1);

subplot(1,2,1)
hold on
for kk = 1:numel(Bn)
    rr = ones(size(Bn(kk).phase));
    
    yy = Bn(kk).phase;
    
    if kk > 20
        yy = unwrap(Bn(kk).phase,pi);
    end
    
    yy = yy + offset_1(kk);
    plot3(Bn(kk).t,Bn(kk).rp*rr,yy/pi,'marker','.','markersize',3,'LineStyle','none')
end
view([70,10])

offset_2 = zeros(size(Bt,2),1);

subplot(1,2,2)
hold on
for kk = 1:numel(Bt)
    rr = ones(size(Bt(kk).phase));
    
    yy = Bt(kk).phase;
    
    if kk>20
        yy = unwrap(Bt(kk).phase,pi);
    end
    
    yy = yy + offset_2(kk);
    plot3(Bt(kk).t,Bt(kk).rp*rr,yy/pi,'marker','.','markersize',3,'LineStyle','none')
end
view([70,10])

% Produce steady-state values:
for kk = 1:numel(Bn)
    Bn_ss.phase(kk)  = mean(Bn(kk).phase);
    Bn_ss.dphase(kk) = std(Bn(kk).phase,1);
    Bt_ss.phase(kk)  = mean(Bt(kk).phase);
    Bt_ss.dphase(kk) = std(Bt(kk).phase,1);
end

% Define LCFS
rLCFS = 6.1; % [cm] based on magnetic field code

figure('color','w'); 
subplot(2,1,1)
hold on
hBn = errorbar(rprobe,Bn_ss.mag,Bn_ss.dmag);
hBt = errorbar(rprobe,Bt_ss.mag,Bt_ss.dmag);
hL  = line([1,1]*rLCFS,[0,3]);
set(gca,'FontName','Times','FontSize',11)
set(hBn,'color','k','lineWidth',2,'marker','sq')
set(hBt,'color','r','lineWidth',2,'marker','o')
set(hL,'color','g','lineWidth',2)
box on
grid on
set(hBn,'marker','sq')
xlim([-1,6.5])
ylabel('$|\widetilde{B}|$ [mT]','Interpreter','Latex','FontSize',13)
xlabel('r [cm]','Interpreter','Latex','FontSize',13)
legendText{1} = ['$\widetilde{B}_r$'];
legendText{2} = ['$\widetilde{B}_z$'];
legendText{3} = ['LUFS'];
hLeg = legend([hBn,hBt,hL],legendText);
set(hLeg,'interpreter','Latex','FontSize',11,'Location','NorthEast')

hold on
box on
grid on

subplot(2,1,2)
hold on
hPn = errorbar(rprobe,unwrap(Bn_ss.phase)/pi,Bn_ss.dphase/pi);
hPt = errorbar(rprobe,unwrap(Bt_ss.phase)/pi,Bt_ss.dphase/pi);
hL  = line([1,1]*rLCFS,[-4,3]);
box on
grid on
set(gca,'FontName','Times','FontSize',11)
set(hPn,'color','k','lineWidth',2,'marker','sq')
set(hPt,'color','r','lineWidth',2,'marker','o')
set(hL,'color','g','lineWidth',2)
xlim([-1,6.5])
ylim([-4,2])
ylabel('Arg($\widetilde{B}$)/$\pi$ [Rad]','Interpreter','Latex','FontSize',13)
xlabel('r [cm]','Interpreter','Latex','FontSize',13)

shotnum = ['shots: ',num2str(shotSeries(1)),'-',num2str(shotSeries(end))];
text(3.5,-3,shotnum,'FontSize',10,'Interpreter','latex')

timeRng = ['time: ',num2str(tStart),'-',num2str(tEnd),' [s]'] ;
text(3.5,-2,timeRng,'FontSize',10,'Interpreter','latex')
