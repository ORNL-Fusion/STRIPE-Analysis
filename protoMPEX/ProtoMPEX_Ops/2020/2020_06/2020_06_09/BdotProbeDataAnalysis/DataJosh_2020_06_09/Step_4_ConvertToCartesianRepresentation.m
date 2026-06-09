% Demonstrate how to read data from structure and convert it to cartesian
% representation of the phasor

clear all
close all

% Load data:
% =========================================================================
load('Bdot_Cap_ProbeData_2020_06_09.mat')
% Experimental data is represented in polar coordinates: magnitude and
% phase
% We want to convert this into its cartesian representation which provides
% the real and imaginary parts:


% The following is the "guts" of this script!!

% Convert to cartesian representation:
% =========================================================================
% =========================================================================
% Projection of the phasor into the real axis:
real_Bz = Bdot.Btang.mag.*cos(Bdot.Btang.phase);
% Projection of the phasor into the imaginary axis:
imag_Bz = Bdot.Btang.mag.*sin(Bdot.Btang.phase);
% =========================================================================
% =========================================================================




% Plot data:
% =========================================================================
% Bdot probe data with respect to radial position and time:

% Real part of phasor:
figure('color','w'); 
surf(Bdot.Btang.RR,Bdot.Btang.TT,real_Bz,'LineStyle','none')
hL = line([-1,6],[1,1]*Bdot.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([-0.5,6])
zlabel('Re$\{B_z\}$ [mT]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('Real part: B-dot probe','interpreter','latex','Fontsize',12)

saveFig = 1;

if saveFig
    saveas(gcf,'Step_4_RealPartBdotData','tiffn')
end

% Iaginary part of phasor:
figure('color','w'); 
surf(Bdot.Btang.RR,Bdot.Btang.TT,imag_Bz,'LineStyle','none')
hL = line([-1,6],[1,1]*Bdot.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([-0.5,6])
zlabel('Im$\{B_z\}$ [mT]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('Imag part: B-dot probe','interpreter','latex','Fontsize',12)

saveFig = 1;

if saveFig
    saveas(gcf,'Step_4_ImagPartBdotData','tiffn')
end
