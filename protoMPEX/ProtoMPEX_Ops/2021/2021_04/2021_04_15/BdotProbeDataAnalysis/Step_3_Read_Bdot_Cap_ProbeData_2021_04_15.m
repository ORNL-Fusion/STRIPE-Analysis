% Demonstrate how to read data from structure:
clear all
close all

% Load data:
load('Bdot_Cap_ProbeData_2020_06_09.mat')

% #########################################################################
% INSTRUCTIONS:
% #########################################################################
% Two variables are loaded: Bdot and Cap
% These are structures that contain all the data relating to the
% experiment on 20_06_09 
% Please the see the contents of the structures Bdot and Cap


% #########################################################################
% Bdot probe data:
% #########################################################################

% Display contents on .mat file in command screen:
Bdot

% Display contents of metadata:
disp(Bdot.dateOfExperiment)
disp(Bdot.comment)

% Display of data:
Bdot.Bnorm

% Bdot probe data with respect to radial position and time:
figure('color','w'); 
surf(Bdot.Btang.RR,Bdot.Btang.TT,Bdot.Btang.mag,'LineStyle','none')
hL = line([-1,6],[1,1]*Bdot.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([-0.5,6])
zlabel('$|\widetilde{B}_z|$ [mT]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('B-dot probe','interpreter','latex','Fontsize',12)

saveFig = 1;

if saveFig
    saveas(gcf,'Step_3_BdotData','tiffn')
end

% Bdot probe data with respect to flux coordinate and time:
figure('color','w'); 
surf(Bdot.Btang.Xi,Bdot.Btang.TT,Bdot.Btang.mag,'LineStyle','none')
hL = line([-1,6],[1,1]*Bdot.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([0,1])
zlabel('$|\widetilde{B}_z|$ [mT]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('B-dot probe','interpreter','latex','Fontsize',12)

saveFig = 1;

if saveFig
    saveas(gcf,'Step_3_BdotData_Xi','tiffn')
end


% #########################################################################
% Cap probe data:
% #########################################################################
% Display contents on .mat file in command screen:
Cap

% Display of data:
Cap.Vrms

% Cap probe data with respect to radial position and time:
figure('color','w'); 
surf(Cap.Vrms.RR,Cap.Vrms.TT,Cap.Vrms.mag,'LineStyle','none')
hL = line([-1,6],[1,1]*Cap.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([-0.5,6])
zlabel('$|\widetilde{V}_{RMS}|$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('Capacitive probe','interpreter','latex','Fontsize',12)

saveFig = 1;

if saveFig
    saveas(gcf,'Step_3_CapData','tiffn')
end

% Cap probe data with respect to flux coordinate and time:
figure('color','w'); 
surf(Cap.Vrms.Xi,Cap.Vrms.TT,Cap.Vrms.mag,'LineStyle','none')
hL = line([-1,6],[1,1]*Cap.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([0,1])
zlabel('$|\widetilde{V}_{RMS}|$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('Capacitive probe','interpreter','latex','Fontsize',12)

saveFig = 1;

if saveFig
    saveas(gcf,'Step_3_CapData_Xi','tiffn')
end
