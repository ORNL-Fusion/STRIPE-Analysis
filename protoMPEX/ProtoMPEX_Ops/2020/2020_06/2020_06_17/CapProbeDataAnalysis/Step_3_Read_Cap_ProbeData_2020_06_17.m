% Demonstrate how to read data from structure:
clear all
close all

% Load data:
load('Cap_ProbeData_2020_06_17.mat')

% #########################################################################
% INSTRUCTIONS:
% #########################################################################
% Two variables are loaded: Cap
% These are structures that contain all the data relating to the
% experiment on 20_06_17 
% Please the see the contents of the structures: Cap

% #########################################################################
% Cap probe data:
% #########################################################################
% Display contents on .mat file in command screen:
Cap

% Display of data:
Cap.Vrms


saveFig = 1;

% Cap probe data with respect to radial position and time:
figure('color','w'); 
surf(Cap.Vrms.RR,Cap.Vrms.TT,Cap.Vrms.mag,'LineStyle','none')
hL = line([-1,13],[1,1]*Cap.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
ylim([4.16,4.25])
xlim([-0.5,13])
zlabel('$|\widetilde{V}_{RMS}|$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('r [cm]','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend(hL,'Start of RF');
hLeg.Location = 'best'
title('Capacitive probe','interpreter','latex','Fontsize',12)

if saveFig
    saveas(gcf,'Step_3_CapData','tiffn')
end

% Cap probe data with respect to flux coordinate and time:
figure('color','w'); 
surf(Cap.Vrms.Xi,Cap.Vrms.TT,Cap.Vrms.mag,'LineStyle','none')
hL  = line([-1,10],[1,1]*Cap.t_rfStart,[0.1,0.1],'Color','r','LineWidth',3)
hXi = line([1,1],[0,5],[0.1,0.1],'Color','r','LineWidth',2,'LineStyle','--')
ylim([4.16,4.25])
xlim([0,10])
zlabel('$|\widetilde{V}_{RMS}|$ [V]','interpreter','latex','Fontsize',12)
ylabel('time [sec]','interpreter','latex','Fontsize',12)
xlabel('$\sqrt{\chi}$','interpreter','latex','Fontsize',12)
set(gca,'PlotBoxAspectRatio',[1 1.5 1],'FontName','Times','FontSize',10)
view([30,60])
colorbar
hLeg = legend([hL,hXi],'Start of RF','$\chi$ = 1');
hLeg.Location = 'best';
hLeg.Interpreter = 'latex';
title('Capacitive probe','interpreter','latex','Fontsize',12)

if saveFig
    saveas(gcf,'Step_3_CapData_Xi','tiffn')
end
