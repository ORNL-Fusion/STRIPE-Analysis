% Step 1: read emissivity data as a function of temperature of the surface
% treated W target

close all
clear all

% Flags
% =========================================================================
saveFig = 1;

% Read excell spreadsheet:
% =========================================================================
fileName = '2019_10_29_Emissivity_Wplate.xlsx';
T = readtable(fileName);

% Assign data:
% =========================================================================
% x = y line:
% ===========
TC{1} = T.TC + 273.15; %[K]

% 0.35 emissivity:
% ================
TC{2} = T.TC_1 + 273.15;
IR{2} = T.EmissivityAt0_35 + 273.15;

% 0.45 IR:
% ================
TC{3} = T.TC_2 + 273.15;
IR{3} = T.EmissivityAt0_45 + 273.15;

% 0.3a IR:
% ================
TC{4} = T.TC_3 + 273.15; 
IR{4} = T.EmissivityAt0_3 + 273.15;

% 0.3b IR:
% ================
TC{5} = T.TC_4 + 273.15;
IR{5} = T.EmissivityAt0_3_1 + 273.15;

% Create figure:
% =========================================================================
% All data points:
% ================
figure('color','w')
hold on
Tmax = 160 + 273.15;
Tmin = 273;
h(1) = line([0,1]*Tmax,[0,1]*Tmax);
h(2) = plot(TC{2},IR{2});
h(3) = plot(TC{3},IR{3});
h(4) = plot(TC{4},IR{4});
h(5) = plot(TC{5},IR{5});

% Formatting:
xlim([Tmin,Tmax])
ylim([Tmin,Tmax])
set(gca,'FontSize',11,'FontName','Times')
box on
set(h(1),'color','k','LineWidth',1)
set(h(2),'color','r','Marker','o'      ,'LineStyle','-')
set(h(3),'color','k','Marker','sq'     ,'LineStyle','-')
set(h(4),'color','g','Marker','*'      ,'LineStyle','-')
set(h(5),'color','m','Marker','diamond','LineStyle','-')

% Legend:
hL = legend(h(2:end),'0.35','0.45','0.3','0.3','');

% Labels:
xlabel('TC [K]','Interpreter','Latex','FontSize',12)
ylabel('IR [K]','Interpreter','Latex','FontSize',12)

% Final figure:
% =============
figure('color','w')
hold on
Tmax = 650 + 273.15;
Tmin = 0;
hF(1) = line([0,1]*Tmax,[0,1]*Tmax);
hF(2) = plot(TC{4},IR{4});
hF(3) = plot(TC{5},IR{5});

% Formatting:
xlim([Tmin,Tmax])
ylim([Tmin,Tmax])
set(gca,'FontSize',11,'FontName','Times')
box on
set(hF(1),'color','k','LineWidth',1)
set(hF(2),'color','k','Marker','o','LineStyle','none','MarkerFaceColor','k')
set(hF(3),'color','k','Marker','o','LineStyle','none','MarkerFaceColor','k')
set(gcf,'Position',[360   314   388   304]);

% Legend:
hL = legend(hF(2),'$\epsilon$ = 0.3');
hL.Location = 'northwest';
hL.Interpreter = 'Latex';
hL.FontSize = 14;

% Labels:
xlabel('TC [K]','Interpreter','Latex','FontSize',12)
ylabel('IR [K]','Interpreter','Latex','FontSize',12)
title('W plate, 90 degree (normal) view','Interpreter','Latex','FontSize',12)

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_1_IRemissivityWplate'];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end