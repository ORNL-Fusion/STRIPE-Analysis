% IR thermographic temperature measurement on a spare AlN window using a
% heat gun and thermocouple for comparison
clear all 
close all

figureName{1} = 'Step_3_AlN_Emissivity_90deg';
figureName{2} = 'Step_3_AlN_Emissivity_12deg';
figureName{3} = 'Step_3_AlN_Emissivity_TwoAngles';

% 90 deg viewing angle:
% =========================================================================
tc_90deg = [27,32,33  ,35  ,36   ,38  ,39  ,41  ,43  ,46  ,48  ];
ir_90deg = [28,33,34.1,36.4,36.2 ,39.1,39.6,42.3,43.6,47.1,48.6];
em_90deg = 0.9;

fontSize1 = 14;
fontSize2 = 14; 

figure('color','w')
hold on
x = linspace(0,50);
set(gca,'FontName','times','Fontsize',fontSize1)
plot(x,x,'k')
plot(tc_90deg,ir_90deg,'rsq','MarkerSize',10)
ylim([0,50])
ylabel('IR [C]','Interpreter','latex','fontsize',fontSize2)
xlim([0,50])
xlabel('Thermocouple [C]','Interpreter','latex','fontsize',fontSize2)
box on
grid on
title('AlN window, $\epsilon$ = 0.9','interpreter','latex','fontsize',fontSize2)
set(gcf,'position',[241.0000  313.0000  370.0000  312.6667])

% Save figure:
% =========================================================================
saveas(gcf,figureName{1},'tiffn')

% 12 deg viewing angle:
% =========================================================================
tc_12deg = [60 59 56.4 55.7 55.2 54.1 52.5 50.3 47.6 46.6 44.7 43.0 41.9 40.8 39.6 36.0 33.6 32.7 32.1 31.0 29.9 29.6 28.3 27.0 26.1];
ir_12deg = [61 59 57.4 56.2 55.4 54.2 52.8 50.5 47.6 46.8 45.1 43.4 42.1 41.1 39.8 36.4 33.3 32.5 31.8 30.5 29.4 29.0 27.5 26.3 25.1];
em_12deg = 0.55;

figure('color','w')
hold on
x = linspace(0,70);
set(gca,'FontName','times','Fontsize',fontSize1)
plot(x,x,'k')
plot(tc_12deg,ir_12deg,'rsq','MarkerSize',10)
ylim([0,70])
ylabel('IR [C]','Interpreter','latex','fontsize',fontSize2)
xlim([0,70])
xlabel('Thermocouple [C]','Interpreter','latex','fontsize',fontSize2)
box on
grid on
title('AlN window, $\epsilon$ = 0.55','interpreter','latex','fontsize',fontSize2)
set(gcf,'position',[241.0000  313.0000  370.0000  312.6667])

% Save figure:
% =========================================================================
saveas(gcf,figureName{2},'tiffn')

% Both 90 and 12 deg viewing angle:
% =========================================================================
fontSizeAxes = 11;
fontSizeLeg = 12;

figure('color','w')
hold on
x = linspace(0,70);
set(gca,'FontName','times','Fontsize',fontSizeAxes)
plot(x,x,'k')
hem(2) = plot(tc_12deg,ir_12deg,'rsq','MarkerSize',10,'MarkerFaceColor','r');
legendText{2} = '$\epsilon$ = 0.55 at 12$^\circ$';
hem(1) = plot(tc_90deg,ir_90deg,'ko' ,'MarkerSize',10,'MarkerFaceColor','k')
legendText{1} = '$\epsilon$ = 0.9 at 90$^\circ$';
ylim([0,70])
ylabel('IR [C]','Interpreter','latex','fontsize',fontSize2)
xlim([0,70])
xlabel('Thermocouple [C]','Interpreter','latex','fontsize',fontSize2)
box on
grid on
title('AlN window','interpreter','latex','fontsize',fontSize2)
hLeg = legend(hem,legendText,'interpreter','latex','fontsize',fontSizeLeg);
hLeg.Location = 'southeast';
hLeg.FontSize = fontSizeLeg;
set(gcf,'position',[241.0000  313.0000  370.0000  312.6667])

% Save figure:
% =========================================================================
saveas(gcf,figureName{3},'tiffn')
