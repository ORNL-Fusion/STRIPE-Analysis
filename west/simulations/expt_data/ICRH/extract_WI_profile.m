clear all
close all
clc
shot = 57877
load(['ICRH_' num2str(shot) '_8s.mat'])
%data.VISTA.S57877.ROI: {'LODIVIN','LODIVOU','UPDIV','IC2L','IC2R'}
a = 5
%data.VISTA.S57877.impurity: {'WI_4009','DI_4341','OII_4415','BII_4122','NII_3995','CII_3921'}
b = 1
%plot visible spectro on the IC antenna
figure()
hold on
errorbar(data.VISTA.S57877.Spectro_ROI{1, a}.position{1,b},0.1.*data.VISTA.S57877.Spectro_ROI{1, a}.mean_profile{1,b},0.1.*data.VISTA.S57877.Spectro_ROI{1, a}.std_profile{1,b},'r^','linewidth',2)
grid on
l=legend('57877','Location','northeast')
ylabel({'Brightness (ions.m^{-2}.s^{-1})'},'FontName','arial','fontsize',16,'FontWeight','bold');
xlabel({'Vertical position along Q2 right limiter: z (cm)'},'FontName','arial','fontsize',16,'FontWeight','bold');
set(gca,'Fontsize',15,'Fontweight','bold','LineWidth',2,'GridLineStyle','--')
set(gcf,'color','w');