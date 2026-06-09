% Script to check content of heat flux data structure:

clear all
close all
clc

% Load file and view comments:
% =========================================================================
load Step_4c_MpexHeatFlux2D_175kW_CoupledPwr_created_13-Nov-2020.mat

for jj = 1:numel(dataSet)
    % Print basic information:
    % =========================================================================
    dataSet{jj}.comment'
    disp(['Seal Type : ', num2str(dataSet{jj}.MPEX_sealType)]);
    disp(['Window length : ', num2str(dataSet{jj}.LengthWindow), ' [m]']);
    disp(['Window Inner radius : ', num2str(dataSet{jj}.RadiusWindow), ' [m]']);
    disp(['RF coupled power : ', num2str(dataSet{jj}.RFpower_Coupled), ' [kW]']);
    disp(['RF net power : ', num2str(dataSet{jj}.RFpower_NET), ' [kW]']);
    disp(['Antenna length : ', num2str(dataSet{jj}.AntennaLength), ' [m]']);

    % Get heat flux
    % =========================================================================
    % Azimuthal angle:
    phi = dataSet{jj}.phi_q2D;
    % Axial coordinate:
    z = dataSet{jj}.z_q2D;
    % Heat flux:
    heatflux = dataSet{jj}.q2D;

    % Plot for confirmation
    % =========================================================================
    figure('color','w')
    image(phi(1,:)*180/pi,z(:,1)*1e2,heatflux/1000,'CDataMapping','scaled')
    colorbar
    colormap(hot)
    caxis([0 1000]);
    set(gca,'XDir','reverse','YDir','normal')
    xlabel('$\theta$ [deg]','interpreter','Latex','FontSize',12)
    ylabel('z [cm]','interpreter','Latex','FontSize',12)
    titleText{jj} = {'MPEX heat flux [kWm$^{-2}$], 175 kW Coupled power, seal: ', num2str(dataSet{jj}.MPEX_sealType)};
    title(titleText{jj},'interpreter','Latex','FontSize',12)
    ylim([0,45])
    
    % Draw antenna:
    lw = 13;
    % Transverse straps:
    L1 = dataSet{jj}.L1;
    L2 = dataSet{jj}.L2;
    Lw = dataSet{jj}.LengthWindow;
    Lcenter = Lw/2;
    line([360,000],(Lw-L1)*[1,1]*1e2,[500,500],'color','k','LineWidth',lw)
    line([360,000],L2*[1,1]*1e2,[500,500],'color','k','LineWidth',lw)
    % Bottom helical strap:
    line([225,135],[L2,(Lw-L1)]*1e2,[500,500],'color','k','LineWidth',lw)
    % HV top side helical strap:
    line([045,000],[L2,Lcenter]*1e2,[500,500],'color','k','LineWidth',lw)
    % GND top side helical strap:
    line([360,315],[Lcenter,(Lw-L1)]*1e2,[500,500],'color','k','LineWidth',lw)
    drawnow
    
    % Annotation:
    text(050,17,'HV','interpreter','Latex','FontSize',14)
    text(340,22,'GND','interpreter','Latex','FontSize',14)

    % Save figure
    % =========================================================================
    figureName = ['Step_5b_CheckContents_MpexHeatFlux_175kWCoupled_',num2str(dataSet{jj}.MPEX_sealType)];
    saveas(gcf,figureName,'tiffn')
end