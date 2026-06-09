% Script to check content of heat flux data structure:

clear all
close all
clc

% Load file and view comments:
% =========================================================================
load Step_4c_MPEX_ICH_HeatFlux2D_155kW_CoupledPwr_created_18-Dec-2020.mat

for jj = 1:numel(dataSet_ICH)
    % Print basic information:
    % =========================================================================
    dataSet_ICH{jj}.comment'
    disp(['Seal Type : ', num2str(dataSet_ICH{jj}.MPEX_sealType)]);
    disp(['Window length : ', num2str(dataSet_ICH{jj}.LengthWindow), ' [m]']);
    disp(['Window Inner radius : ', num2str(dataSet_ICH{jj}.RadiusWindow), ' [m]']);
    disp(['RF coupled power : ', num2str(dataSet_ICH{jj}.RFpower_Coupled), ' [kW]']);
    disp(['Antenna length : ', num2str(dataSet_ICH{jj}.AntennaLength), ' [m]']);

    % Get heat flux
    % =========================================================================
    % Azimuthal angle:
    phi = dataSet_ICH{jj}.phi_q2D;
    % Axial coordinate:
    z = dataSet_ICH{jj}.z_q2D;
    % Heat flux:
    heatflux = dataSet_ICH{jj}.q2D;

    % Plot for confirmation
    % =========================================================================
    figure('color','w')
    image(phi(1,:)*180/pi,z(:,1)*1e2,heatflux/1000,'CDataMapping','scaled')
    colorbar
    colormap(hot)
    caxis([0 1800]);
    set(gca,'XDir','reverse','YDir','normal')
    xlabel('$\theta$ [deg]','interpreter','Latex','FontSize',12)
    ylabel('z [cm]','interpreter','Latex','FontSize',12)
    titleText{jj} = {'MPEX ICH window heat flux [kWm$^{-2}$], 175 kW Coupled power, seal: ', num2str(dataSet_ICH{jj}.MPEX_sealType)};
    title(titleText{jj},'interpreter','Latex','FontSize',12)
    ylim([0,45])
    
    % Draw Left handed antenna:
    if jj == 1
        phi_offset = 0;
    elseif jj == 2
        phi_offset = 90;
    end
    lw = 13;
    % Transverse straps:
    L1 = dataSet_ICH{jj}.L1;
    L2 = dataSet_ICH{jj}.L2;
    Lw = dataSet_ICH{jj}.LengthWindow;
    Lcenter = Lw/2;
    line([360,000],(Lw-L1)*[1,1]*1e2,[500,500],'color','k','LineWidth',lw)
    line([360,000],L2*[1,1]*1e2,[500,500],'color','k','LineWidth',lw)
    % Bottom helical strap:
    line([135,225] + phi_offset,[L2,(Lw-L1)]*1e2,[500,500],'color','k','LineWidth',lw)
    % HV top side helical strap:
    line([000,045] + phi_offset,[Lcenter,(Lw-L1)]*1e2,[500,500],'color','k','LineWidth',lw)
    % GND top side helical strap:
    line([315,360] + phi_offset,[L2,Lcenter]*1e2,[500,500],'color','k','LineWidth',lw)
    drawnow
    
    % Annotation:
    text(050,17,'HV','interpreter','Latex','FontSize',14)
    text(340,22,'GND','interpreter','Latex','FontSize',14)

    % Save figure
    % =========================================================================
    figureName = ['Step_5b_CheckContents_Mpex_ICH_HeatFlux_155kWCoupled_',num2str(dataSet_ICH{jj}.MPEX_sealType)];
    saveas(gcf,figureName,'tiffn')
end