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
    heatflux{jj} = dataSet{jj}.q2D;

    % Plot for confirmation
    % =========================================================================
    figure('color','w')
    image(phi(1,:)*180/pi,z(:,1)*1e2,heatflux{jj}/1000,'CDataMapping','scaled')
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

%% Create a table for DES seal type:

% Select data set:
jj = 1;

% The radius of the window is 0.0635 m

% Data to be put into table format:
% dataSet{1}.q2D
% dataSet{1}.phi_q2D
% dataSet{1}.z_q2D

% q2D is resolved in:
% - phi with 100 elements
% - z with 60 elements

% Save data to Excel spreadsheet:
% =========================================================================
% q2D:
% -------------------------------------------------------------------------
FileName = 'DES_HeatFluxMap.xlsx';
xlswrite(FileName,dataSet{jj}.q2D);

% z_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_z.xlsx';
xlswrite(FileName,num2cell(dataSet{jj}.z_q2D));

% phi_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_phi.xlsx';
xlswrite(FileName,num2cell(dataSet{jj}.phi_q2D));


% Save data to text file:
% ==========================================================================
% q2D:
% -------------------------------------------------------------------------
FileName = 'DES_HeatFluxMap.txt';
dlmwrite(FileName,dataSet{jj}.q2D);

% z_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_z.txt';
dlmwrite(FileName,num2cell(dataSet{jj}.z_q2D));

% phi_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_phi.txt';
dlmwrite(FileName,num2cell(dataSet{jj}.phi_q2D));

%% Transform data into ANSYS CFX file format:

% Loop over the two seal designs:
for ii = 1:2
    % Determine the size of the heat flux 2D map:
    [nRows,nCols] = size(dataSet{ii}.q2D);
    
    % Print information about data to CLI:
    textMessage1 = ['Assembling data for ',dataSet{ii}.MPEX_sealType,' seal type ...'];
    disp(textMessage1);
    
    textMessage2 = ['No of rows: ',num2str(nRows)];
    disp(textMessage2);
    
    textMessage3 = ['No of columns: ',num2str(nCols)];
    disp(textMessage3);
    
    % Initialize variable to hold data:
    HeatLoad{ii} = zeros(nRows*nCols,4);
    
    % Assemble data in CFX file format:
    for cc = 1:(nCols)
        rng = (nRows*(cc-1) + 1):nRows*cc;
        
        disp(num2str(rng(end)));
        
        % Window radius:
        HeatLoad{ii}(rng,1) = dataSet{ii}.RadiusWindow;
        % Azimuthal angle "phi":
        HeatLoad{ii}(rng,2) = dataSet{ii}.phi_q2D(:,cc); 
        % z coordinate:
        HeatLoad{ii}(rng,3) = dataSet{ii}.z_q2D(:,cc);        
        % Heat flux:
        HeatLoad{ii}(rng,4) = dataSet{ii}.q2D(:,cc);

    end
    
    figure('color','w')
    plot3(HeatLoad{ii}(:,2),HeatLoad{ii}(:,3),HeatLoad{ii}(:,4),'k.');
    xlabel('phi [Rad]','interpreter','latex');
    ylabel('z [m]','interpreter','latex')
    zlabel('heat flux [MWm$^{-3}$]','interpreter','latex')
    title(textMessage1,'interpreter','latex')
    
    ylim([0,0.5])
    view([60,60])
end

csvwrite("ProtoMPEX_CFX_200kW_fix_DES.csv",HeatLoad{1});
csvwrite("ProtoMPEX_CFX_200kW_fix_MCS.csv",HeatLoad{2});
