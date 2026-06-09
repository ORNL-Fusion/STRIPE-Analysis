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
    heatflux{jj} = dataSet_ICH{jj}.q2D;

    % Plot for confirmation
    % =========================================================================
    figure('color','w')
    image(phi(1,:)*180/pi,z(:,1)*1e2,heatflux{jj}/1000,'CDataMapping','scaled')
    colorbar
    colormap(hot)
    caxis([0 1800]);
    set(gca,'XDir','reverse','YDir','normal')
    xlabel('$\theta$ [deg]','interpreter','Latex','FontSize',12)
    ylabel('z [cm]','interpreter','Latex','FontSize',12)
    titleText{jj} = {'MPEX ICH window heat flux [kWm$^{-2}$], 155 kW Coupled power, seal: ', num2str(dataSet_ICH{jj}.MPEX_sealType)};
    title(titleText{jj},'interpreter','Latex','FontSize',12)
    ylim([0,45])
    
    % Draw Left handed antenna:
    % ======================================================================
    if jj == 1
        phi_offset = 0;
    elseif jj == 2
        phi_offset = 90;
    end
    lw = 13;
    L1 = dataSet_ICH{jj}.L1;
    L2 = dataSet_ICH{jj}.L2;
    Lw = dataSet_ICH{jj}.LengthWindow;
    Lcenter = Lw/2;
    
    % Transverse straps:
    % Target facing:
    line([360,000],(Lw-L1)*[1,1]*1e2,[500,500],'color','k','LineWidth',lw)
    % Dump facing:
    line([360,000],L2*[1,1]*1e2,[500,500],'color','k','LineWidth',lw)
    
    % Bottom helical strap:
    line([135,225] + phi_offset,[L2,(Lw-L1)]*1e2,[500,500],'color','k','LineWidth',lw)
    
    % HV top side helical strap:
    line([000,045] + phi_offset,[Lcenter,(Lw-L1)]*1e2,[500,500],'color','k','LineWidth',lw)
    % GND top side helical strap:
    line([315,360] + phi_offset,[L2,Lcenter]*1e2,[500,500],'color','k','LineWidth',lw)
    drawnow
    % ======================================================================

    % Annotation:
    text(050,17,'HV','interpreter','Latex','FontSize',14)
    text(340,22,'GND','interpreter','Latex','FontSize',14)

    % Save figure
    % =========================================================================
    figureName = ['Step_5b_CheckContents_MPEX_ICH_HeatFlux_155kWCoupled_',num2str(dataSet_ICH{jj}.MPEX_sealType)];
    saveas(gcf,figureName,'tiffn')
end

%% Create a table for DES seal type:

% Select data set:
jj = 1;

% The radius of the window is 0.0635 m

% Data to be put into table format:
% dataSet_ICH{1}.q2D
% dataSet_ICH{1}.phi_q2D
% dataSet_ICH{1}.z_q2D

% q2D is resolved in:
% - phi with 100 elements
% - z with 60 elements

% Save data to Excel spreadsheet:
% =========================================================================
% q2D:
% -------------------------------------------------------------------------
FileName = 'DES_HeatFluxMap.xlsx';
xlswrite(FileName,dataSet_ICH{jj}.q2D);

% z_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_z.xlsx';
xlswrite(FileName,num2cell(dataSet_ICH{jj}.z_q2D));

% phi_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_phi.xlsx';
xlswrite(FileName,num2cell(dataSet_ICH{jj}.phi_q2D));


% Save data to text file:
% ==========================================================================
% q2D:
% -------------------------------------------------------------------------
FileName = 'DES_HeatFluxMap.txt';
dlmwrite(FileName,dataSet_ICH{jj}.q2D);

% z_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_z.txt';
dlmwrite(FileName,num2cell(dataSet_ICH{jj}.z_q2D));

% phi_q2D:
% -------------------------------------------------------------------------
FileName = 'DES_phi.txt';
dlmwrite(FileName,num2cell(dataSet_ICH{jj}.phi_q2D));

%% Transform data into ANSYS CFX file format:

% Loop over the two seal designs:
for ii = 1:2
    % Determine the size of the heat flux 2D map:
    [nRows,nCols] = size(dataSet_ICH{ii}.q2D);
    
    % Print information about data to CLI:
    textMessage1 = ['Assembling data for ',dataSet_ICH{ii}.MPEX_sealType,' seal type ...'];
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
        HeatLoad{ii}(rng,1) = dataSet_ICH{ii}.RadiusWindow;
        % Azimuthal angle "phi":
        HeatLoad{ii}(rng,2) = dataSet_ICH{ii}.phi_q2D(:,cc); 
        % z coordinate:
        HeatLoad{ii}(rng,3) = dataSet_ICH{ii}.z_q2D(:,cc);        
        % Heat flux:
        HeatLoad{ii}(rng,4) = dataSet_ICH{ii}.q2D(:,cc);

    end
    
    figure('color','w')
    plot3(HeatLoad{ii}(:,2),HeatLoad{ii}(:,3),HeatLoad{ii}(:,4)*1e-6,'k.');
    xlabel('phi [Rad]','interpreter','latex');
    ylabel('z [m]','interpreter','latex')
    zlabel('heat flux [MWm$^{-3}$]','interpreter','latex')
    title(textMessage1,'interpreter','latex')
    
    ylim([0,0.5])
    view([60,60])
end

csvwrite("MPEX_ICH_155kW_CoupledRfPwr_DES_Inlet.csv",HeatLoad{1})
csvwrite("MPEX_ICH_155kW_CoupledRfPwr_DES_Outlet.csv",HeatLoad{2})
