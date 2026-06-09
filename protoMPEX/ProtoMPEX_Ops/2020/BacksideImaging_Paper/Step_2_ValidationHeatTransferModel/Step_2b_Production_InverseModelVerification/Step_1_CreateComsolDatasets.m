% Step 1: Create COMSOL datasets based on text files
% This requires assembling them into .mat files

clear all
close all
clc

% Data extraction options:
% =========================================================================
saveData = 0;
extractData = 0;
saveFig = 1;
datasets = [1:8];
targetDir = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\BacksideImaging_Paper\Step_2_ValidationHeatTransferModel\Step_2b_Production_InverseModelVerification\fwdSol_comsol_data\';

for dd = datasets
    % Create filename:
    % =====================================================================
    fileName{1} = ['T2D_BackSurface_case_' ,num2str(datasets(dd)),'.txt'];
    fileName{2} = ['q2D_FrontSurface_case_',num2str(datasets(dd)),'.txt'];
    matFileName = ['fwdSol_comsol_case_'   ,num2str(datasets(dd)),'.mat']; 
    
    if extractData
        % Read data:
        % =================================================================
        d{1} = readcell([targetDir,fileName{1}]);
        d{2} = readcell([targetDir,fileName{2}]);

        % Spatial coordinates:
        % =================================================================
        xx = str2num(cell2mat(d{1}(2)));
        yy = str2num(cell2mat(d{1}(3)));
        NS = 150;
        tt = linspace(0,1.5,NS);

        % Extract data:
        % =================================================================
        kk = 5;
        for jj = 1:NS
            for ii = 1:numel(xx)
                T(ii,:,jj) = str2num(cell2mat(d{1}(kk)))';
                q(ii,:,jj) = str2num(cell2mat(d{2}(kk)))';
                kk = kk + 1;
            end
            kk = kk + 1;
        end
        
        % Save data:
        % =================================================================
        if saveData
            varList  = {'xx','yy','tt','T','q'};
            save(matfileName,varList{:});
        end
    else
        load(matfileName);
    end
    
    % Plot data:
    % =====================================================================
    if 0
        figure('color','w'); 
        for jj = 1:1:150
            mesh(xx,yy,T(:,:,jj));
            caxis([30,1000])
            zlim([30,1000])
            view([30,10])
            drawnow
        end
    end

    if 0
        figure('color','w'); 
        for jj = 1:1:150
            mesh(xx,yy,q(:,:,jj));
            caxis([0,5]*1E6)
            zlim([0,5]*1E6)
            view([30,10])
            drawnow
        end
    end

    figure('color','w'); 
    hfig(1) = gcf;
    plot(tt,permute(T(50,50,:),[3,1,2]),'LineWidth',2)
    ylim([300,1000])
    grid on
    box on
    title(['COMSOL case ',num2str(dd)])
    set(gca,'fontName','Times','fontSize',13)
    xlabel('t [ms]')
    ylabel('T [K]')   
    
    figure('color','w'); 
    hfig(2) = gcf;    
    plot(tt,permute(q(50,50,:)*1E-6,[3,1,2]),'LineWidth',2)
    ylim([0,30])
    grid on
    box on
    title(['COMSOL case ',num2str(dd)])
    set(gca,'fontName','Times','fontSize',13)
    ylabel('[MWm$^{-2}$]','interpreter','latex')
    xlabel('[s]','interpreter','latex')    
    
    if saveFig
        figureName{1} = ['Step_1_T_case_',num2str(dd)];
        figureName{2} = ['Step_1_q_case_',num2str(dd)];
        
        % PDF figure:
        exportgraphics(hfig(1),[figureName{1},'.pdf'],'Resolution',300) 
        exportgraphics(hfig(2),[figureName{2},'.pdf'],'Resolution',300) 
        
        % TIFF figure:
        exportgraphics(hfig(1),[figureName{1},'.tiff'],'Resolution',600) 
        exportgraphics(hfig(2),[figureName{2},'.tiff'],'Resolution',600) 
    end
    
    
end


