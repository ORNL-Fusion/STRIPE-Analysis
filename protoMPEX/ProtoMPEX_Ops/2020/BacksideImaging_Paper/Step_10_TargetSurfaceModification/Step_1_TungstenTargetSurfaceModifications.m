% Step 1: Create figure showing the surface modification on the front side
% of the target due to plasma exposure and deposition

close all
clear all

% Flags
% =========================================================================
saveFig = 1;

% Read image:
% =========================================================================
fileName = 'TungstenTargetSurfaceModification.jpg';
imdata = imread(fileName);

% Create figure:
% =========================================================================
figure('color','w');
image(imdata)
axis('image')
set(gca,'XTick',[])
set(gca,'YTick',[])
ylim([250,1350])

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_1_TargetSurfaceModification'];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end