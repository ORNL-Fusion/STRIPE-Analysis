% Step_0c: isometric view of electron heating section:

clear all
close all

saveFig = 1;

% Load raw figure:
% =========================================================================
f = imread('ElectronHeatingChamber.jpg');

% Plot cropped image:
% =========================================================================
figure('color','w')
hI = image(f);
axis image
set(gca,'XTickLabel',[],'YTickLabel',[])

% Arrows:
% =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 11;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 5;
fields.HeadStyle = 'vback2';

% Ellipsoidal reflector:
fields.String = 'G';
x =      [200,350];
y =      [600,700]- 450;
hta = myTextArrow(gca,x,y,fields);

% 28 GHz waveguide::
fields.String = 'D';
x =      [200,500];
y =      [600,750]- 550;
hta = myTextArrow(gca,x,y,fields);


% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0c_IsometricEchChamber';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end