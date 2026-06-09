% Step_0d: IR Target view

clear all
close all

saveFig = 1;

% Load raw figure:
% =========================================================================
f = imread('IR_TargetView.png');

% Plot cropped image:
% =========================================================================
figure('color','w')
hold on
rng1 = 50:500; 
rng2 = 125:625;
hI = image(f(rng1,rng2,:));
axis image
set(gca,'XTickLabel',[],'YTickLabel',[],'YDir','reverse')

x0 = 256;
y0 = 235;
plot(x0,y0,'k.','MarkerSize',15)

% % Arrows:
% % =========================================================================
% Text arrows fields:
fields.HorizontalAlignment = 'center';
fields.Color = 'k';
fields.FontSize = 13;
fields.Interpreter = 'Latex';
fields.HeadLength = 5;
fields.HeadWidth = 7;
fields.HeadStyle = 'vback2';
fields.LineWidth = 2;

% Target alligned coordinate system:
theta = 44;
fields.String = [];
x =      [x0,x0 + 40*cosd(theta)] - 2;
y =      [y0,y0 + 40*sind(theta)] - 18;
hta = myTextArrow(gca,x,y,fields);

theta = 35;
fields.String = [];
x =      [x0,x0 + 37*cosd(theta-90)] - 2;
y =      [y0,y0 + 37*sind(theta-90)] - 18;
hta = myTextArrow(gca,x,y,fields);

hT(1)= text(275, 193,'$\hat{y}_*$');
hT(2)= text(285, 265,'$\hat{x}_*$');
set(hT,'interpreter','Latex','color','k','FontSize',15)

% Target label:
fields.String = ['W Target'];
x =      [256,256];
y =      [75,161] ;
fields.Color = 'r';
fields.FontSize = 15;
hta = myTextArrow(gca,x,y,fields);

% Target plate boundary:
% =========================================================================


% Save figure:
% =========================================================================
if saveFig
    figureName = 'Step_0d_IR_TargetView';

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',600,'ContentType', 'vector') 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end