% Image of the stainless steel mirror:
clear all
close all

% Figure name:
% =========================================================================
figureName = 'Step_3_StainlessSteelMirror';

% Load raw figure:
% =========================================================================
f = imread('StainlessSteelMirror.jpg');

% Define cropped region:
% =========================================================================
dx1 = 600;
dx2 = 0;
x1 = size(f,1);
x2 = size(f,2);
rng1 = (1+dx1):(x1-dx1) - 150;
rng2 = (1+dx2):(x2-dx2);

% Plot cropped image:
% =========================================================================
figure('color','w')
hI = image(f(rng1,rng2,:));
axis image
set(gca,'XTickLabel',[],'YTickLabel',[])

% Arrows:
% Conflat flange:
fields.String = '2 3/4 Conflat flange';
x =      [645 ,700 ];
y =      [1000,900];
hta = myTextArrow(gca,x,y,fields);
set(hta,'interpreter','Latex','FontSize',15)

% Stainless-steel mirror:
fields.String = {'Stainless-steel';'mirror'};
x =      [392,500];
y =      [800,600 ];
hta = myTextArrow(gca,x,y,fields);
set(hta,'interpreter','Latex','FontSize',15)

% Save figure:
% =========================================================================
saveas(gcf,figureName,'tiffn')