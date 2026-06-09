% Extract data for Juergen
% we need to take the ratio between D Balmer lines for shots 6510, 6511 and
% 6514. The data for the Balmer lines gamma, beta and alpha where taken
% with the fast camera. the data for this camera is located in the:
% "Protompex_data\Visible_Cameras\2015_12_11c" directory

close all
clear all

% Gather RF power data from the MDSserver
mdsconnect('mpexserver')
Shots = [6510,6511,6514];
%Shots = [6506,6507,6508];
address{1} = '\MPEX::TOP.FSCOPE:TUBE08:PMT_VOLT';
address{2} = '\MPEX::TOP.MPEX1:SPARE'; % 28 GHz signal
address{3} = '\MPEX::TOP.MACHOPS1:TRANS_I'; % Helicon signal
address{4} = '\MPEX::TOP.T_ZERO'; % Tzero

[Tube,t_tube] = my_mdsvalue_v2(Shots,address(1));
[RF_28GHZ,t_28] = my_mdsvalue_v2(Shots,address(2));
[RF_helicon,t_rf] = my_mdsvalue_v2(Shots,address(3));
[T_zero,~] = my_mdsvalue_v2(Shots,address(4));

% Gather Visible light emission from the fast cameras
Address = 'Z:\Visible_Cameras\2015_12_11c';
home = cd;
cd(Address);
% The times and associated files that we need to import:
% 6510: Gamma : 1:31 pm : 6:31 pm = slomo_1449858670.mov : 1/3e4
% 6511: Beta  : 1:35 pm : 6:35 pm = slomo_1449858930.mov : 1/3e4
% 6514: Alpha : 1:52 pm : 6:52 pm = slomo_1449859948.mov : 1/6e4

% Note that the data is imported as unsigned 8 bit intergers
V{1} = importdata('slomo_1449858670.mov'); % Gamma
V{2} = importdata('slomo_1449858930.mov'); % Beta
V{3} = importdata('slomo_1449859948.mov'); % Alpha

% Alpha: 6506: 1:16 pm : 6:16 pm = slomo_1449857754.mov : 1/4e4
% Beta : 6507: 1:20 pm : 6:20 pm = slomo_1449858029.mov : 1/2e4 
% Gamma: 6508: 1:23 pm : 6:23 pm = slomo_1449858195.mov : 1/2e4

%V{1} = importdata('slomo_1449858195.mov'); % Gamma
%V{2} = importdata('slomo_1449858029.mov'); % Beta
%V{3} = importdata('slomo_1449857754.mov'); % Alpha

cd(home)
%%
CentralPix_x = 325;
dPix_x = 20; Pix_x = [CentralPix_x-dPix_x : CentralPix_x + dPix_x];
%CentralPix_y = 680; % Saturates alpha at this position
CentralPix_y = 650;
dPix_y = 20; Pix_y = [CentralPix_y-dPix_y : CentralPix_y + dPix_y];

for s = 1:400;
    Gamma(:,:,s) = 0.5*double(V{1}(:,:,3,s)); 
%     GammaTarget(s) = mean(mean(Gamma(Pix_x,Pix_y,s)));
    GammaTarget(s) = Gamma(CentralPix_x,CentralPix_y,s);
    
    Beta(:,:,s) = 0.5*double(V{2}(:,:,3,s)); 
    Beta2(:,:,s) = 0.5*double(V{2}(:,:,2,s)); 
    BetaTarget(s) = Beta(CentralPix_x,CentralPix_y,s) + Beta2(CentralPix_x,CentralPix_y,s);    
%     BetaTarget(s) = mean(mean(Beta(Pix_x,Pix_y,s) + Beta2(Pix_x,Pix_y,s)));
    
    Alpha(:,:,s) = double(V{3}(:,:,1,s)); 
    AlphaTarget(s) = Alpha(CentralPix_x,CentralPix_y,s);
%     AlphaTarget(s) = mean(mean(Alpha(Pix_x,Pix_y,s)));
%     AlphaTarget2(s)= mean(mean(Alpha(Pix_x + 150,Pix_y,s)));
end

%% 
close all

% movie diraction 0.4 seconds and 400 frames; hence we have dt = 0.4/400
T = 0.4;
dt =  T/400;
tStart = 4.0638;
tVid = tStart:dt:(tStart + T - dt);

figure;
n = find(tVid>=4.295); n = n(1);
surf(Alpha(:,:,n),'LineStyle','none')
hold on
plot3(CentralPix_y,CentralPix_x,Alpha(CentralPix_x,CentralPix_y,n),'ko')

tStartPlot = 4.1;
tEndPlot = 4.35;

% Smoothing data:

figure; 
subplot(3,1,1);hold on; 
h(1) = plot(tVid,GammaTarget,'bl');
h(2) = plot(tVid,BetaTarget,'g'); 
h(3) = plot(tVid,AlphaTarget,'r');
L = legend(h,'$D_\gamma$','$D_\beta$','$D_\alpha$','location','NorthWest'); set(L,'interpreter','Latex','FontSize',12,'box','off')
xlim([tStartPlot,tEndPlot]); ylim([0,300])
ax(1) = gca; Pos{1} = get(ax(1),'position');
set(ax(1),'Position',[0.1300    0.11 + 2*0.2614    0.7750    0.2614],'XTick',[])
ylabel('$$ D_{\alpha\beta\gamma} $$ [arb.]','Interpreter','Latex','FontSize',15)

subplot(3,1,2); hold on
h(1) = plot(tVid,GammaTarget./BetaTarget,'k','LineWidth',2)
h(2) = plot(tVid,GammaTarget./AlphaTarget,'k','LineWidth',1)
L = legend(h,'\gamma/\beta','\gamma/\alpha','location','NorthWest'); set(L,'FontSize',12,'box','off')
ylabel('$$ Ratios $$','Interpreter','Latex','FontSize',13)
set(gca,'YTick',[0,0.5,1,1.5])
xlim([tStartPlot,tEndPlot]); ylim([0,2])
ax(2) = gca; Pos{2} = get(ax(2),'position');
set(ax(2),'Position',[0.1300    0.1100 + 0.2614    0.7750    0.2614],'XTick',[])

subplot(3,1,3); hold on
ECH_factor = 10.6; % 10.6 kW/v
h(1) = plot(t_rf{1},(104/0.67)*RF_helicon{1},'k','LineWidth',2)
h(2) = plot(t_28{1},ECH_factor*sgolay_t(RF_28GHZ{1},3,101),'k')
L = legend(h,'Helicon','ECH','location','NorthWest'); set(L,'FontSize',10,'box','off')
ylabel('$$ Power $$ $$ [kW] $$','Interpreter','Latex','FontSize',13)
set(gca,'YTick',[0,50,100])
xlim([tStartPlot,tEndPlot]); ylim([0,150])
ax(3) = gca; Pos{3} = get(ax(3),'position');
set(ax(3),'Position',[0.1300    0.1100    0.7750    0.2614])
xlabel('$Time$ $[s]$','Interpreter','Latex','FontSize',12) 

set(gcf,'color','w')
set(ax,'box','on')
