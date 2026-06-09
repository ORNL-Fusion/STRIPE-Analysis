% Read visible camera videos and count speckles due to Xrays:
% =========================================================================

clear all
close all
clc

% Import data:
% =========================================================================

% Visible camera:
% -------------------------------------------------------------------------
% Shot with Xray specles:
shotName_1 = 'shot_BW_14056.mov';
shotlist = [14056];

% % t = 4.35s
% shotName_1 = 'slomo_14001.mov';
% shotlist = [14001];
% 
% % t = 4.30s
% shotName_1 = 'slomo_13996.mov';
% shotlist = [13996];

% % t = 4.20s
% shotName_1 = 'slomo_13991.mov';
% shotlist = [13991];

d1 = importdata(shotName_1);
f{1} = permute(d1,[1 2 4 3]);


Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

% 28 GHz pulse:
% -------------------------------------------------------------------------
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1));

% 13.56 MHz pulse:
% -------------------------------------------------------------------------
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1));

% DC magnet currents:
% -------------------------------------------------------------------------
DA{1} = [RootAddress,'PS1_I'];
[PS1,t_ps1]   = my_mdsvalue_v2(shotlist,DA(1));
DA{1} = [RootAddress,'PS1_I'];
[PS1,t_ps1]   = my_mdsvalue_v2(shotlist,DA(1));

% Electron heating neutral gas pressure:
% -------------------------------------------------------------------------
address{1} = '\MPEX::TOP.MACHOPS1:PG3'
[PG6,t_pg6]   = my_mdsvalue_v2(shotlist,address(1));

%% Post process visible light camera videos:
% =========================================================================
close all

% Extract intensity, mean intensity and max intensity profiles:
for jj = 1:numel(f)
    for ii = 1:size(f{jj},3)
       
        % Light intensity:        
        inten{jj}(:,:,ii) = double(f{jj}(:,:,ii,1));
                       
        % Select data range:
            rng_r = 1:20:size(inten{jj},1);
            if jj == 2
                rng_c = 600:10:size(inten{jj},2);
            else
                rng_c = 1:10:size(inten{jj},2);
            end
            
         % Mean intensity:
         mean_inten{jj}(ii) = mean(mean(inten{jj}(rng_r,rng_c,ii)));
         
         % Max intensity:
         max_inten{jj}(ii)  = max(  max(inten{jj}(rng_r,rng_c,ii)));
    end
end

% Time coordinate:
dt = 1e-3;
t1 = ([1:1:size(inten{1},3)]')*dt;

% Clear up memory:
clear f

%% Extract speckles:
% =========================================================================
for jj = 1:numel(max_inten)
    for ii = 1:size(inten{jj},3)
        % Calculate gradient of intensity:
        [Fx,Fy] = gradient(inten{jj}(:,:,ii));
        
        % Speckle gradient:
        speckleGradient{jj}(:,:,ii) = Fx;
        
        % Reject area where probe blows up:
        rngx = 1:300;
        rngy = 300:400;
        speckleGradient{jj}(rngx,rngy,ii) = 0;
        
        % Extract speckles:
        [nx{jj}{ii},ny{jj}{ii},vals{jj}{ii}] = find(speckleGradient{jj}(:,:,ii) >= 15);
        
        % Count speckles:
        speckleCount{jj}(ii) = numel(nx{jj}{ii});
    end
end

%% Plot data:
close all 

% Speckle count vrs time:
% =========================================================================
figure('color','w')
hTile = tiledlayout(2,1);

% Font size:
fontSize.axes = 11;
fontSize.label = 12;
fontSize.legend = 12;
fontSize.title = 12;

nexttile
hold on
colororder({'k','k'})

% 13.56 MHz power:
yyaxis right
ax(1) = gca;
rfSignal = abs(RF{1});
Nrf = numel(RF{1});
rng = round(0.5*Nrf):Nrf;
rfSignal = 100*rfSignal/max(rfSignal(rng)); % [kW]
rfSignal(find(t_rf{1}(1:end-1)>4.444)) = 0;
plot(ax(1),t_rf{1}(1:end-1),rfSignal,'g','lineWidth',2)

% 28 GHz power: 
yyaxis right
echSignal = abs(movmean(ECH{1},7));
echSignal = 20*echSignal/max(echSignal);
plot(ax(1),t_ech{1}(1:end-1),echSignal,'r-','lineWidth',2)
ylim(ax(1),[0,120])

% Speckle:
yyaxis left
ax(2) = gca;
t0 = 4 + 3e-3;
plot(ax(2),t1 + t0,speckleCount{1},'k','lineWidth',2)
ylim(ax(2),[0,800])

nexttile
hold on
colororder({'k','k'})

% Gas pressure:
yyaxis right
ax(3) = gca;
neutralGasPressure = PG6{1}*2*0.1333;
rng = find(t_pg6{1}(1:end-1) > 3,1);
plot(ax(3),t_pg6{1}(1:end-1),neutralGasPressure - neutralGasPressure(rng),'m','lineWidth',2)
ylim(ax(3),[0,0.7])

% Magnetic current:
yyaxis left
ax(4) = gca;
plot(ax(4),t_ps1{1}(1:end-1),PS1{1},'c','lineWidth',2)
ylim(ax(4),[0,7])

% Full figure formatting:
xlim(ax,[4.1,4.7])
xlabel(ax,'time [s]','interpreter','latex','fontSize',fontSize.label)
box(ax,'on')
set(ax,'fontName','Times','fontSize',fontSize.axes)
for ii = 1:numel(ax)
    ax(ii).YColor = 'k';
end


% Save figure:
saveFig = 1;
figureName = ['CentralChamberXrays_2017_04_11_',shotName_1(1:end-4)];
if saveFig
    saveas(gcf,figureName,'tiffn')
end

return
%%

figure('color','w')
hold on
plot(t1,mean_inten{1},'r','LineWidth',2)
plot(t1,max_inten{1},'r','LineWidth',1)
% plot(t1,mean_inten{2},'g','LineWidth',2)
% plot(t1,max_inten{2},'g','LineWidth',1)

figure('color','w')
hold on
plotType = 2;
switch plotType
    case 1     
        plot(t1,mean_inten{1},'r','LineWidth',2)
%         plot(t1,mean_inten{2},'g','LineWidth',2)
        plot(t1,speckleCount{1},'r')
%         plot(t1,speckleCount{2},'g')
    case 2   
        plot(mean_inten{1},'r','LineWidth',2)
        plot(mean_inten{2},'g','LineWidth',2)
        plot(speckleCount{1},'r')
        plot(speckleCount{2},'g')
end



return


for ii = 1:size(inten{jj},3)
    surf(inten{jj}(:,:,ii),'lineStyle','none')
    view([0,90])
%     caxis([])
    pause(0.01)
    drawnow
end



dt = 1e-3;
t1 = ([1:1:size(f{1},3)]')*dt;

return

plotTitle{1} = 'ECH, w/o skimmer probe'
plotTitle{2} = 'ECH, w/ skimmer probe'

figure('color','w')
for jj = 1:2;
    subplot(2,1,jj)
    hold on
    hinten(jj) = plot(t1+4.003,max_inten{jj}/75,'lineWidth',3)
    % plot(t_rf{1}(1:end-1),-RF{jj})
    hech(jj) = plot(t_ech{1}(1:end-1),ECH{jj},'lineWidth',3)
%     plot(t1+4.003,sum_inten{jj},'lineWidth',3)
    hps1(jj) = plot(t_ps1{1}(1:end-1),PS1{jj}/4,'lineWidth',3)
    hpg(jj) = plot(t_pg6{1}(1:end-1),PG6{jj},'lineWidth',3)
    legend([hinten(jj),hech(jj),hps1(jj),hpg(jj)],'fVC','28 GHz','PS1','PG 6.5')
    grid on
    xlim([4.1,4.6])
    ylim([0,4])
    hold off
    box on
end


