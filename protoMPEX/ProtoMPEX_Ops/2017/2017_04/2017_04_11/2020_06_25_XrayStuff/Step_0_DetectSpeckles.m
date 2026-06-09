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
shotName_1 = 'slomo_1491928683.mov';
d1 = importdata(shotName_1);
f{1} = permute(d1,[1 2 4 3]);
clear d1

if 0
% Shot with DLP in:
shotName_2 = 'slomo_1491940126.mov';
d2 = importdata(shotName_2);
f{2} = permute(d2,[1 2 4 3]);
clear d2
end

shotlist = [13809,13840]; % DLP out and DLP in 
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
for jj = 1:2
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
for jj = 1:2
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
rfSignal = 100*rfSignal/max(rfSignal); % [kW]
rfSignal(find(t_rf{1}(1:end-1)>4.34)) = 0;
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
ylim(ax(2),[0,300])

nexttile
hold on
colororder({'k','k'})

% Gas pressure:
yyaxis right
ax(3) = gca;
plot(ax(3),t_pg6{1}(1:end-1),PG6{1}*2*0.1333,'m','lineWidth',2)
ylim(ax(3),[0,0.7])

% Magnetic current:
yyaxis left
ax(4) = gca;
plot(ax(4),t_ps1{1}(1:end-1),PS1{1},'c','lineWidth',2)
ylim(ax(4),[0,7])

% Full figure formatting:
xlim(ax,[4.1,4.5])
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