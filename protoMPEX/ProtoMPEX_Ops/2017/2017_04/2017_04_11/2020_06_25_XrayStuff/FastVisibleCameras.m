% Read videos:

clear all
close all

% Shot with Xray specles:
% d1 = importdata('slomo_1491928683.mov');
% f{1} = permute(d1,[1 2 4 3]);

d1 = importdata('shot_BW_14056.mov');
f{1} = permute(d1,[1 2 4 3]);

% Shot with DLP in:
d2 = importdata('slomo_1491940126.mov');
f{2} = permute(d2,[1 2 4 3]);
clear d1 d2

% Load RF pulse:
shotlist = [13809,13840]; % DLP out and DLP in 
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'PWR_28GHZ']; % Isx
[ECH,t_ech]   = my_mdsvalue_v2(shotlist,DA(1))
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shotlist,DA(1))
DA{1} = [RootAddress,'PS1_I'];
[PS1,t_ps1]   = my_mdsvalue_v2(shotlist,DA(1))
DA{1} = [RootAddress,'PS1_I'];
[PS1,t_ps1]   = my_mdsvalue_v2(shotlist,DA(1))
address{1} = '\MPEX::TOP.MACHOPS1:PG3'
[PG6,t_pg6]   = my_mdsvalue_v2(shotlist,address(1))

for jj = 1:2;
    for ii = 1:size(f{jj},3)
            inten{jj}(:,:,ii) = f{jj}(:,:,ii,2);
            rng_r = 1:10:size(inten{jj},1);
            if jj == 2
                rng_c = 600:1:size(inten{jj},2);
            else
                rng_c = 1:1:size(inten{jj},2);
            end
            sum_inten{jj}(ii) = sum(sum(inten{jj}(rng_r,rng_c,ii)))/(size(f{jj},1)*size(f{jj},2));
            max_inten{jj}(ii) = double(max(max(inten{jj}(rng_r,rng_c,ii))));
    end
end

dt = 1e-3;
t1 = ([1:1:size(f{1},3)]')*dt;

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

saveFig = 1;
figureName = 'CentralChamberXrays_2017_04_11';
if saveFig
    saveas(gcf,figureName,'tiffn')
end
