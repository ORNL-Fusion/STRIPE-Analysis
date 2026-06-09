% 2017_11_09, Bdot probe scan at 4kA

close all
clear all
LoadData = 1; 

if LoadData
shot = 16300 + [];

Data{1} = ['\MPEX::TOP.MACHOPS1:MN_CURRENT']; % Vp0
Data{2} = ['\MPEX::TOP.MACHOPS1:EA_CURRENT']; % Vp90
Data{3} = ['\MPEX::TOP.MACHOPS1:GEN_RF_PWR']; % Vmag

% Data{1} = ['\MPEX::TOP.MACHOPS1:CENT_LP']; % Vp0
% Data{2} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1']; % Vp90
% Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2']; % Vmag

Data{4} = ['\MPEX::TOP.MACHOPS1:CENT_LP'];

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
DA{1} = [RootAddress,'RF_FWD_PWR'];
[RF,t_rf]   = my_mdsvalue_v2(shot,DA(1))
DA{1} = [RootAddress,'TARGET_LP'];
[Isat,t_isat]   = my_mdsvalue_v2(shot,DA(1))

Title{1} = 'Vp0';
Title{2} = 'Vp90';
Title{3} = 'Vmag0';
Title{4} = 'None';

figure;
for ch = 1:4 % For all channels on digitizer
    [f{ch},t{ch}] = my_mdsvalue_v2(shot,Data(ch));
    % f{channel}{shot}
    subplot(2,2,ch);
    hold on
    for s = 1:length(shot)
            rng = find(t{ch}{s}>= 4.12 & t{ch}{s}<= 4.5);
            t{ch}{s} = t{ch}{s}(rng);
            if isempty(rng)
                rng = find(t{ch}{s+1}>= 4.12 & t{ch}{s+1}<= 4.5);
                t{ch}{s} = t{ch}{s+1}(rng);
            end
            f{ch}{s} = sgolay_t(f{ch}{s}(rng),3,11);
            plot(t{ch}{s},f{ch}{s})
    end
    ylim([0,2]);
    title(Title{ch})
end

[VmagRatio,P,TBdot] = BdotPhaseDetect_v1(shot,Data);

%     save('Bdot_2017_11_09')
else
%     load('Bdot_2017_11_09')
end

%%
return
%%
close all
tStart = 4.1;
tEnd  = 4.5;
 
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

figure; hold on
for s = 1:length(shot)
tBdot{s} = TBdot{s}(1:length(VmagRatio{s}));    
hA(s) = plot(tBdot{s},VmagRatio{s},C{s},'marker','none','LineStyle','-');
end
xlim([tStart,tEnd])
legend(hA,num2str(shot'))
ylim([0,10])

figure; hold on
for s = 1:length(shot)
hP(s) = plot(tBdot{s}(1:length(P{s})),P{s},C{s},'marker','.','LineStyle','none');
end
xlim([tStart,tEnd])
title('Phase')
ylim([-180,180])
legend(hP,num2str(shot'))

figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.12 & tBdot{s}<=4.39);
    if 0
        plot3(tBdot{s}(rng),R(s)*ones(size(P{s}(rng))),unwrap(P{s}(rng)*pi/180)*180/pi,'marker','.','LineStyle','none')
    else
        plot3(tBdot{s}(rng),R(s)*ones(size(P{s}(rng))),P{s}(rng),'marker','.','LineStyle','none')
    end
end
xlim([tStart,tEnd])
ylim([0,15])
% zlim([-10,10])

figure; hold on
for s = 1:length(shot)
    rng = find(tBdot{s}>=4.12 & tBdot{s}<=4.39);
    plot3(tBdot{s}(rng),R(s)*ones(size(VmagRatio{s}(rng))),VmagRatio{s}(rng))
end
xlim([tStart,tEnd])
ylim([0,15])
zlim([0,7])

%%
close all
figure; hold on

n = find(R == 0);
rng = find(tBdot{n}>=4.1 & tBdot{n}<=4.5);
%plot(t_rf{n}(rng),RF{n}(rng).^2,'r')
% RF normalized:
rf_norm = -RF{n}./max(-0.97*RF{n});

% n = find(R == 8); 
% rng = find(tBdot{n}>=4.16 & tBdot{n}<=4.447);
% plot(tBdot{n}(rng),sgolay_t(VmagRatio{n}(rng),3,51),'k')

n = find(R == 6); % core
rng = find(tBdot{n}>=4.16 & tBdot{n}<=4.444);
tcore = tBdot{n}(rng);
Acore = sgolay_t(VmagRatio{n}(rng)./rf_norm(rng),3,51);
h(1) = plot(tcore,Acore,'k')

n = find(R == 11); % Edge
rng = find(tBdot{n}>=4.16 & tBdot{n}<=4.438);
tedge = tBdot{n}(rng);
Aedge = sgolay_t(VmagRatio{n}(rng)./rf_norm(rng),3,51);
h(2) = plot(tedge,Aedge,'r')

text(4.4,1.1,'Core','Color','k')
text(4.4,0.25,'Edge','Color','r')

ylim([0,1.5])
xlim([4.1,4.5])
set(gcf,'color','w')
set(gca,'PlotBoxAspectRatio',[1 0.3 1],'color','w')
box on

% Save data for usage elsewhere
if 0
save('BdotData_CoreEdge_2017_06_16','tedge','Aedge','tcore','Acore')
end