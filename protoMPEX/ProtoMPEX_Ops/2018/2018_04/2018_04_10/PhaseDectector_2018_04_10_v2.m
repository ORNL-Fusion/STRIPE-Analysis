% 2017_10_12, Bdot probe scan at 4kA

close all
clear all
LoadData = 1; 

if LoadData
shot = 21000 + [64,65,66];
shot = 21000 + [68,69,70,71];
shot = 21000 + [75,76];
shot = 21000 + [79,80];

%Start of radial scan
shot = 21000 + [80,81,82,84,85];
shot = 21000 + [84,85]; % shot 86 does not work
shot = 21000 + [87,88];
shot = 21000 + [88,89,90];
shot = 21000 + [92,93];

shot = 21000 + [84  ,  85,87  , 88  ,89 ,90 ];
R    =         [-2.5,-2.0,-1.0, -0.5,0.0,0.5];
% shot = 21000 + [95];

% 180 degree phase test
% shot = 21000 + [82,83];

% coil S1: 9.5 , B_phi_z
Data{1} = ['\MPEX::TOP.MACHOPS1:COIL1_1']; % Vp0
Data{2} = ['\MPEX::TOP.MACHOPS1:COIL1_2']; % Vp90
Data{3} = ['\MPEX::TOP.MACHOPS1:COIL1_3']; % Vmag
Title{1} = 'S1,Vp0';
Title{2} = 'S1,Vp90';
Title{3} = 'S1,Vmag0';

% coil S2: 9.5 , Br
Data{4} = ['\MPEX::TOP.MACHOPS1:COIL1_4']; % Vp0
Data{5} = ['\MPEX::TOP.MACHOPS1:COIL1_5']; % Vp90
Data{6} = ['\MPEX::TOP.MACHOPS1:COIL1_6']; % Vmag
Title{4} = 'S2,Vp0';
Title{5} = 'S2,Vp90';
Title{6} = 'S2,Vmag0';

% coil S3: 10.5, B_phi_z
Data{7} = ['\MPEX::TOP.MACHOPS1:COIL1_7']; % Vp0
Data{8} = ['\MPEX::TOP.MACHOPS1:COIL1_8']; % Vp90
Data{9} = ['\MPEX::TOP.MACHOPS1:COIL1_9']; % Vmag
Title{7} = 'S3,Vp0';
Title{8} = 'S3,Vp90';
Title{9} = 'S3,Vmag0';

% coil S4: 10.5, Br
Data{10} = ['\MPEX::TOP.MACHOPS1:COIL1_10']; % Vp0
Data{11} = ['\MPEX::TOP.MACHOPS1:COIL1_11']; % Vp90
Data{12} = ['\MPEX::TOP.MACHOPS1:COIL1_12']; % Vmag
Title{10} = 'S4,Vp0';
Title{11} = 'S4,Vp90';
Title{12} = 'S4,Vmag0';

Stem = '\MPEX::TOP.'; Branch = 'MACHOPS1:'; RootAddress = [Stem,Branch];
% DA{1} = [RootAddress,'RF_FWD_PWR'];
% [RF,t_rf]   = my_mdsvalue_v2(shot,DA(1))
% DA{1} = [RootAddress,'EA_CURRENT'];
% [Isat,t_isat]   = my_mdsvalue_v2(shot,DA(1));

% Title{1} = 'Vp0';
% Title{2} = 'Vp90';
% Title{3} = 'Vmag0';
% Title{4} = 'None';

figure;
    for ch = 1:12 % For all channels on digitizer
        [f{ch},t{ch}] = my_mdsvalue_v2(shot,Data(ch));
        % f{channel}{shot}
        subplot(4,3,ch);
        hold on
        for s = 1:length(shot)
                rng = find(t{ch}{s}>= 4.1 & t{ch}{s}<= 4.7);
                t{ch}{s} = t{ch}{s}(rng);
                if isempty(rng)
                    rng = find(t{ch}{s+1}>= 4.1 & t{ch}{s+1}<= 4.7);
                    t{ch}{s} = t{ch}{s+1}(rng);
                end
                f{ch}{s} = f{ch}{s}(rng);%sgolay_t(f{ch}{s}(rng),3,11);
                h(ch,s) = plot(t{ch}{s},f{ch}{s});
        end
        ylim([0,2]);
        title(Title{ch})
%         legend(h(ch,:),num2str(shot))
        xlim([4,4.8])
    end

[VmagRatio_S1,P_S1,tBdot_S1] = BdotPhaseDetect_v3(shot,Data(1:3));
[VmagRatio_S2,P_S2,tBdot_S2] = BdotPhaseDetect_v3(shot,Data(4:6));
[VmagRatio_S3,P_S3,tBdot_S3] = BdotPhaseDetect_v3(shot,Data(7:9));
[VmagRatio_S4,P_S4,tBdot_S4] = BdotPhaseDetect_v3(shot,Data(10:12));

% for c = 1:4
%         
% end

%     save('Bdot_2017_10_12')
else
%     load('Bdot_2017_10_12')
end

% return
%%
tStart = 4.1;
tEnd   = 4.7;
close all
figure
N = ceil(sqrt(length(shot)));
subplot(2,2,1); 
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(VmagRatio_S1{p})),tBdot_S1{p}(1:length(VmagRatio_S1{p})),VmagRatio_S1{p})
ylim([tStart,tEnd])
zlim([0,7])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end
subplot(2,2,2); 
hold on
for p = 1:length(shot)
h2(p) = plot3(R(p)*ones(size(VmagRatio_S2{p})),tBdot_S2{p}(1:length(VmagRatio_S2{p})),VmagRatio_S2{p})
ylim([tStart,tEnd])
zlim([0,7])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end
subplot(2,2,4); 
hold on
for p = 1:length(shot)
h2(p) = plot3(R(p)*ones(size(VmagRatio_S3{p})),tBdot_S4{p}(1:length(VmagRatio_S4{p})),VmagRatio_S4{p})
ylim([tStart,tEnd])
zlim([0,15])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end

figure
N = ceil(sqrt(length(shot)));
subplot(2,2,1); 
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(P_S1{p})),tBdot_S1{p}(1:length(P_S1{p})),P_S1{p}-0*mean(P_S1{p}(1:100)))
ylim([tStart,tEnd])
zlim([-250,250])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end
subplot(2,2,2); 
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(P_S2{p})),tBdot_S2{p}(1:length(P_S2{p})),P_S2{p}-0*mean(P_S2{p}(1:100)))
ylim([tStart,tEnd])
zlim([-250,250])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end
subplot(2,2,4); 
hold on
for p = 1:length(shot)
plot3(R(p)*ones(size(P_S4{p})),tBdot_S4{p}(1:length(P_S4{p})),P_S4{p}-0*mean(P_S4{p}(1:100)))
ylim([tStart,tEnd])
zlim([-250,250])
title([num2str(shot(p)),' ,Vmag'])
grid on
view([-60,40])
end


return
%%



figure
N = ceil(sqrt(length(shot)));
for p = 1:length(shot)
subplot(N,N,p); hold on
plot(tBdot_S1{p}(1:length(VmagRatio_S1{p})),P_S1{p}-mean(P_S1{p}(1:100)),'k.')
plot(tBdot_S2{p}(1:length(VmagRatio_S2{p})),P_S2{p}-mean(P_S2{p}(1:100)),'r.')
plot(tBdot_S4{p}(1:length(VmagRatio_S4{p})),P_S4{p}-mean(P_S4{p}(1:100)),'g.')

xlim([tStart,tEnd])
ylim([-180,180])
title([num2str(shot(p)),' ,Phase'])
grid on
end

return
%%
close all
figure
N = ceil(sqrt(length(shot)));
for p = 1:length(shot)
subplot(N,N,p)
plot(tBdot{p}(1:length(VmagRatio{p})),VmagRatio{p},'r.')
xlim([tStart,tEnd])
ylim([0,10])
title(num2str(shot(p)))
grid on
end

figure
N = ceil(sqrt(length(shot)));
for p = 1:length(shot)
subplot(N,N,p)
plot(tBdot{p}(1:length(VmagRatio{p})),P{p},'k.')
xlim([tStart,tEnd])
ylim([-180,180])
title(num2str(shot(p)))
grid on
end
%%
close all
tStart = 4.1;
tEnd  = 4.5;
 
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

figure; hold on
for s = 1:length(shot)
tBdot{s} = tBdot{s}(1:length(VmagRatio{s}));    
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