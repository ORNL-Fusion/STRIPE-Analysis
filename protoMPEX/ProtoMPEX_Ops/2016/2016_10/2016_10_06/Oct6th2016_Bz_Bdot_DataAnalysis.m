% October 6th Bdot probe XP

% =========================================================================
% B_z data analysis:
% =========================================================================

close all
clear all

% Sequence for Bz scan:
% shotlist = 1e4 + 700 + ...
%     [35,  38, 39,   41, 42,   43,   44, 45,   46, 47,   48, 49,   50, 51, 52,   53,  54, 55,  56, 57,  58, 59,   60,  61, 62, 63, 64,  65, 66, 67,  68];
% R = [8 , 8.5,  9,  9.5, 10, 10.5, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14, 14.5, 8.5,  8, 7.5,  7, 6.5,  6,  5.5, 5.5,  5,  5,  5, 4.5,  4,  4, 3.5];

shotlist = 1e4 + 700 + ...
    [35,  38, 39,   41, 42,   43,   44, 45,   46, 47,   48, 49,   50, 51, 52,   53,  54, 55,  56, 57, 59,   60, 63, 64,  65, 66];
R = [8 , 8.5,  9,  9.5, 10, 10.5, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14, 14.5, 8.5,  8, 7.5,  7,  6,  5.5,  5,  5, 4.5,  4];

%% Phase detector data:
for h = 1
Data{1} = ['\MPEX::TOP.MACHOPS1:INT_4MM_1'];
Data{2} = ['\MPEX::TOP.MACHOPS1:INT_4MM_2'];
Data{3} = ['\MPEX::TOP.MACHOPS1:INT_2MM_1'];
Data{4} = ['\MPEX::TOP.MACHOPS1:INT_2MM_2'];
[VmagRatio,P,tBdot] = BdotPhaseDetect_v1(shotlist,Data);
end
%% Preview Ne and Te data
for h = 1
if 1
    Stem = '\MPEX::TOP.';
    Branch = 'MACHOPS1:';
    RootAddress = [Stem,Branch];
    DataAddress{1} = [RootAddress,'LP_V_RAMP']; % V
    DataAddress{2} = [RootAddress,'TARGET_LP']; % I

    Config.tStart = 4.23; % [s]
    Config.tEnd = 4.32;
    Config.V_Att = 2;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
    Config.I_Att = 2;  % Output voltage of DLP box (Current) = I_att*Digitized data
    Config.V_cal = [12.05,0.205];   % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2)
    Config.I_cal = [-142.5, 1*0.015]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)
    Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
    Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
    Config.AMU = 2; % Ion mass in AMU
    Config.L_tip = 1.2/1000;
    Config.D_tip = 0.254/1000; % [m]
    Config.FitFunction = 2; 
    Config.AreaType = 1; % Cylindrical + cap
    [Ni,Te,time,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5(Config,shotlist,DataAddress);

    figure; 
    for s = 1:16
        subplot(4,4,s); hold on; yyaxis('left')
        plot(time{s},0.5*(Ni{s}{1} + Ni{s}{2}),'bl','LineWidth',1)

        ylim([0,4e19]); xlim([4.2,4.32]); box on
        line([4.2,4.32],[3e19,3e19],'color','bl')

        yyaxis('right'); hold on; 
        plot(time{s},Te{s},'r','LineWidth',1)
        ylim([0,4]); xlim([4.2,4.32]); box on
        line([4.2,4.32],[1.8,1.8],'color','r')
        title(['R = ',num2str(R(s))])
    end

    figure
    for s = 17:28
        subplot(3,4,s-16); hold on; yyaxis('left')
        plot(time{s},0.5*(Ni{s}{1} + Ni{s}{2}),'bl','LineWidth',1)
        ylim([0,4e19]); xlim([4.2,4.32]); box on
        line([4.2,4.32],[3e19,3e19],'color','bl')

        yyaxis('right'); hold on; 
        plot(time{s},Te{s},'r','LineWidth',1)
        ylim([0,4]); xlim([4.2,4.32]); box on
        line([4.2,4.32],[1.8,1.8],'color','r')
        title(['R = ',num2str(R(s))])
    end

    % COMMENT on data:
    % The data shows that ne was no that quiescent during these shots. this
    % could have been because we where operating at lower gas flow:
    % 1st puff 1.12 V and 2nd puff: 0.3 V
    % from experience we know that the best operation occurs with around
    % 1st puff: 1.2 to 1.23 V and 2nd puff: 0.7 to 0.93 V
    % Data reproduceable from R 6 to 10 cm
end
end
%% Preview Bdot probe data:
for h = 1
tStart = 4.26;
tEnd  = 4.32;
t_unwrap = 4.26;

% Phase data:
figure; for s = 1:16
    subplot(4,4,s); hold on
    yyaxis('left')
    plot(tBdot{s},P{s})
    ylim(gca,[-180,180]);
    set(gca,'YTick',-180:90:180)
    
    yyaxis('right')
    rng = find(tBdot{s}>=t_unwrap);
    plot(tBdot{s}(rng),unwrap(P{s}(rng)*pi/180)*180/pi,'r')
    xlim([tStart,tEnd])
    title(['R = ',num2str(R(s)),', ',num2str(shotlist(s))])
end
figure; for s = 17:26
    subplot(4,4,s-16); hold on
    yyaxis('left')
    plot(tBdot{s},P{s})
    ylim(gca,[-180,180]);
    set(gca,'YTick',-180:90:180)
    
    yyaxis('right')
    rng = find(tBdot{s}>=t_unwrap);
    plot(tBdot{s}(rng),unwrap(P{s}(rng)*pi/180)*180/pi,'r')
    xlim([tStart,tEnd])
    title(['R = ',num2str(R(s))])
%     ylim(gca,[-2000,2000])
end
% Magnitude data:
figure; for s = 1:16
        if shotlist(s) == 10749
        continue
    else
        subplot(4,4,s); hold on
        plot(tBdot{s},VmagRatio{s})
        ylim(gca,[0,8]);
        title(['R = ',num2str(R(s)),', ',num2str(shotlist(s))],'Fontsize',9)
        end
end
figure; for s = 17:26
    if shotlist(s) == 10766
        continue
    else
    subplot(4,4,s-16); hold on
    plot(tBdot{s},VmagRatio{s})
    ylim(gca,[0,8]);
    title(['R = ',num2str(R(s)),', ',num2str(shotlist(s))],'Fontsize',9)
    end
end
end
%% Plot3 data
for h = 1;
tStart = 4.23;
tEnd  = 4.32;
rng = find(tBdot{1}>=tStart & tBdot{1}<=tEnd); 

badshots = 10700 + [49,66];

for s = 1:length(shotlist)
    unP{s} = unwrap(P{s}*pi/180)*180/pi;
end

Offset = zeros(size(P{1}));
for k = 1:length(unP)
            Offset(k) =  round(mean(unP{k}(rng))/360)*360; % [Deg]
end

figure; hold on
for s = 1:length(shotlist)
    if sum(shotlist(s) ==  badshots)
        continue
    else
    hP(s) = plot3(tBdot{s},R(s)*ones(size(P{s})), unP{s} - Offset(s));
    end
end
title('Phase')
zlim([-360,360])
xlim([tStart,tEnd])
view([30,30])

figure; hold on
for s = 1:length(shotlist)
    if sum(shotlist(s) ==  badshots)
        continue
    else
    hM(s) = plot3(tBdot{s},R(s)*ones(size(VmagRatio{s})),sgolay_t(VmagRatio{s},3,11));
    end
end
title('Magnitude')
zlim([0,6])
xlim([tStart,tEnd])
view([30,30])

end
%% steady state data:
for h = 1;
    close all
tStart = 4.285;
tEnd  = 4.3;
xOffset = -0*7.5;

rng = find(tBdot{1}>=tStart & tBdot{1}<=tEnd); 
for s = 1:length(shotlist)
    if sum(shotlist(s) ==  badshots)
        p_Bt(s) = NaN;
        dp_Bt(s) =NaN;
        a_Bt(s) = NaN;
        da_Bt(s) = NaN;
    else
        p_Bt(s) = mean(unP{s}(rng) - Offset(s));
        dp_Bt(s) = std(unP{s}(rng) - Offset(s),1);
        a_Bt(s) = mean(VmagRatio{s}(rng));
        da_Bt(s) = std(VmagRatio{s}(rng));
    end
end
figure; subplot(2,1,1);hold on
errorbar((R + xOffset),p_Bt,dp_Bt,'k.')
title('Phase of B_{\theta}')
xlim([0,15]+ xOffset)
ylim([-180,180]); set(gca,'YTick',-180:60:180); grid on; box on
ax(1) = gca;

subplot(2,1,2); hold on
errorbar((R + xOffset),a_Bt,da_Bt,'k.')
title('Magnitude of B_{\theta}')
xlim([0,15] + xOffset); set(gca,'XTick',-15:3:15); grid on; box on
ylim([0,5])
ax(2) = gca;

set(gcf,'color','w')
set(ax,'PlotBoxAspectRatio',[1 0.8 1])
end