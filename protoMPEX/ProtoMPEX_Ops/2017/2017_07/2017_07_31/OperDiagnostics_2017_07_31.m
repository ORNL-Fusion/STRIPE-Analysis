clear all
mdsconnect('mpexserver')

% Add shot numnbers
% Add ECH traces

% =========================================================================
% First few shots on July 20th 2017
shotlist = [15000 + [841]];

% close all
TR2_filt =         ones(size(shotlist)); %[1    ,1    ,1    ,1    ,1    ,1    ,1    ];
TR2_ffrm =         ones(size(shotlist))*15; %[15   , 15  ,15   ,15   ,15   ,15   ,15   ];
Calshot  =         ones(size(shotlist))*15405; %[15405,15460,15460,15460,15460,15460,15460];
PGcenter =         [];
tubes    = []; % choose 4 tubes at most
% =========================================================================

% Dalpha at 6.5, 
% SXR signals, 

% Gather data
% TR1:
address{1} = '\MPEX::TOP.MACHOPS1:TRANS_I'; % TR2 current
[TR1,t_tr2] = my_mdsvalue_v2(shotlist,address(1));

% RF trace (raw signal level)
address{2} = '\MPEX::TOP.MACHOPS1:RF_FWD_PWR'; % RF power trace
[RF,t_rf] = my_mdsvalue_v2(shotlist,address(2));
address{2} = '\MPEX::TOP.MACHOPS1:PWR_28GHZ'; % 28 GHz power trace
[ECH,t_ech] = my_mdsvalue_v2(shotlist,address(2));

% PS1 and PS2
address{3} = '\MPEX::TOP.MACHOPS1:PS1_I'; % PS1 current
[PS1,t_ps1] = my_mdsvalue_v2(shotlist,address(3));
address{4} = '\MPEX::TOP.MACHOPS1:PS2_I'; % PS1 current
[PS2,t_ps2] = my_mdsvalue_v2(shotlist,address(4));

% Mass flow controller voltage
address{5} = '\MPEX::TOP.MACHOPS1:MFC_FLOW_D2'; % Mass flow controller voltage
[MFC,t_mfc] = my_mdsvalue_v2(shotlist,address(5));

% PG data
address{6} = '\MPEX::TOP.MACHOPS1:PG1'; % PG9.5
address{7} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{8} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{9} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5

for s = 1:length(shotlist)
    for p = 1:4
        [PG{p}{s}   ,t_pg{p}{s} ] = my_mdsvalue_v2(shotlist(s),address(5+p));
        [PGcal{p}{s},t_cal{p}{s}] = my_mdsvalue_v2(Calshot(s) ,address(5+p)); 
    end
        PGc{1}{s} = (cell2mat(PG{1}{s})-cell2mat(PGcal{1}{s}))*2;
        t_pgc{1}{s} = cell2mat(t_pg{1}{s});
        PGc{2}{s} = (cell2mat(PG{2}{s})-cell2mat(PGcal{2}{s}))*2;
        t_pgc{2}{s} = cell2mat(t_pg{2}{s});
        PGc{3}{s} = (cell2mat(PG{3}{s})-cell2mat(PGcal{3}{s}))*2;
        t_pgc{3}{s} = cell2mat(t_pg{3}{s});        
        PGc{4}{s} = (cell2mat(PG{4}{s})-cell2mat(PGcal{4}{s}))*10;
        t_pgc{4}{s} = cell2mat(t_pg{4}{s});       
end

% Filterscope data:

%% Plot data
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

figure;
subplot(2,2,1); hold on
for s = 1:length(shotlist)
    if TR2_filt(s)
        fs = sgolay_t(TR1{s},3,TR2_ffrm(s));
        plot(t_tr2{s}(1:length(TR1{s})),fs*1000,C{s})
    else
        plot(t_tr2{s}(1:length(TR1{s})),TR1{s}*1000,C{s})
    end
end
ylim([100,300])
xlim([4.1,4.5])
grid on
box on
title('TR2')

subplot(2,2,3); hold on
for s = 1:length(shotlist)
    h(s) = plot(t_rf{s}(1:length(RF{s})),RF{s}.^2,C{s})
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}/10,C{s})

end
ylim([0,1])
xlim([4.1,4.5])
grid on
box on
title('RF')
legend(h,{num2str(shotlist')},'location','NorthWest')

subplot(2,2,2); hold on
for s = 1:length(shotlist)
    h(1) = plot(t_ps1{s}(1:length(PS1{s})),PS1{s},C{s})
    h(2) = plot(t_ps2{s}(1:length(PS2{s})),PS2{s},C{s})
end
plot(t_rf{1}(1:length(RF{1})),4 + 0.2*RF{1}.^2,C{1})
set(h(1),'LineWidth',2)
ylim([3,5])
xlim([4,5])
grid on
box on
title('PS1 [kA]')
legend(h,'PS1','PS2')

subplot(2,2,4); hold on
for s = 1:length(shotlist)
    plot(t_mfc{s}(1:length(MFC{s})),MFC{s},C{s})
end
plot(t_rf{1}(1:length(RF{1})),0.2*RF{1}.^2,C{1})
ylim([0,1.2])
xlim([3.85,5])
grid on
box on
title('MFC {D2}')

set(gcf,'color','w')

% Gas data
figure;
T = {'9.5','2.5','6.5','4.5'};
for p = 1:4
subplot(2,2,p); hold on
for s = 1:length(shotlist)
    h(s) = plot(t_pgc{p}{s}(1:length(PGc{p}{s})),PGc{p}{s},C{s});
    set(h(s),'LineWidth',2)
    title(['PG ',T{p}])
    box on
    grid on
    ylabel('mTorr')
    xlim([4,4.7])
    switch T{p}
        case '9.5'
            ylim([-0.5,3])
        case '2.5'
            ylim([-0.5,15])
        case '6.5'
            ylim([-0.5,5])
        case '4.5'
            ylim([-0.5,15])  
    end
end
end
legend(h,{num2str(shotlist')},'location','NorthWest')
set(gcf,'color','w')



