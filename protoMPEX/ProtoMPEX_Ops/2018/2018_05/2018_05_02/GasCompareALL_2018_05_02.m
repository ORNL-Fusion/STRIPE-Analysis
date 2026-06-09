

mdsconnect('mpexserver')
if 0
    close all
clear all
% =========================================================================
% 2018_04_30
% Really good shot showing low neutral pressure at spool 6.5
% TR2 = 180 A, PS1/2 = 4.5 kA, spool 1.5 puffer at 3.5 v level
% shotlist = 21000 + [269]

% Changing timming of gas puff relative to RF pulse
shotlist = 21000 + [717  ,719  ,720  ,721  ,722  ,723  ,724  ,725  ,727  ];
dt       =         [3.985,3.980,3.975,3.970,3.965,3.955,3.945,3.935,3.935];

% Comparing two traces at differnt times, 
shotlist = 21000 + [719  ,727  ];
dt       =         [3.985,3.935]; % the difference is 50 ms
% the data suggests that we need to scan the effect of the gas every 50 ms
% to study the effect of the neutral gas puff timing relative to the RF

% =========================================================================
% 2018_05_01
% Changing the timming of the gas puff:
shotlist = 21000 + [739,742,743,744,745];
dt       =         [0  ,0  ,100,200,300];

shotlist = 21000 + [742,743,744,745,830];
dt       =         [0  ,100,200,300,400];

% Changing the amplitude of the puff:
% shotlist = 21000 + [746,747,748,826,828];
% dt       =         [3.5,5  ,5  ,5  ,6.5];

% Closing the central chamber pump:
% shotlist = 21000 + [742,831,832:833];
% dt       =         [];

end

PressureInPascal = 1;
RemoveOffsetStartRF = 1;

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % PG9.5
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5

if strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG3')
     f = 2;
     T = 'PG6.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG1')
     f = 2;
     T = 'PG9.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG4')
     f = 10;
     T = 'PG4.5';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG2')
     f = 2;
     T = 'PG2.5';
 end
% f = 2; % use 10 for PG4 and 2 for all other

[PG1,t1] = my_mdsvalue_v2(shotlist,address(1));
[PG2,t2] = my_mdsvalue_v2(shotlist,address(2));
[PG3,t3] = my_mdsvalue_v2(shotlist,address(3));
[PG4,t4] = my_mdsvalue_v2(shotlist,address(4));

% Calibration shot:
Calshot = 14971;
Calshot = 17266; % 2017_11_01
Calshot = 17539; % 2017_11_10
Calshot = 17830; % 2017_11_16
Calshot = 19006; % 2018_01_16
Calshot = 19752; % 2018_01_16
Calshot = 20493; % 2018_03_26, 1 sec RF pulse
Calshot = 20839; % 2018_04_03, 1 sec RF pulse
Calshot = 20887; % 2018_04_05, 2.5 kA case, 1 sec RF pulse
Calshot = 20908; % 2018_04_05, 2.0 kA case, 1 sec RF pulse
% Calshot = 14120;
% 
[PG1cal,t0cal] = my_mdsvalue_v2(Calshot,address(1));
[PG2cal,t1cal] = my_mdsvalue_v2(Calshot,address(4));
[PG3cal,t2cal] = my_mdsvalue_v2(Calshot,address(3));
[PG4cal,t4cal] = my_mdsvalue_v2(Calshot,address(4));


%%
figure;
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

if PressureInPascal
    fctr = 0.1333;
else
    fctr = 1;
end

L1 = length(PG1{1});
L2 = length(PG1cal{1});
if L1>L2
    Ldata = L2;
else
    Ldata = L1;
end

for s = 1 :length(shotlist)
    P{1}{s} = ( PG1{s}(1:Ldata)-PG1cal{1}(1:Ldata) )*2*fctr;
    P{2}{s} = ( PG2{s}(1:Ldata)-PG2cal{1}(1:Ldata) )*2*fctr;
    P{3}{s} = ( PG3{s}(1:Ldata)-PG3cal{1}(1:Ldata) )*2*fctr;
    P{4}{s} = ( PG4{s}(1:Ldata)-PG4cal{1}(1:Ldata) )*10*fctr;
    
end
T = {'9.5','2.5','6.5','4.5'};

for p = 1:4
subplot(2,2,p); hold on
for s = 1:length(shotlist)
    
    if s == 31
        s
    end
    Time = t1{s}(1:length(P{p}{s}));
    
    if RemoveOffsetStartRF
        rngMean = find(Time>3.6 & Time <4.0);
        rngMean = find(Time>2 & Time <3.5);

        Offset{p}{s} = mean(P{p}{s}(rngMean));
    else 
        Offset{p}{s} = 0;
    end
    
    h(s) = plot(Time,P{p}{s} - Offset{p}{s},C{s});
    set(h(s),'LineWidth',2)
    title(['PG ',T{p}])
    box on
    grid on
    if PressureInPascal
        ylabel('[Pa]')
    else
        ylabel('mTorr')
    end
    xlim([3.5,6])
    switch T{p}
        case '9.5'
            ylim([-0.5,4]*fctr)
        case '2.5'
            ylim([-0.5,25]*fctr)
        case '6.5'
            ylim([-0.5,8]*fctr)
        case '4.5'
            ylim([-0.5,35]*fctr)  
    end
end
end
legend(h,{num2str(shotlist')},'location','NorthWest')
set(gcf,'color','w')
