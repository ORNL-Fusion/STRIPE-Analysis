close all
clear all

mdsconnect('mpexserver')

shotlist = [11000 + [67,83,597,618,662,777,842]]; 

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
Calshot = 11832;

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
        Offset = mean(P{p}{s}(rngMean));
    else 
        Offset = 0;
    end
    
    h(s) = plot(Time,P{p}{s} - Offset,C{s});
    set(h(s),'LineWidth',2)
    title(['PG ',T{p}])
    box on
    grid on
    if PressureInPascal
        ylabel('[Pa]')
    else
        ylabel('mTorr')
    end
    xlim([3.8,6])
    switch T{p}
        case '9.5'
            ylim([-0.5,4]*fctr)
        case '2.5'
            ylim([-0.5,25]*fctr)
        case '6.5'
            ylim([-1.5,20]*fctr)
        case '4.5'
            ylim([-0.5,35]*fctr)  
    end
end
end
legend(h,{num2str(shotlist')},'location','NorthWest')
set(gcf,'color','w')
