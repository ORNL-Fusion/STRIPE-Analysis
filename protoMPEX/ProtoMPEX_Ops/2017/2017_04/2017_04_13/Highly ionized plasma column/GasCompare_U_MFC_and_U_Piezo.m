mdsconnect('mpexserver')
if 1
    close all
clear all
% =========================================================================
% final best shots for paper with calibration
shotlist = [13962,21269,17536];
Calshot =  [13968,14120,14120];
% =========================================================================
end

PressureInPascal = 1;
RemoveOffsetStartRF = 1;

address{1} = '\MPEX::TOP.MACHOPS1:PG1'; % PG9.5
address{2} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5

if strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG3')
     f = 2;
     T = 'P3';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG1')
     f = 2;
     T = 'P4';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG4')
     f = 10;
     T = 'P2';
elseif strcmp(address{1},'\MPEX::TOP.MACHOPS1:PG2')
     f = 2;
     T = 'P1';
 end
% f = 2; % use 10 for PG4 and 2 for all other

[PG1,t1] = my_mdsvalue_v2(shotlist,address(1));
[PG2,t2] = my_mdsvalue_v2(shotlist,address(2));
[PG3,t3] = my_mdsvalue_v2(shotlist,address(3));
[PG4,t4] = my_mdsvalue_v2(shotlist,address(4));

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

for s = 1 :length(shotlist)
            L1 = length(PG1{1});
        L2 = length(PG1cal{1});
        if L1>L2
            Ldata = L2;
        else
            Ldata = L1;
        end

        if length(Calshot)>1
            k = s;
        else
            k = 1;
        end
            P{1}{s} = ( PG1{s}(1:Ldata)-PG1cal{k}(1:Ldata) )*2*fctr;
            P{2}{s} = ( PG2{s}(1:Ldata)-PG2cal{k}(1:Ldata) )*2*fctr;
            P{3}{s} = ( PG3{s}(1:Ldata)-PG3cal{k}(1:Ldata) )*2*fctr;
            P{4}{s} = ( PG4{s}(1:Ldata)-PG4cal{k}(1:Ldata) )*10*fctr;
end
T = {'4','1','3','2'};

for p = 1:4
subplot(2,2,p); hold on
for s = 1:length(shotlist)
    
    if s == 31
        s
    end
    Time = t1{s}(1:length(P{p}{s}));
    
    if RemoveOffsetStartRF
        rngMean = find(Time>3.6 & Time <4.0);
%         rngMean = find(Time>2 & Time <3.5);

        Offset{p}{s} = mean(P{p}{s}(rngMean));
    else 
        Offset{p}{s} = 0;
    end
    
    h(s) = plot(Time,P{p}{s} - Offset{p}{s},C{s});
    set(h(s),'LineWidth',2)
    title(['P',T{p}])
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
            ylim([-0.5,5]*fctr)
        case '4.5'
            ylim([-0.5,35]*fctr)  
    end
end
end
legend(h,{num2str(shotlist')},'location','NorthWest')
set(gcf,'color','w')

figure; 
hold on
p = 3;
for s = length(shotlist):-1:1
Time = t1{s}(1:length(P{p}{s}));
if RemoveOffsetStartRF
        rngMean = find(Time>3.8 & Time <4.0);
%         rngMean = find(Time>2 & Time <3.5);

        Offset{p}{s} = mean(P{p}{s}(rngMean));
    else 
        Offset{p}{s} = 0;
end
    h(s) = plot(Time,P{p}{s} - Offset{p}{s},C{s});
    set(h(s),'LineWidth',2)
    title(['P',T{p}])
    box on
    grid on
    if PressureInPascal
        ylabel('[Pa]')
    else
        ylabel('mTorr')
    end
    
end
legend(h,'G2','G3','location','NorthWest')
xlim([4,4.8])
ylim([-0.5,3]*fctr)
xlabel('time [s]')
set(gcf,'position',[565  366  355  251],'color','w')