% close all
% clear all

% shotlist = 14000 + [104,367];
% shotlist = [14000 +[105:107]]
mdsconnect('mpexserver')

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
Calshot = 14713;
[PG1cal,t0cal] = my_mdsvalue_v2(Calshot,address(1));
[PG2cal,t1cal] = my_mdsvalue_v2(Calshot,address(4));
[PG3cal,t2cal] = my_mdsvalue_v2(Calshot,address(3));
[PG4cal,t4cal] = my_mdsvalue_v2(Calshot,address(4));


%%
figure;
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

for s = 1 :length(shotlist)
    P{1}{s} = (PG1{s}-PG1cal{1})*2;
    P{2}{s} = (PG2{s}-PG2cal{1})*2;
    P{3}{s} = (PG3{s}-PG3cal{1})*2;
    P{4}{s} = (PG4{s}-PG4cal{1})*10;
    
%     P{1}{s} = (PG1{s})*2;
%     P{2}{s} = (PG2{s})*2;
%     P{3}{s} = (PG3{s})*2;
%     P{4}{s} = (PG4{s})*10;
end
T = {'9.5','2.5','6.5','4.5'};

for p = 1:4
subplot(2,2,p); hold on
for s = 1:length(shotlist)
    h(s) = plot(t1{s}(1:length(P{p}{s})),P{p}{s},C{s});
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
