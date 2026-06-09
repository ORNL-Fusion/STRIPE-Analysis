% Load dataset 6 from 2018_04_13
clear all
close all

DataName = 'DataSet_1_to_6_IFP_2018_04_13_RawData';
load(DataName)
struct2table(AddressTable)
RawData

%% Data subset
Shotlist = 21000 + [238,239,240,241,245,246,247];
v        =         [10 ,9  ,8  ,8  ,7  ,6  ,  5];
[~,b] = sort(v);

RawData_2 = GetRawDataSubset(Shotlist,RawData);
% Ion current from probe:
I = RawData_2.I;
t_I   = RawData_2.t_I;
% Neutral gas pressure:
for s = b
    P{1}{s} = RawData_2.PG1{s}*2/7.5; 
    P{2}{s} = RawData_2.PG2{s}*2/7.5;
    P{3}{s} = RawData_2.PG3{s}*2/7.5;
    P{4}{s} = RawData_2.PG4{s}*2/7.5;
end
t_P = RawData_2.t_PG1{1}(1:end-1);

figure;

for s = 1:length(v)
    subplot(3,3,s)
    plot(t_I{b(s)}(1:end-1),I{b(s)})
    ylim([-0.5,0.5])
    xlim([4,5])
    title(['R: ',num2str(v(b(s)))])
end

Config.tStart          = 4.17;
Config.tEnd            = 4.7 ;
Config.FitFunction     = 2   ;
Config.Center_V        = 1   ;
Config.Center_I        = 1   ;
Config.FilterDataInput = 1   ;
Config.SGF             = 11  ;
Config.TimeMode        = 1   ;
Config.AMU             = 2   ;
Config.AreaType        = 2   ;
Config.L_tip           = 0   ;
Config.D_tip           = ((2^(1/4))*1.22)/1000; % [m]
Config.V_Att           = 1   ;
Config.I_Att           = 1   ;
Config.V_cal           = [(0.46e-3)^-1,0]; % Voltage output of DLP
Config.I_cal           = [-1,0]; % Current output of DLP

DLPData = DLP_fit_V6(Config,RawData_2);
Isat = DLPData.Isat_m;
t_Isat = DLPData.time;

%% Steady-state values:
% Ion flux data
close all
figure
A_fp = 0.25*sqrt(2)*pi*(1.22/1000)^2;

for s = b
    t1{s} = 4.62;
    t2{s} = 4.64;
end

t1{2} = 4.54;
t2{2} = 4.56;

for s = b
    subplot(3,3,s); hold on
    
    % For each sweep, select only the good fits
    GoodFits{s} = DLPData.GlitchFlag{s} == 0 & DLPData.StdResNorm{s}<=0.25...
        & DLPData.Isat_m{s}>0 & DLPData.Isat_m{s}>0 & DLPData.Te{s}<=10;
    
    IsatGood{s} = sgolay_t(Isat{s}(GoodFits{s}),3,7);
    t_IsatGood{s} = t_Isat{s}(GoodFits{s});
    RFpwr(s) = v(s)*8;
    
    plot(t_I{s}(1:end-1),abs(I{s}))
    plot(t_IsatGood{s},IsatGood{s})

    ylim([0,0.4])
    xlim([4,5])
    title(['R: ',num2str(v(s)),' shot: ',num2str(DLPData.shot(s))])
        
    rngT{s}{1} = [t_IsatGood{s}>=4.20 & t_IsatGood{s}<=4.23]; 
    rngT{s}{2} = [t_IsatGood{s}>=4.5 & t_IsatGood{s}<=4.6];  
    rngT{s}{3} = [t_IsatGood{s}>= t1{s} & t_IsatGood{s}<= t2{s}]; 
    
   for p = 1:3
         t_Ifp{p}(s)  = mean(t_IsatGood{s}(rngT{s}{p}));
         Ifp{p}(s)    = mean(IsatGood{s}(rngT{s}{p}));
         dIfp{p}(s)   = std(IsatGood{s}(rngT{s}{p}),1,2);
         flux{p}(s)  = Ifp{p}(s)/(e_c*A_fp);
         dflux{p}(s) = dIfp{p}(s)/(e_c*A_fp);
        
       plot(t_IsatGood{s}(rngT{s}{3}),IsatGood{s}(rngT{s}{3}),'ro')
    end
end

figure; 
errorbar(RFpwr,(1e-24)*flux{3},(1e-24)*dflux{3})
xlim([30,100])
ylim([0,2])

%% Neutral gas data

A_fp = 0.25*sqrt(2)*pi*(1.2/1000)^2;
Ac = pi*(2.0/100)^2;

for s = b
    t1{s} = 4.62;
    t2{s} = 4.64;
end

for s = b
    for m = 1:4
        rngMean = find(t_P>3.6 & t_P <4.0);
        Offset{m}{s} = mean(P{m}{s}(rngMean));
        Pg{m}{s} = P{m}{s} - Offset{m}{s};
        
    rngT{s}{1} = [t_IsatGood{s}>=4.20 & t_IsatGood{s}<=4.23]; 
    rngT{s}{2} = [t_IsatGood{s}>=4.5 & t_IsatGood{s}<=4.6];  
    rngT{s}{3} = [t_IsatGood{s}>= t1{s} & t_IsatGood{s}<= t2{s}]; 
        
       for p = 1:3
         Pg_m{p}{m}(s)  = mean(interp1(t_P,Pg{m}{s},t_IsatGood{s}(rngT{s}{p})));
         dPg_m{p}{m}(s) = std(interp1(t_P,Pg{m}{s},t_IsatGood{s}(rngT{s}{p})),1,2);

       end
    end
end

figure; hold on
for s = b;
    h(s) = plot(t_P,Pg{1}{s});
end
legend(h,num2str(v'))


figure;
yyaxis right
ax(1) = gca;
plot(RFpwr([1:3,5:end]),Pg_m{3}{1}([1:3,5:end]),'ko-')
ylim(ax(1),[0,0.28])
xlim(ax(1),[30,90])
ax(1).YColor = 'k';
t = ylabel('$P_{Target}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

yyaxis left
ax(2) = gca;
errorbar(RFpwr,(1e-24)*flux{3},2*(1e-24)*dflux{3},'rsq-')
xlim(ax(2),[30,90])
ylim(ax(2),[0,2])
ax(2).YColor = 'r';
t = ylabel('${\Gamma}_+$ $[m^{-2}s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

t = xlabel('RF power $[kW]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')
box on
set(gcf,'position',[360.3333  389.6667  560.0000  228.0000])
set(gcf,'position',[2.3333  191.0000  458.6667  284.0000])



figure;
yyaxis right
ax(1) = gca;
plot(RFpwr([1:3,5:end]),Pg_m{3}{1}([1:3,5:end]),'ko-')
ylim(ax(1),[0,0.28])
xlim(ax(1),[30,90])
ax(1).YColor = 'k';
t = ylabel('$P_{Target}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

yyaxis left
ax(2) = gca;
R = 2.1/100;
A = pi*R^2;
Sion = 0.3*flux{3}*A;
dSion = 0.3*dflux{3}*A;
errorbar(RFpwr,Sion,3*dSion,'rsq-')
xlim(ax(2),[30,90])
ylim(ax(2),[0,8]*1e20)
ax(2).YColor = 'r';
t = ylabel('$S_+$ $[s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

t = xlabel('RF power $[kW]$ '); t.Interpreter = 'latex'; t.FontSize = 13;

set(gcf,'color','w')
box on
set(gcf,'position',[360.3333  389.6667  560.0000  228.0000])
set(gcf,'position',[2.3333  191.0000  458.6667  284.0000])

%%
figure
[a,b] = sort(RFpwr);
subplot(2,2,1); hold on
for p = 3
h(p) = errorbar(RFpwr(b),flux{p}(b),dflux{p}(b));
end
xlim([4,11])
t = title('${\Gamma}_+$ $[m^{-2}s^{-1}]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,2); hold on
for p = 3
h(p) = errorbar(RFpwr(b),pg1_ss{p}(b),dpg1_ss{p}(b),'ko-')
end
xlim([4,11])
ylim([0,0.3])
t = title('$P_{Target}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,3); hold on
for p = 3
h2 = plot(RFpwr(b),pg2_ss{p}(b),'ko-')
h4 = plot(RFpwr(b),pg4_ss{p}(b),'ro-')
end
legend([h2,h4],'PG 2.5','PG 4.5')
xlim([4,11])
ylim([0,20])
t = title('$P_{Source}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

subplot(2,2,4); hold on
for p = 3
h(p) = plot(RFpwr(b),pg3_ss{p}(b),'ko-')
end
xlim([4,11])
ylim([-0.5,2.5])
t = title('$P_{CC}$ $[Pa]$ '); t.Interpreter = 'latex'; t.FontSize = 13;
box on

set(gcf,'color','w')