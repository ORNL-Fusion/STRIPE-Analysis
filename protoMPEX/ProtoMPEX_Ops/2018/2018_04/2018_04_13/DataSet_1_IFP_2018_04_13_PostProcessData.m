% Load dataset 6 from 2018_04_13
clear all
close all

DataName = 'DataSet_1_to_6_IFP_2018_04_13_RawData';
load(DataName)
struct2table(AddressTable)
RawData

%% Data subset

Shotlist = [21000 + [225,227,228,230,231,232,233,234,235,236,237]];
r        =          [5.5,5.0,4.5,4.0,3.5,6.5,7.0,7.5,8.0,8.5,6.0] ;
r        =          [5.5,5.0,4.5,4.0,3.5,6.5,7.5,7.0,8.0,8.5,6.0] ;

RawData_1 = GetRawDataSubset(Shotlist,RawData);

figure;
[~,b] = sort(r);

for s = 1:length(r)
    subplot(4,3,s)
    plot(RawData_1.t_I{b(s)}(1:end-1),RawData_1.I{b(s)})
    ylim([-0.4,0.4])
    xlim([4,5])
    title(['R: ',num2str(r(b(s)))])
end

Config.tStart          = 4.15;
Config.tEnd            = 4.8 ;
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

DLPData = DLP_fit_V6(Config,RawData_1);

%%
figure; hold on
for s = 1:length(r)
    plot(DLPData.Isat{s}{1},'k')
    plot(DLPData.Isat{s}{2},'r')
    plot(DLPData.Ip{s},'g')
end
ylim([-0.5,0.5])

figure
A_fp = 0.25*sqrt(2)*pi*(1.22/1000)^2;
for s = 1:length(r)
    subplot(4,3,s); hold on
    GoodFits{s} = DLPData.GlitchFlag{b(s)} == 0 & DLPData.StdResNorm{b(s)}<=0.2...
        & DLPData.Isat_m{b(s)}>0 & DLPData.Isat_m{b(s)}>0 & DLPData.Te{b(s)}<=20;
    f = sgolay_t(DLPData.Isat_m{b(s)}(GoodFits{s}),3,7);
%     f = DLPData.Isat_m{b(s)}(GoodFits{s});
    t_f = DLPData.time{b(s)}(GoodFits{s});
    plot(t_f,f)
    ylim([0,0.4])
    xlim([4,5])
    title(['R: ',num2str(r(b(s))),' shot: ',num2str(DLPData.shot(b(s)))])
    
    SteadyStateRng = find(t_f>4.5 & t_f<4.6);
    Is(b(s)) = mean(f(SteadyStateRng));
    dIs(b(s))= std(f(SteadyStateRng),1);
    tIs(b(s))= mean(t_f(SteadyStateRng));
    plot(t_f(SteadyStateRng),f(SteadyStateRng),'r')

    plot(tIs(b(s)),Is(b(s)),'ro')
end

% Integrate ion flux over the entire cross sectional area
r_F = r(b)-6.0;
F = Is(b)/(e_c*A_fp);
dF = dIs(b)/(e_c*A_fp);
np = r_F>=0;
Sp = 2*pi*trapz(r_F(np)/100,F(np).*r_F(np)/100)
nn = r_F<=0;
Sn = 2*pi*trapz(-r_F(nn)/100,F(nn).*r_F(nn)/100)
S  = mean([Sp,Sn])
dS = std([Sp,Sn],1)

figure; 
errorbar(r_F,F,dF,'ko-')

title(['S_+ = ',num2str(Sp,2),' [s^{-1}]'])
xlim([-5,5])
xlabel('R [cm]')
ylim([0,15e23])
ylabel('\Gamma [m^{-2}s^{-1}]')
grid on
set(findobj('-property','YTick'),'box','on')
set(findobj('-property','NextPlot'),'color','w')
set(gcf,'Position',[360  350  570  265])

%% 
FitShow = 0;
FitPreview
