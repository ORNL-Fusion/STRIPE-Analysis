% Load dataset 3 from 2018_04_13
clear all
close all

DataName = 'DataSet_1_to_6_IFP_2018_04_13_RawData';
load(DataName)
struct2table(AddressTable)
RawData

%% Data subset
Shotlist = 21000 + [248,249,251,252,253,254,255,256,257,258];
r        =         [6.0,6.5,7.0,7.5,8.0,8.5,5.5,5.0,4.5,4.0];

RawData_3 = GetRawDataSubset(Shotlist,RawData);

figure;
[~,b] = sort(r);

for s = 1:length(r)
    subplot(4,3,s)
    plot(RawData_3.t_I{b(s)}(1:end-1),RawData_3.I{b(s)})
    ylim([-0.2,0.2])
    xlim([4,5])
    title(['R: ',num2str(r(b(s)))])
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

DLPData = DLP_fit_V6(Config,RawData_3);

%%
close all
figure
A_fp = 0.25*sqrt(2)*pi*(1.22/1000)^2;
for s = 1:length(r)
    subplot(4,3,s); hold on
%     GoodFits{s} = DLPData.GlitchFlag{b(s)} == 0 & DLPData.StdResNorm{b(s)}<=0.7...
%         & DLPData.Isat_m{b(s)}>0 & DLPData.Isat_m{b(s)}>0 & DLPData.Te{b(s)}<=10;
    GoodFits{s} = DLPData.StdResNorm{b(s)}<=0.7...
        & DLPData.Isat_m{b(s)}>0 & DLPData.Isat_m{b(s)}>0 & DLPData.Te{b(s)}<=10;
    
    f = sgolay_t(DLPData.Isat_m{b(s)}(GoodFits{s}),3,7);
    t_f = DLPData.time{b(s)}(GoodFits{s});
    plot(t_f,f)
    ylim([0,0.2])
    xlim([4,5])
    title(['R: ',num2str(r(b(s))),' shot: ',num2str(DLPData.shot(b(s)))])
    
    SteadyStateRng = find(t_f>4.6 & t_f<4.62);
    Is(b(s)) = mean(f(SteadyStateRng));
    dIs(b(s))= std(f(SteadyStateRng),1);
    tIs(b(s))= mean(t_f(SteadyStateRng));
    plot(t_f(SteadyStateRng),f(SteadyStateRng),'r')

    plot(tIs(b(s)),Is(b(s)),'ro')
end

% Integrate ion flux over the entire cross sectional area
r_F = r(b)-6.25;
F = Is(b)/(e_c*A_fp);
dF = dIs(b)/(e_c*A_fp);
np = r_F>=0;
nn = r_F<=0;
Sp = 2*pi*trapz(r_F(np)/100,F(np).*r_F(np)/100)
Sn = 2*pi*trapz(-r_F(nn)/100,F(nn).*r_F(nn)/100)
S  = mean([Sp,Sn])
dS = std([Sp,Sn],1)

figure; 
errorbar(r_F,F,dF,'ko-')

title(['S_+ = ',num2str(S,3),' [s^{-1}]'])
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

%% Plot neutral gas data:
% close all
figure; 
D = RawData_6;

s = 1;
a = 0.13333;
t = D.t_PG1{s}(1:end-1);
rngX = [3.9,5];

p0(1) = mean(D.PG1{s}(find(t>3.8   & t<4.01)));
p{1}  = a*2*(D.PG1{s}-p0(1));
p0(2) = mean(D.PG2{s}(find(t>3.8 & t<3.9)));
p{2}  = a*2*(D.PG2{s}-p0(2));
p0(3) = mean(D.PG3{s}(find(t>3.8 & t<3.9)));
p{3}  = a*2*(D.PG3{s}-p0(3));
p0(4) = mean(D.PG4{s}(find(t>3.8 & t<3.9)));
p{4}  = a*10*(D.PG4{s}-p0(4));

subplot(2,2,1); hold on
h(1) = plot(t,p{1})
xlim(rngX)
ylim([-0.01,1])
title('P4, Target')
legend(h,['Shot: ' num2str(D.shot(s))])

subplot(2,2,2); hold on
plot(t,p{2})
xlim(rngX)
ylim([-0.01,1])
title('P1, Source')

subplot(2,2,3); hold on
plot(t,p{3})
% plot(t,p{1},'r')
xlim(rngX)
ylim([0,1])
title('P3, Heating section')

subplot(2,2,4); hold on
plot(t,p{4})
xlim(rngX)
ylim([0,1])
title('P2, Source')

set(findobj('-property','YTick'),'box','on')
set(findobj('-property','NextPlot'),'color','w')
set(findobj(gcf,'Type','line'),'LineWidth',3)


