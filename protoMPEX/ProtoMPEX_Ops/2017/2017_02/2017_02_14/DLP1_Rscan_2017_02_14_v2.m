% this script was created April 26th 2017
% We are plotting the radial scan performed at DLP 1.5

% The measurements here reported, were performed Feb 14th 2017. According
% to the operator's logbook, the probe at spool 1.5 was Mach-DLP #1 and the
% probe at spool 4.5 was Mach-DLP #2.

clear all 
close all

CMPT = 0;

if CMPT == 1
% load timing table
D = importdata('DLP1_table_2017_02_14.xlsx');
T = D.data; % table of values
N = length(T(:,1));
% col 1: shot number
% col 2: radial position
% col 3: tStart
% col 4: tEnd
% col 5: t1
% col 6: t2
x  = T(:,2);
t1 = T(:,5);
t2 = T(:,6);
shotlist = 13200 + T(:,1);

% Data address:
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
% Use isolated voltage sweep source for DLP 1.5
DataAddress{1} = [RootAddress,'INT_4MM_1']; % V
DataAddress{2} = [RootAddress,'INT_4MM_2']; % I
Config.V_Att = 1;  % Output voltage of DLP box (Voltage) = V_Att*Digitized data 
Config.I_Att = 1;  % Output voltage of DLP box (Current) = I_att*Digitized data
Config.V_cal = [(0.46e-3)^-1,0]; % Voltage output of DLP = V_cal(1)*Output voltage of DLP box + V_cal(2) %  2.1739e+03
Config.I_cal = [-1,0]; % Current output of DLP = (I_cal(2) + Output voltage of DLP box)/Ical(1)

% Configuration variables
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
DLP = 1.5;
% Tip length
% On July 11th, we removed the DLP-Mach #2 from spool 1.5 and was
% inspected. The measured probe tip length was between 1.5 and 1.7 mm with an
% uncertainty of +- 0.2 mm
%Config.L_tip = 2/1000;
Config.L_tip = 1.8/1000;
Config.D_tip = 0.254/1000; % [m]
Config.FitFunction = 2; 
Config.AreaType = 1; % 1: Cylindrical + cap

for s = 1:N
    Config.tStart = T(s,3);
    Config.tEnd = T(s,4);
    [a,b,c,Ifit,Ip,Vp,tm,Vsweep,Isweep] = DLP_fit_V5_2017_02_14(Config,shotlist(s),DataAddress);
    Ni{s} = 0.5*(a{1}{1} + a{1}{2});
    Ni_u{s} = a{1}{1};
    Ni_l{s} = a{1}{2};
    Te{s} = b{1};
    time{s} = c{1};
end
    save('DLP1_2017_02_14.mat')
else
    load('DLP1_2017_02_14.mat')
end

%%
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
% Preview each shot
figure; 
for s = 1:N
    subplot(5,4,s);
    hold on
    rng = find(time{s}>=t1(s) & time{s}<=t2(s));
    plot(time{s},Ni{s},C{s})
    plot(time{s}(rng),Ni{s}(rng),C{s},'LineWidth',2)    
    ylim([0,2e20])
    xlim([4.15,4.35])
    title(['r = ',num2str(x(s)),', shot = ',num2str(shotlist(s))])
    
    ne_m(s) = mean(Ni{s}(rng));
    dne_m(s) = std(Ni{s}(rng),1,2);

end
figure; 
for s = 1:N
    subplot(5,4,s);
    hold on
    rng = find(time{s}>=t1(s) & time{s}<=t2(s));
    plot(time{s},Te{s},C{s})
    plot(time{s}(rng),Te{s}(rng),C{s},'LineWidth',2)    
    ylim([0,10])
    xlim([4.15,4.35])
    title(['r = ',num2str(x(s)),', shot = ',num2str(shotlist(s))])
    
    Te_m(s) = mean(Te{s}(rng));
    dTe_m(s) = std(Te{s}(rng),1,2);

end

% For this experiment, we used DLP 1.5 and we used the Mach probe # 1 as
% described in page 70 and 80 of Journal #2 (2016-2017).
% Based on this, we need to move 8.5 to 8.75 cm into the chamber from the fully
% retracted position to reach the center of the vacuum chamber
xoffset = 8.75;

[f,fn] = sort(x);
figure; 
subplot(2,1,1);
errorbar(x(fn)-xoffset,ne_m(fn),dne_m(fn),'ko')
title('n_e on DLP 1.5')
ylim([0,1e20])
xlim([-5,5])
subplot(2,1,2);
errorbar(x(fn)-xoffset,Te_m(fn),dTe_m(fn),'ko')
title('T_e on DLP 1.5')
ylim([0,4])
xlim([-5,5])
set(gcf,'color','w')
box on

% Create table to export to EXCEL
if 1
G = [shotlist(fn),x(fn)-xoffset,ne_m(fn)',dne_m(fn)',Te_m(fn)',dTe_m(fn)'];
F = {'Shot','R [cm]','n_e [m^-3]','dNe [m^-3]','T_e [eV]','dTe [eV]'};
FileName = 'NeTe_Spool_1_5_2017_02_14.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end