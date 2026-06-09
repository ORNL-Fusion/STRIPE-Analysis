clear all
close all

DownloadData = 0;
SaveData = 0;
DataName = 'PS2_Scan_2018_09_04_RawData';

switch DownloadData
    case 1
% To download data from the server, we use the following functions:
% - MakeAddressTable.m
% - GetRawDataFromServer.m
% =========================================================================
% DEFINE ADDRESSES
% =========================================================================
ProtoMPEX_AddressInfo_2018_09_06
% =========================================================================
% COMPOSE THE ADDRESS TABLE WITH SHOTS
% =========================================================================
    %--------------------------------------------------------------------------------------------------------------------------------------
    % Define SHOTLIST and associated METADATA
    MetaData.shot         =    [22507,22504,22503,22493,22511];
    MetaData.PS2          =    [2600 ,2900 ,3200 ,3500 ,3800 ];
    MetaData.DataSetDescription     = ['PS2 scan using 28 GHz ECH experiments '];
    %--------------------------------------------------------------------------------------------------------------------------------------
    % Define ADDRESSES
    Address.V   =        {aV{1}};
    Address.I   =        {aI{1}};
    Address.PG1 =        {aPG{1}};
    Address.PG2 =        {aPG{2}};
    Address.PG3 =        {aPG{3}};
    Address.PG4 =        {aPG{4}};
    Address.RF_fwd =     {aRF_fwd{1}};
    Address.ECH =        {aECH{1}};   
    Address.t_zero =     {aT0{1}};
    Address.TR2 =        {aTR2{1}};
    Address.PS1 =        {aPS1{1}};
    Address.PS2 =        {aPS2{1}};    
% =========================================================================
% MAKE ADDRESS TABLE
% =========================================================================
   % INPUTS: shotlist (array), Address {cell}
   % OUTPUT: AddressTable
    AddressTable = MakeAddressTable(MetaData.shot,Address);
   
% =========================================================================
% GATHER DATA FROM SERVER GIVEN THE ADDRESS TABLE
% =========================================================================
    % Required inputs: AddressTable, MetaData
    % Outputs: RawData
    [RawData] = GetRawDataFromServer(AddressTable,MetaData);

    struct2table(AddressTable)
    RawData

if SaveData == 1
save(DataName,'AddressTable','RawData')
end

case 0
    % In this case, we do not need to get the data from the server. we
    % need to reload the data from the local memory
    load(DataName)
    struct2table(AddressTable)
    RawData
end
% #########################################################################
% #########################################################################
% End of Gather Raw Data
% #########################################################################
% #########################################################################

%% Calculate DLP data

AttType = 'Vx1,Ix1'; % Attenuation on the digitized signals
SweepType = 'iso'; % "iso" for isolated sweep, "niso" for non-isolated sweep
DLPType = '12MP';
Config.tStart = 4.15; % [s]
Config.tEnd = 4.65;
Config.FitFunction = 2; 
Config.Center_V = 1; % Remove offset on V: 1 (yes) 0(no)
Config.Center_I = 1; % Remove offset on I: 1 (yes) 0(no)
Config.FilterDataInput = 1; % Filter input data with savitsky Golay filter order 3 frame size 7
Config.SGF = 7;
Config.TimeMode = 2; % Effective time of sweep, (1) start of ramp or (2) mean time of ramp
Config.AMU = 2; % Ion mass in AMU
Config.AreaType = 1; % 1: Cylindrical + cap
ProtoMPEX_DLPInfo_2018_09_06

DLPData = DLP_fit_V6(Config,RawData);

%% Assign DLP data to variables
ni = DLPData.Ni;
Te = DLPData.Te;
time = DLPData.time;
Ifit = DLPData.Ifit;
Ip = DLPData.Ip;
Vp = DLPData.Vp;
tm = DLPData.tm;
Vsweep = DLPData.Vsweep;
Isweep = DLPData.Isweep;
GlitchFlag = DLPData.GlitchFlag;
StdResNorm = DLPData.StdResNorm;

for s = 1:length(RawData.shot)
    Ni{s} = 0.5*(ni{s}{1} + ni{s}{2});
end
t_ech = RawData.t_ECH;
ECH = RawData.ECH;
RF = RawData.RF_fwd;
t_rf = RawData.t_RF_fwd;
t = RawData.t_zero;

%% Filter DLP data
% Details of the data filtering
ApplyFilter.Status = 1;
ApplyFilter.Order = 3;
ApplyFilter.Frame = 11;

% Collect shots that have enought points to be filtered and produce
% filteres data sets
k = 0;
for s = 1:length(RawData.shot)
    GoodFits{s} = GlitchFlag{s} == 0 & StdResNorm{s}<=0.15 & Ni{s}>0 & Ni{s}<1e21 & Te{s}<=19;
    GoodShots_Boolean(s) = length(Ni{s}(GoodFits{s}))>ApplyFilter.Frame;
    Ni{s} = Ni{s}(GoodFits{s});
    Te{s} = Te{s}(GoodFits{s});
    time{s} = time{s}(GoodFits{s});
    if GoodShots_Boolean(s) == 1
            time_Filtered{s} = time{s};
            Ni_Filtered{s} = sgolay_t(Ni{s},ApplyFilter.Order,ApplyFilter.Frame);
            Te_Filtered{s} = sgolay_t(Te{s},ApplyFilter.Order,ApplyFilter.Frame);
    end
end

%% Plot DLP data
TimePlotStart = 4.15;
TimePlotEnd = 4.7;
NiMaxLim = 50e18; 
TeMaxLim = 15; 
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

figure; 
subplot(2,2,1); hold on
for s = find(GoodShots_Boolean==1)
          
    switch ApplyFilter.Status
        case 0
            hNe(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
        case 1
            hNe(s) = plot(time_Filtered{s},Ni_Filtered{s},C{s},'lineWidth',2);
    end
    hold on
    plot(t_ech{s}(1:length(ECH{s})),ECH{s}*1e19,C{s},'lineWidth',0.5)
    plot(t_rf{s}(1:length(RF{s})),(RF{s}.^2)*1e19,C{s},'lineWidth',0.5)

end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(hNe(1),['DLP ',num2str(DLP)],'location','SouthEast')
ylim([0,NiMaxLim])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,2); hold on
for s = find(GoodShots_Boolean==1)
    
   switch ApplyFilter.Status
        case 0
            hTe(s) = plot(time{s},Te{s},C{s},'lineWidth',2);
        case 1
            hTe(s) = plot(time_Filtered{s},Te_Filtered{s},C{s},'lineWidth',2);
    end
    L{s} = [num2str(RawData.shot(s)),' ,t=',num2str(t{s}(10:14))];

end
legend(hTe(GoodShots_Boolean),L(GoodShots_Boolean),'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,TeMaxLim])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,3); hold on
for s = 1:length(RawData.shot)
    h(s) = plot(time{s},Ni{s},C{s},'lineWidth',2);
end
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
legend(h,['DLP ',num2str(DLP)],'location','SouthEast')
ylim([0,NiMaxLim])
xlim([TimePlotStart,TimePlotEnd])

subplot(2,2,4); hold on
for s = 1:length(RawData.shot)
    h(s) = plot(time{s},Te{s},C{s},'lineWidth',2);
    L{s} = [num2str(RawData.shot(s)),' ,t=',num2str(t{s}(10:14))];
end
legend(h,L,'location','NorthEast')
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,TeMaxLim])
xlim([TimePlotStart,TimePlotEnd])

set(findobj('-Property','YTick'),'box','on')
set(gcf,'color','w')

%% Plot steady-state part:
close all
% ECH on part : 4.485 to 4.5
% ECH off part: 4.53 to 4.55
PS2 = RawData.MetaData.PS2;

% Get steady state part for each shot
for s = 1:length(Te)
    % ECH OFF
    rng = find(time{s}>=4.53 & time{s}<=4.6);
    ne0(s) = mean(Ni{s}(rng));
    dne0(s) = std(Ni{s}(rng),1);
    Te0(s) = mean(Te{s}(rng));
    dTe0(s) = std(Te{s}(rng),1);
    
    % ECH ON
    rng = find(time{s}>=4.48 & time{s}<=4.5);
    ne1(s) = mean(Ni{s}(rng));
    dne1(s) = std(Ni{s}(rng),1);
    Te1(s) = mean(Te{s}(rng));
    dTe1(s) = std(Te{s}(rng),1);
    
end

figure; set(gcf,'color','w')
% Ne ----------------------------------------------------------------------
ax(1) = subplot(2,2,1); box on;  hold on; grid on
h(1) = errorbar(PS2,ne0,dne0)
h(1).Color = 'k'

h(2) = errorbar(PS2,ne1,dne1)
h(2).Color = 'r'

ylabel(ax(1),'$n_e$ $[m^{-3}]$','Interpreter','Latex','FontSize',13)
ylim([0,4e19])

% Te ----------------------------------------------------------------------
ax(2) = subplot(2,2,2); box on; hold on; grid on
h(3) = errorbar(PS2,Te0,dTe0)
h(3).Color = 'k'

h(4) = errorbar(PS2,Te1,dTe1)
h(4).Color = 'r'
ylim([0,7])
ylabel(ax(2),'$T_e$ $[eV]$','Interpreter','Latex','FontSize',13)
legend([h(3),h(4)],'28 GHz OFF','28 GHz ON')

% Pe ----------------------------------------------------------------------
ax(3) = subplot(2,2,3); box on; hold on; grid on
h(5) = errorbar(PS2,e_c*ne0.*Te0,e_c*dne0.*dTe0)
h(5).Color = 'k'

h(6) = errorbar(PS2,e_c*ne1.*Te1,e_c*dne1.*dTe1)
h(6).Color = 'r'
ylabel(ax(3),'$P_e$ $[Pa]$','Interpreter','Latex','FontSize',13)
ylim([0,15])

for s = 1:3
    xlim(ax(s),[2500,4000]);
    xlabel(ax(s),'PS2 [I]','Interpreter','Latex','FontName','Times')
end

for s = 1:6
    h(s).LineWidth = 1;
end


%% Convert DLP data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
[~,b] = sort(RawData.shot);

% Input data
for s = 1:length(RawData.shot)
    
DataToWrite.FileName =  ['DLP_Shot_',num2str(RawData.shot(s))];
DataToWrite.Column{1}.Heading = 'time [s]';
DataToWrite.Column{1}.Data    = time{s};
DataToWrite.Column{2}.Heading = 'ne un-filtered [m^-3]';
DataToWrite.Column{2}.Data    = Ni{s};
DataToWrite.Column{3}.Heading = 'Te un-filtered [eV]';
DataToWrite.Column{3}.Data    = Te{s};
DataToWrite.Column{4}.Heading = 'ne filtered [m^-3]';
DataToWrite.Column{4}.Data    = Ni_Filtered{s};
DataToWrite.Column{5}.Heading = 'Te filtered [eV]';
DataToWrite.Column{5}.Data    = Te_Filtered{s};

if 0
% Write data to an Excel file
DataToExcel(DataToWrite)
end
end
