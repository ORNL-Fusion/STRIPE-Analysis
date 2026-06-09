clear all
close all

DownloadData = 0;
SaveData = 0;
DataName = 'DataSet_Nischal_2018_04_13_RawData';

switch DownloadData
    case 1
% To download data from the server, we use the following functions:
% - MakeAddressTable.m
% - GetRawDataFromServer.m
% =========================================================================
% DEFINE ADDRESSES
% =========================================================================
ProtoMPEX_AddressInfo_2018_06_04
% =========================================================================
% COMPOSE THE ADDRESS TABLE WITH SHOTS
% =========================================================================
    %--------------------------------------------------------------------------------------------------------------------------------------
    % Define SHOTLIST and associated METADATA
    MetaData.shot = [21200 +   [60  ,65  ,66  ,67  ,68  ,69  ,70  ]];
    MetaData.r    =            [0   ,0   ,0   ,0   ,0   ,0   ,0   ] ;
    MetaData.PS1  =            [4.5 ,4.0 ,3.5 ,3.0 ,2.5 ,2.0 ,2.0 ]*1e3 ;
    MetaData.DataSetDescription     = ['IFP XP 2018_04_13: PS1 current scan'];
    %--------------------------------------------------------------------------------------------------------------------------------------
    % Define ADDRESSES
    Address.V   =        {aV{1}};
    Address.I   =        {aI{2}};
    Address.PG1 =        {aPG{1}};
    Address.PG2 =        {aPG{2}};
    Address.PG3 =        {aPG{3}};
    Address.PG4 =        {aPG{4}};
    Address.RF_fwd =     {aRF_fwd{1}};
%     Address.ECH =        {aECH{1}};   
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


%% Preview data
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

TR2 = RawData.TR2;
PS1 = RawData.PS1;

figure; 
subplot(3,2,1);hold on
for s = 1:length(PS1)
    h(s) = plot(RawData.t_I{s}(1:end-1),RawData.I{s},C{s});
end
ylim([-0.4,0.4])
xlim([4.5,4.51])

subplot(3,2,2);hold on
for s = 1:length(PS1)
    h(s) = plot(RawData.t_PS1{s}(1:end-1),sgolay_t(RawData.PS1{s},3,89)*1000,C{s});
end
legend(h,num2str(RawData.shot'))
ylim([0,5000])

subplot(3,2,3);hold on
for s = 1:length(PS1)
    rng = find(RawData.t_PG1{s}(1:end-1)>4.1);
    h(s) = plot(RawData.t_PG1{s}(rng),(RawData.PG1{s}(rng)-mean(RawData.PG1{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData.shot'))
xlim([4,5])
ylim([-0.01,0.3])

subplot(3,2,4);hold on
for s = 1:length(PS1)
    rng = find(RawData.t_PG2{s}(1:end-1)>3.95);
    h(s) = plot(RawData.t_PG2{s}(rng),(RawData.PG2{s}(rng)-mean(RawData.PG2{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData.shot'))
xlim([4,5])
ylim([-0.01,1])

subplot(3,2,5);hold on
for s = 1:length(PS1)
    rng = find(RawData.t_PG3{s}(1:end-1)>3.85);
    h(s) = plot(RawData.t_PG3{s}(rng),(RawData.PG3{s}(rng)-1*mean(RawData.PG3{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData.shot'))
xlim([4,5])
ylim([-0.01,0.3])

set(findobj('-property','YTick'),'box','on')
set(gcf,'color','w')
