clear all
close all

DownloadData = 1;
SaveData = 1;
DataName = 'DataSet_1_to_8_IFP_2018_03_27_RawData';

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
    MetaData.shot = [20500 +   [8:9, 10,11,14:16, 17:26, 32,33, 36,37, 38:44,46:48, 49:56, 126:136] ];
%     MetaData.r    = 
    MetaData.DataSetDescription     = ['DataSet_1_to_8_IFP_2018_03_27'];
    %--------------------------------------------------------------------------------------------------------------------------------------
    % Define ADDRESSES
    Address.V   =        {aV{1}};
    Address.I   =        {aI{2}};
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

% return

%% Preview data
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

Shotlist_6 = 20500 + [36:37]; r_6 = [1,2];
% Shotlist_6 = 20500 + [126:136];
% Shotlist_6 = 20500 + [38:44,46:48]; r_6 = [1,1.5,2, 2.5, 3, 0, -0.5, -1, -1.5, -2];

RawData_6 = GetRawDataSubset(Shotlist_6,RawData);

figure; 
subplot(3,2,1);hold on
for s = 1:length(r_6)
    h(s) = plot(RawData_6.t_I{s}(1:end-1),RawData_6.I{s},C{s});
end
ylim([-0.4,0.4])
xlim([4.5,4.51])

subplot(3,2,2);hold on
for s = 1:length(r_6)
    h(s) = plot(RawData_6.t_PS1{s}(1:end-1),sgolay_t(RawData_6.PS1{s},3,89)*1000,C{s});
end
legend(h,num2str(RawData_6.shot'))
ylim([0,5000])

subplot(3,2,3);hold on
for s = 1:length(r_6)
    rng = find(RawData_6.t_PG1{s}(1:end-1)>4.1);
    h(s) = plot(RawData_6.t_PG1{s}(rng),(RawData_6.PG1{s}(rng)-mean(RawData_6.PG1{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData_6.shot'))
xlim([4,5])
ylim([-0.01,0.3])

subplot(3,2,4);hold on
for s = 1:length(r_6)
    rng = find(RawData_6.t_PG2{s}(1:end-1)>3.95);
    h(s) = plot(RawData_6.t_PG2{s}(rng),(RawData_6.PG2{s}(rng)-mean(RawData_6.PG2{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData_6.shot'))
xlim([4,5])
ylim([-0.01,1])

subplot(3,2,5);hold on
for s = 1:length(r_6)
    rng = find(RawData_6.t_PG3{s}(1:end-1)>3.85);
    h(s) = plot(RawData_6.t_PG3{s}(rng),(RawData_6.PG3{s}(rng)-1*mean(RawData_6.PG3{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData_6.shot'))
xlim([4,5])
ylim([-0.01,0.3])

set(findobj('-property','YTick'),'box','on')
set(gcf,'color','w')
