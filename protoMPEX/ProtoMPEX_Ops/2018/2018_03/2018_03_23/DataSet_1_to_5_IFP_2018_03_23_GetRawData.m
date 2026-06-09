clear all
close all

DownloadData = 0;
SaveData = 0;
DataName = 'DataSet_1_to_5_IFP_2018_03_23_RawData';

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
    MetaData.shot = [20000 +   [381:384,386:389,399:411,458:466]];
%     MetaData.r =               [0 ] ;
    MetaData.DataSetDescription     = ['DateSet 1 to 5: IFP XP with high density'];
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


%% Preview data
close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
% #########################################################################
% Dataset 5
% It looks like we may have not increased TR2 in this shot sequence
% #########################################################################
Shotlist = [20400 + [59 ,60 ,61 ,62 ,64 ]];
TR2       =         [120,140,160,180,100] ;

RawData5 = GetRawDataSubset(Shotlist,RawData);

figure; 
subplot(3,2,1);hold on
for s = 1:length(TR2)
    h(s) = plot(RawData5.t_I{s}(1:end-1),RawData5.I{s},C{s});
end
ylim([-0.4,0.4])
xlim([4.5,4.51])

subplot(3,2,2);hold on
for s = 1:length(TR2)
    h(s) = plot(RawData5.t_TR2{s}(1:end-1),sgolay_t(RawData5.TR2{s},3,89)*1000,C{s});
end
legend(h,num2str(RawData5.shot'))
ylim([0,200])

subplot(3,2,3);hold on
for s = 1:length(TR2)
    rng = find(RawData5.t_PG1{s}(1:end-1)>4.1);
    h(s) = plot(RawData5.t_PG1{s}(rng),(RawData5.PG1{s}(rng)-mean(RawData5.PG1{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData5.shot'))
xlim([4,5])
ylim([-0.01,0.3])

subplot(3,2,4);hold on
for s = 1:length(TR2)
    rng = find(RawData5.t_PG2{s}(1:end-1)>3.95);
    h(s) = plot(RawData5.t_PG2{s}(rng),(RawData5.PG2{s}(rng)-mean(RawData5.PG2{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData5.shot'))
xlim([4,5])
ylim([-0.01,1])

subplot(3,2,5);hold on
for s = 1:length(TR2)
    rng = find(RawData5.t_PG3{s}(1:end-1)>3.85);
    h(s) = plot(RawData5.t_PG3{s}(rng),(RawData5.PG3{s}(rng)-1*mean(RawData5.PG3{s}(rng(1:10))))*2/7.5,C{s});
end
% legend(h,num2str(RawData5.shot'))
xlim([4,5])
ylim([-0.01,0.3])

set(findobj('-property','YTick'),'box','on')
set(gcf,'color','w')
