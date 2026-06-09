clear all
close all

DownloadData = 0;
SaveData = 1;
DataName = 'RawData_GasFlowTimingScan_2018_05_01';

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
    MetaData.shot = 21000 + [742,743,744,745,830];
    MetaData.dt   =         [0  ,100,200,300,400];
    MetaData.DataSetDescription     = ['Systematically changed the gas timing relative to the RF pulse'];
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

