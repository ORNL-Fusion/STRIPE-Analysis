% Step 1: Gather data raw DLP data from ECH power scan with 120 kW helicon
% plasma

clc 
clear all
close all

fetchDataFromServer = 0;

if fetchDataFromServer
    % Define shotlist and associated metadata:
    % =========================================================================   
    MetaData.shot      = [29671,29669,29665,29668,29661,29663,29662,29673,29655,29657];
    MetaData.pwr28GHz  = [5.92 ,26.12,36.84,38.29,45.41,50.98,52.61,54.07,70.07,70.48];
    MetaData.DataSetDescription = ['DLP 12.5, hot spot during 2nd harmonic EBW heating power scan XP, 2020_03_30'];
    MetaData.DLP_tipLength      = 0.002;
    MetaData.DLP_tipDiam        = 0.000254;
    MetaData.DLP_tipUnits       = 'meters';
    MetaData.pwr28GHzUnits      = 'kW';
    
    % -------------------------------------------------------------------------
    % 3 - Define the MDSplus channels to use and associated addresses
    % -------------------------------------------------------------------------
    % First, we load the list of addresses
    ProtoMPEX_AddressInfo_2018_06_04

    % DLP 12.5:
    Address.V       = {aV{1}};
    Address.I       = {aI{3}};

    % RF fwd power:
    Address.RF_13MHz = {aRF_fwd{1}};
    
    % 28 GHz:
    Address.RF_28GHz = {aECH{1}};

    % -------------------------------------------------------------------------
    % 4 - Make a table for that contains the shots and their addresses
    % -------------------------------------------------------------------------
    % To make the address table we use the following function: MakeAddressTable.m
    % INPUTS: ShotList (array), ChannelAddress (cell)
    % OUTPUT: AddressTable; the table is plotted on the command line
    AddressTable = MakeAddressTable(MetaData.shot,Address);

    % -------------------------------------------------------------------------
    % 5 - Download data from the server using the address table
    % -------------------------------------------------------------------------
    % This step requres the function: GetRawDataFromServer.m
    % INPUTS: AddressTable, MetaData
    % OUTPUT: RawData; a description of the raw data is plotted on the command
    % line
    [RawData] = GetRawDataFromServer(AddressTable,MetaData);
    RawData

    fileName = ['Step_1_GetDlpRawData_EchPwrScan_20200330.mat'];
    save(fileName,'RawData')
else
    fileName = ['Step_1_GetDlpRawData_EchPwrScan_20200330.mat'];
    load(fileName);
end

%%
% Preview data:
figure('color','w');
hold on
for ss = 1:numel(RawData.shot)
    plot(RawData.t_RF_13MHz{ss}(1:end-1),RawData.RF_13MHz{ss})
end
box on

figure('color','w');
hold on
for ss = 1:numel(RawData.shot)
    plot(RawData.t_RF_28GHz{ss}(1:end-1),RawData.RF_28GHz{ss})
end
box on

figure('color','w');
hold on
for ss = 1:2
    Z = RawData.I{ss};
    plot(RawData.t_I{ss}(1:end-1),Z)
end
xlim([4.1,4.8])
box on

