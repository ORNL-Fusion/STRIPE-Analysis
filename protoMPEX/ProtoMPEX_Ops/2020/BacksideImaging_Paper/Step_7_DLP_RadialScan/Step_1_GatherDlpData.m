% Step 1: Gather data raw DLP data from radial profile during 2nd harmonic EBW heating
% experiments:

clc 
clear all
close all

fetchDataFromServer = 0;

if fetchDataFromServer
    % Define shotlist and associated metadata:
    % =========================================================================
    MetaData.shot      = 29400 + [52  ,62  ,63  ,64  ,65 ,66 ,67 ,68 ,69 ,70 ,71  ];
    MetaData.rgauge    =         [-3.0,-2.0,-2.5,-1.0,0.0,1.0,0.5,2.0,3.0,3.5,-3.5];
    MetaData.DataSetDescription = ['DLP 12.5 radial scan during 2nd harmonic EBW heating XP, 2020_02_27'];
    MetaData.DLP_tipLength      = 0.002;
    MetaData.DLP_tipDiam        = 0.000254;
    MetaData.DLP_tipUnits       = 'meters';
    MetaData.rgaugeUnits        = 'cm';
    
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

    fileName = ['Step_1_GetDlpRawData.mat'];
    save(fileName,'RawData')
else
    fileName = ['Step_1_GetDlpRawData.mat'];
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

