% Analyze the data associated with the radial plasma density measured
% during a MPEX-limiter configuration:

clc 
clear all
close all

fetchDataFromServer = 0;

if fetchDataFromServer
    % Define shotlist and associated metadata:
    % =========================================================================
    MetaData.shot      = [28308       ,28424       ,20908];
    MetaData.rgauge    = [10.0        ,10.0        ]; % [cm]
    MetaData.limMode   = {'Window'    ,'MPEX'      }
    MetaData.date      = ['2019_11_18','2019_11_26'];
    MetaData.rfPwr     = [129         ,126         ]; % [kW]
    MetaData.DataSetDescription     = ['DLP 12.5 on-axis shots scan for two magnetic configurations'];
    Metadata.PG_CalShot = 'shot 20908 is used to calibrate the PGs';
    % -------------------------------------------------------------------------
    % 3 - Define the MDSplus channels to use and associated addresses
    % -------------------------------------------------------------------------
    % First, we load the list of addresses
    ProtoMPEX_AddressInfo_2018_06_04

    % PS2:
    Address.PS2_I   = {aPS2_I{1}};
    % PS3:
    Address.PS3_V   = {aPS3_V{1}};
    % DLP 12.5:
    Address.V       = {aV{1}};
    Address.I       = {aI{3}};
    % Pressure gauges:
    Address.PG1     = {aPG{1}};
    Address.PG2     = {aPG{2}};
    Address.PG3     = {aPG{3}};
    Address.PG4     = {aPG{4}};
    % RF fwd power:
    Address.RF      = {aRF_fwd{1}};

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

    fileName = ['Step_9b_RawData.mat'];
    save(fileName,'RawData')
else
    fileName = ['Step_9b_RawData.mat'];
    load(fileName);
end

%%
% Preview data:
figure;
hold on
for ss = 1:numel(RawData.shot)
    plot(RawData.t_PS2_I{ss}(1:end-1),RawData.PS2_I{ss})
end

figure;
hold on
for ss = 1:2
    Z = RawData.I{ss};
    plot(RawData.t_I{ss}(1:end-1),Z)
end
xlim([4.1,4.8])