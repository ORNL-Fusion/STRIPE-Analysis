% Create DLPinfo table automatically by reading the operators logbook
% stored in the MDSplus tree

close all
clear all

% The following script searches over the entire "Data Analysis" folder
% structure and collects all the XPS performed 

RetrieveXP_List = 0;
% 1: It takes 34 seconds to retrieve all the data
% 0: It takes 0.2 seconds to load previoulsy retrieved data

if RetrieveXP_List
RootAddress = 'C:\Users\nfc\Documents\Proto-MPEX Data Analysis';
a = dir(RootAddress);
x = 0;
p = 0;
for sx = 3:length(a);
    if strcmp(a(sx).name(1:3),'201')
        x = x + 1;
        Year{x} = a(sx).name;
        b = dir([RootAddress,'\',Year{x}]);
        y = 0;
        for sy = 3:length(b);
           if strcmp(b(sy).name(1:3),'201')
              y = y + 1;
              Month{x}{y} = b(sy).name;
                c = dir([RootAddress,'\',Year{x},'\',Month{x}{y}]);
                z = 0;
                for sz = 3:length(c);
                    if strcmp(c(sz).name(1:3),'201')
                       if  length(c(sz).name)>10
                           continue
                       end
                       p = p + 1;
                       z = z + 1;
                       Day{x}{y}{z} = c(sz).name;
                       AddressXYZ{x}{y}{z} = [RootAddress,'\',Year{x},'\',Month{x}{y},'\',Day{x}{y}{z}];
                       AddressP{p} = AddressXYZ{x}{y}{z};
                       Date{p} = Day{x}{y}{z};
                       Topic{p} = cell2mat(importdata([AddressXYZ{x}{y}{z},'\Topic_',Day{x}{y}{z},'.txt'],'',3));
                       if isstruct(Topic{p}) | isempty(Topic{p})
                           p
                       end
                    end
                end
           end
        end   
     end
end
    save('ProtoMPEX_XP_list','AddressP', 'AddressXYZ','Date','Topic','RootAddress')
else
    load('ProtoMPEX_XP_list')
end

%% Gather all DLP relevant shots
% From 2016_12_19 to 2017_07_20

RetrieveOperShotSummary = 0;
% 1: To retrieve data from spreadsheets takes 27 seconds
% 0: To load from previously retrieved data: 0.5 seconds


if RetrieveOperShotSummary
    StartDate = '2016_12_19';
    Start_nxp = find(strcmp(Date,StartDate)==1);
for nxp = Start_nxp:length(Date)
    
    try
    % #####################################################################
    % Retrieve Operator's logbook
    % #####################################################################
    [num,txt,OperShotSummary{nxp}] = xlsread([AddressP{nxp},'\',Date{nxp},'_OperShotSummary']);
    
    %======================================================================
    % Retrieve the shotlist for the current date
    %======================================================================
    % Shot #, find location of column
    [~,c] = find(strcmp(OperShotSummary{nxp},'Shot #')==1);
    ShotList_All{nxp} = (num(:,c(1)));
    %======================================================================
    % Collect the shots that can be DLP analyzed
    % Only select shots that have: FP=1, HMJ=1, Y~=OUT
    %======================================================================
    % Shots that we need for DLP analysis are those which:
    % 1- Are full pulses FP
    % 2 - Have the probe inserted in the plasma
    % 3 - Are high density discharges

    % FP, Full pulse
    [r,c] = find(strcmp(OperShotSummary{nxp},'FP')==1);
    ShotList_FP = ShotList_All{nxp}(num(:,c(1)) == 1);
    % Y_cm. When Probe is out, we generally write "OUT" in this column and
    % also on the "Gauge" column. a number on either of each of this column
    % generally indicates that the probe is in
                if nxp == 85;
                    nxp
                end
    [r,c] = find(strcmp(OperShotSummary{nxp},'Y [cm]')==1);
%     ShotList_ProbeIn = ShotList_All{nxp}( (~isnan(num(:,c(1))) | ~isnan(num(:,c(1)+1))) );
    ShotList_ProbeIn = ShotList_All{nxp}( ~isnan(num(:,c(1))) );
    % HMJ, Helicon mode jump. It really refers to a helicon plasma
    [r,c] = find(strcmp(OperShotSummary{nxp},'HMJ')==1);
    ShotList_HMJ = ShotList_All{nxp}(num(:,c(1)) == 1);
    % Sift through the shot list and collect those which are DLP suitable:
    [ShotList_FP_ProbeIn] = intersect(ShotList_FP,ShotList_ProbeIn);
    [ShotList_FP_ProneIn_HMJ] = intersect(ShotList_HMJ,ShotList_FP_ProbeIn);
    [ShotDLP{nxp},ia{nxp},ib{nxp}] = intersect(ShotList_All{nxp},ShotList_FP_ProneIn_HMJ);
    N_DLP(nxp) = length(ShotDLP{nxp});
    % nxp is the index for each XP, there are 113 XPs
    % For each XP there are shots which can DLP analyzed: ShotDLP{nxp}
    % where ShotDLP{nxp} = ShotList_All{nxp}( ia{nxp} )
    
    %======================================================================
    % Collect information from the DLP suitable shots:
    %======================================================================
    % Comments-------------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'Note')==1);
    Comment{nxp} =  OperShotSummary{nxp}(r+ia{nxp},c);
    % MDS+ V and I channels------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'MDS+')==1);
    MDS_I{nxp} =  upper(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)));
    MDS_V{nxp} =  upper(OperShotSummary{nxp}(r(1)+ia{nxp},c(2))); % use  cell2mat(MDS_V{nxp}(1)) to convert to char
    % Spool # -------------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'Spool #')==1);
    if isempty(r)
        [r,c] = find(strcmp(OperShotSummary{nxp},'Spool')==1);
    end        
    Spool{nxp} =  cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1))); 
    % Probe L and D -------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},' L [mm]')==1);
    Ltip{nxp} = cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)));   % In mm
    Dtip{nxp} = cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)+1)); % In mm
    % Probe Drive ---------------------------------------------------------
    % Elijah, Elijah (non-isolated), isol.
    [r,c] = find(strcmp(OperShotSummary{nxp},'Drive')==1);
    DriveTypeRaw{nxp} =  OperShotSummary{nxp}(r(1)+ia{nxp},c(1));
    for s = 1:length(DriveTypeRaw{nxp})     
%             if s == 37 && nxp == 88;s,end
        if strcmp(DriveTypeRaw{nxp}{s}(1:4),'Elij') 
            DriveType{nxp}{s} = 'niso';
        elseif strcmp(DriveTypeRaw{nxp}{s}(1:4),'isol')
            DriveType{nxp}{s} = 'iso';
        end     
    end
    % Vatt and Iatt -------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'V')==1);
    Vatt{nxp} =  cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)));
    Iatt{nxp} =  cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)+1));
    % Y_cm ----------------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'Y [cm]')==1);
    Y{nxp} =  OperShotSummary{nxp}(r(1)+ia{nxp},c(1));
    % ICH_FP --------------------------------------------------------------
    [r0,cICH] = find(strcmp(OperShotSummary{nxp},'ICH')==1);
    [r,c] = find(strcmp(OperShotSummary{nxp},'FP')==1);
    ICH_FP{nxp} =  OperShotSummary{nxp}(r(3)+ia{nxp},c(c>cICH));
    % ECH_FP --------------------------------------------------------------
    [r0,cECH] = find(strcmp(OperShotSummary{nxp},'ECH (28 GHz)')==1);
    [r,c] = find(strcmp(OperShotSummary{nxp},'FP')==1);
    ECH_FP{nxp} =  OperShotSummary{nxp}(r(3)+ia{nxp},c(c>cECH & c<cICH));
    % HMJ  ----------------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'HMJ')==1);
    HMJ{nxp} =  cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)));
    % RF pulse length -----------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'T [ms]')==1);
    if isempty(r) 
        RfPulseLength{nxp} =  150*ones(size(ia{nxp}));
    else
        RfPulseLength{nxp} =  cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)));
    end
    % Gauge ---------------------------------------------------------------
    [r,c] = find(strcmp(OperShotSummary{nxp},'Gauge')==1);
    Gauge{nxp} =  cell2mat(OperShotSummary{nxp}(r(1)+ia{nxp},c(1)));
    
    catch
        continue
    end
end
    save('DLPinfo_Raw_2016_12_19_to_2017_07_20','ShotList_All','ShotDLP',...
        'ia','N_DLP','Comment','MDS_I','MDS_V','Spool','Ltip','Dtip',...
        'DriveTypeRaw','DriveType','Vatt','Iatt','Y','ICH_FP','ECH_FP',...
        'Start_nxp','StartDate','HMJ','RfPulseLength','OperShotSummary','Gauge')
else
    load('DLPinfo_Raw_2016_12_19_to_2017_07_20')
end

figure
plot(N_DLP)

%% Test put together a DLPinfo spreadsheet

Dheader = {'Year','Month','Day','Shot','Spool','AxisType','R_cm','Gauge',...
    'HMJ','Tstart','Tend','Gstart','Gend','Vcal_1','Vcal_2','Vatt',...
    'Ical_1','Ical_2','Iatt','SweepType','Ltip_mm','Dtip_mm','MDS_V',...
    'MDS_I','ICH_FP','ECH_FP'};	

tstart = 4.18; % DLP analysis start time
RFtstart = 4.15; % RF pulse start time
for nxp = Start_nxp:length(ShotDLP)
    D{nxp} = Dheader;

    for ns = 1:length(ShotDLP{nxp})
        
     if strcmp(DriveType{nxp}{ns},'niso')
         vcal1 = 12.05;
         vcal2 = 0.205;
         ical1 = -142.5;
         ical2 = 0.015;
     elseif strcmp(DriveType{nxp}{ns},'iso')
         vcal1 = (0.46e-3)^-1;
         vcal2 = 0;
         ical1 = -1;
         ical2 = 0; 
     end
     if strcmpi(ICH_FP{nxp}{ns},'x') || isnan(ICH_FP{nxp}{ns})
         ichfp = 0;
     elseif strcmpi(ICH_FP{nxp}{ns},'y') || strcmp(ICH_FP{nxp}{ns},'1') || ICH_FP{nxp}{ns} == 1
         ichfp = 1;
     end   
     if strcmpi(ECH_FP{nxp}{ns},'x') || isnan(ECH_FP{nxp}{ns})
         echfp = 0;
     elseif strcmpi(ECH_FP{nxp}{ns},'y') || strcmp(ECH_FP{nxp}{ns},'1') || ECH_FP{nxp}{ns} == 1
         echfp = 1;
     end
     
     tend = tstart + (RfPulseLength{nxp}(ns)*1e-3) - 10e-3;
     Gstart = tend - 40e-3;
     Gend   = tend - 10e-3;
     
         D{nxp}(ns+1,:) = {Date{nxp}(1:4),Date{nxp}(6:7),Date{nxp}(9:10),...
         ShotDLP{nxp}(ns),Spool{nxp}(ns),[],cell2mat(Y{nxp}(ns)),...
         Gauge{nxp}(ns),HMJ{nxp}(ns),tstart,tend,Gstart,Gend,vcal1,vcal2,Vatt{nxp}(ns),...
         ical1,ical2,Iatt{nxp}(ns),DriveType{nxp}{ns},Ltip{nxp}(ns),...
         Dtip{nxp}(ns),MDS_V{nxp}{ns},MDS_I{nxp}{ns},ichfp,...
         echfp};     
    end
end

T = D{Start_nxp}(1,:);
for nxp = Start_nxp:length(ShotDLP)
    T = [T;D{nxp}(2:end,:)];
end

SaveDLPinfoTable = 1;

if SaveDLPinfoTable
DLPinfoName = ['DLPinfo_',num2str(Date{Start_nxp}),'_to_',num2str(Date{end}),'.xlsx']; 
writetable(table(T),DLPinfoName,'WriteVariableNames',0)
end

%% General notes:
% tasks remaining:
% none

% =========================================================================
% Useful code lines:
% =========================================================================

% 1 - Find keyword, use the following:-------------------------------------
% for nxp = 58:113;
%    k = strfind(Topic{nxp},'Long');
%    if ~isempty(k
%        nxp % output the XP number
%     end
% end

% 2 - Vector operations on Cells: -----------------------------------------
% Use the following to compute a function to all the elements of a cell
% cellfun('size',Topic,2); where 'size' can be replaced by a function @f
% anothe example: multiply each cell in ShotDLP by 3 using the anonynous
% function provided @(x) 3*x
% cellfun(@(x) 3*x, ShotDLP, 'UniformOutput', false)

% Plot Excel tables on command window: ------------------------------------
% use cell2table(D{nxp}) to plot spreadsheets

% =========================================================================
% Changes needed in OperShotSummary
% =========================================================================

% At the end of each operational day, we need to run a code like this one
% and compute and make sure all the required data is available

% 2017_02_14: nxp = 70,  DLP 6.5 is used on-axis and sweep, the Y column has an "x"
% this means that the script in this file does not include these shots as
% DLP relevant

% 2017_02_14: nxp = 70,  Data has no time trace, DLP code fails

% 2017-06-09: nxp = 99,  RF pulse length in these shots is reflected as "X"

% Use binary drop down lists in cells to restrict inputs to 0 1 or X 
% - Use sheets to produce the drop lists for all binary data types

% =========================================================================
% Additions to PreviewDLP code:
% =========================================================================
% - Include the shot comments and the topic of that day's experiment
