% Fetch and load data from MPEX server:
% =========================================================================

clear all
close all

% Color coding:
% B_NORM -> Yellow -> S1
% B_TANG -> Red    -> S2

fetchDataFromServer = 0;

if fetchDataFromServer
    % Radial scan:
    % =====================================================================
    shotlist =  [30000 + [148,149,150,151,153,154,156,157,159,160,161,162]];
    r        =           [1.5,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5] ;
    arr      =           ['T','D','D','D','D','D','D','D','D','D','D','D'] ;

    % Attenuators:
    % =====================================================================
    REF_IN_ATT = 0;
    AD8302_ATT = 20;

    % Splitter:
    % =====================================================================
    REF_SPLT_ATT = 3.52;
    BDOT_SPLT_ATT = 3.52;
    
    % Address for the bdot probe data:
    % =====================================================================
    % B_NORM -> S1:
    DataAddress{1} = ['\MPEX::TOP.MACHOPS1:BDOT_3']; % VMAG_1
    DataAddress{2} = ['\MPEX::TOP.MACHOPS1:BDOT_4']; % VPH00_1
    DataAddress{3} = ['\MPEX::TOP.MACHOPS1:BDOT_5']; % VPH90_1

    % B_TANG -> S2:
    DataAddress{4} = ['\MPEX::TOP.MACHOPS1:BDOT_6']; % VMAG_2
    DataAddress{5} = ['\MPEX::TOP.MACHOPS1:BDOT_7']; % VPH00_2
    DataAddress{6} = ['\MPEX::TOP.MACHOPS1:BDOT_8']; % VPH90_2

    % Reference -> REF_IN:
    DataAddress{7} = ['\MPEX::TOP.MACHOPS1:BDOT_9']; % VREF_0
    
    % DLP 12.5 ion sat
    DataAddress{8} = ['\MPEX::TOP.MACHOPS1:LP_1'];

    % Fwd RF power trace:
    DataAddress{9} = ['\MPEX::TOP.MACHOPS1:RF_FWD_PWR'];
    
    % Fwd RF power trace:
    DataAddress{10} = ['\MPEX::TOP.MACHOPS1:RF_REF_PWR'];
    
    % Load data from server:
    % =====================================================================
    for ii = 1:numel(shotlist)
            for ch = 1:10
                [f{ch},t{ch}] = my_mdsvalue_v2(shotlist(ii),DataAddress(ch));   
            end

            % Assign raw data to variables:
            % S1:
            VMAG_1{ii}  = f{1}{1};
            VPH00_1{ii} = f{2}{1};
            VPH90_1{ii} = f{3}{1};

            % S2:
            VMAG_2{ii}  = f{4}{1};
            VPH00_2{ii} = f{5}{1};
            VPH90_2{ii} = f{6}{1};

            % REF_IN:
            VREF_0{ii}  = f{7}{1};

            % time:
            t_V{ii} = t{1}{1}(1:end-1);
            
            % Ne:
            isat{ii}   = f{8}{1};
            t_isat{ii} = t{8}{1}(1:end-1);
            
            % RF_FWD:
            FWD{ii}   = f{9}{1};
            t_FWD{ii} = t{9}{1}(1:end-1);
            
            % RF_REF:
            REF{ii}   = f{10}{1};
            t_REF{ii} = t{10}{1}(1:end-1);
    end
    
    varList = {'shotlist','r','arr',...
               'REF_IN_ATT','AD8302_ATT',...
               'REF_SPLT_ATT','BDOT_SPLT_ATT',...
               'DataAddress',...
               'VMAG_1','VPH00_1','VPH90_1',...
               'VMAG_2','VPH00_2','VPH90_2',...
               'VREF_0','t_V',...
               'isat','t_isat',...
               'FWD','t_FWD',...
               'REF','t_REF'};

           
    save('Step_1_GetData_BdotProbe_2020_06_09.mat',varList{:})
else
    load('Step_1_GetData_BdotProbe_2020_06_09.mat')
end

% Plasma radius:
rp = 7.5 - (r+1.5);

% Extract absolute value of Reference signal:
% =========================================================================
% Attenuation ratios:
att1 = 10^(-REF_IN_ATT/20);
att2 = 10^(-REF_SPLT_ATT/20);
att3 = 10^(-55/20);
att4 = att1*att2*att3;

% Reference signal voltage at port S0 (REF_IN):
a0 = 0.0055;
b0 = 0.7645;
for ii = 1:numel(shotlist)
    S0{ii} = ((VREF_0{ii}-a0)/b0);
end

% Reference signal voltage at matching box DC output port:
for ii = 1:numel(shotlist)
    S0_MatchingBox{ii} = S0{ii}/(att1*att2);
end

% Compute amplitude ratio R1 and R2:
% =========================================================================
for ii = 1:numel(shotlist)
    % Calibration factors:
    a = 1.23;
    b = 0.61;
    
    % S1:
    A = (VMAG_1{ii} - a)/b;
    R1{ii} = (10.^A); % R1 = S1/(S0*att5) where att5 = 10^(-AD8302_ATT/20)
        
    % S2:
    A = (VMAG_2{ii} - a)/b;
    R2{ii} = (10.^A); % R2 = S2/(S0*att) where att = 10^(-AD8302_ATT/20)
end

% Extract absolute values of S1 and S2:
% =========================================================================
att5 = 10^(-AD8302_ATT/20);
att6 = 10^(-BDOT_SPLT_ATT/20);

for ii = 1:numel(shotlist)
    S1{ii} = R1{ii}.*S0{ii}*att5;
    S2{ii} = R2{ii}.*S0{ii}*att5;
    
    % Bdot probe output voltages:
    V_BDOT_NORM{ii} = S1{ii}/att6;
    V_BDOT_TANG{ii} = S2{ii}/att6;
end

% Compute phase:
% =========================================================================
for ii = 1:numel(S1) 
    % Loop for each element of V00:
    N = numel(VPH00_1{ii});   
    for jj = 1:N
       % Determine quadrant:
       if VPH90_1{ii}(jj)<1
           % Upper half
           phase_1{ii}(jj) = (VPH00_1{ii}(jj) - (-0.0248))/(+0.0108);
       else
           % Lower half
           phase_1{ii}(jj) = (VPH00_1{ii}(jj) - (+3.9   ))/(-0.0109);
       end       
    end
end

for ii = 1:numel(S2) 
    % Loop for each element of V00:
    N = numel(VPH00_2{ii});   
    for jj = 1:N
       % Determine quadrant:
       if VPH90_2{ii}(jj)<1
           % Upper half
           phase_2{ii}(jj) = (VPH00_2{ii}(jj) - (-0.0248))/(+0.0108);
       else
           % Lower half
           phase_2{ii}(jj) = (VPH00_2{ii}(jj) - (+3.9   ))/(-0.0109);
       end       
    end
end

%%
% Plot data:
% =========================================================================
close all

figure
hold on
for ii = 1:numel(S1)
    h(1) = plot(t_V{ii},S0{ii})
    h(2) = plot(t_V{ii},S1{ii})
    h(3) = plot(t_V{ii},S2{ii})
end
legend(h,'S0','S1','S2')
xlim([4.1,4.7])

figure
subplot(2,1,1)
hold on
for ii = 1:numel(S1);
    h(1) = plot(t_V{ii},(phase_1{ii}*1/180))
    h(2) = plot(t_V{ii},(phase_2{ii}*1/180))
end
grid on
xlim([4.1,4.7])
legend(h,'S1','S2')

subplot(2,1,2)
hold on
for ii = 1:numel(S1);
    h(1) = plot(t_V{ii},unwrap(phase_1{ii}*pi/180)/pi)
    h(2) = plot(t_V{ii},unwrap(phase_2{ii}*pi/180)/pi)
end
grid on
xlim([4.1,4.7])
legend(h,'S1','S2')


figure
hold on
for ii = 1:numel(S1);
    h = plot3(t_V{ii},rp(ii)*ones(size(S1{ii})),(phase_2{ii}*1/180))
    h.LineStyle = 'none';
    h.Marker = '.';
end
xlim([4.2,4.4])
view([20,30])

figure
hold on
for ii = 1:numel(S1);
    h = plot3(t_V{ii},rp(ii)*ones(size(S1{ii})),S2{ii})
    h.LineStyle = 'none';
    h.Marker = '.';
end
xlim([4.2,4.4])
view([20,30])


