% ==================
% This code is designed to retrive data from Mach Probe, plot current vs time
% from MP_n1 and MP_n2 and calculate the current ratio
% Isense A goes to MP tip 1 and Isense B goes to MP tip 2
% First written on Feb 14, 2017

% MP1_1 points towards Target, MP1_2 towards dump
% Unless otherwise stated, red is for probe facing dump and black to probe
% facing target
% ==================

close all
clear all

shotlist = 13000 + [186,187,253]; % Feb 14th 2017, MP 4.5 on-axis 
x = [0,0,0];

% Acquiring Isat current from MP tip 1 and tip 2
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

DataAddress{1} = [RootAddress,'INT_2MM_1']; % I_senseA Channel 09
DataAddress{2} = [RootAddress,'INT_2MM_2']; % I_senseA Channel 10
DataAddress{3} = [RootAddress,'TARGET_LP']; % I
DataAddress{4} = [RootAddress,'LP_V_RAMP']; % V

% Mach Probe tip lengths

MP_type = 'MP1'; % MP1 or MP2'

switch MP_type
    case 'MP1'
        L_mp1 = 5.5; % [mm]
        L_mp2 = 5.5; %[mm]
    case 'MP2'
        L_mp1 = 4.0; % [mm]
        L_mp2 = 5.0; %[mm]
end

% Standard probe tip configuration
% --------------------------------
% MPI_1  = probe tip facing target 
% MPI_2 = probe tip facing dump

[f_1,tf1]   = my_mdsvalue_v2(shotlist,DataAddress(1)); % [V] signal from digitizer from Tip A
[f_2,tf2]   = my_mdsvalue_v2(shotlist,DataAddress(2)); % [V] signal from digitizer from Tip B
[Isat,tisat] = my_mdsvalue_v2(shotlist,DataAddress(3));
[V,tv] = my_mdsvalue_v2(shotlist,DataAddress(4));


% Calibration using 10 Ohm resitor
R = 10.0; % Ohms

%%
%close all
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
tStart = 4.18; % [s]
tEnd = 4.36;
    
% Plot currents
for s = 1:length(shotlist)
%     subplot(5,4,s);
    hold on
    rng = 1:35e3; % No time trace on Feb 14th, use index instead
    
    % Select data range:
    MPI_1{s} = f_1{s}(rng); % to smooth use: sgolay_t(f_1{s}(rng),3,41);
    MPI_2{s} = f_2{s}(rng);
    % Clean data fron NaNs
    MPI_1{s}(isnan(MPI_1{s})) = 0;
    MPI_2{s}(isnan(MPI_2{s})) = 0;
    % Clear data from glitches, any voltage greater than 5 V is considered
    % a glitch
    [n1,~] = find(abs(MPI_1{s})>5);MPI_1{s}(n1) = 0;
    [n2,~] = find(abs(MPI_2{s})>5);MPI_2{s}(n2) = 0;
    
    plot(-1*MPI_1{s}*1000/R,'k')
    plot(-1*MPI_2{s}*1000/R,'r')
%     plot(abs(Isat{s})*1e3,'g:')
%     xlabel('t[s]')
    ylabel('[mA]')
    ylim([-50,400])
    xlim([1,4]*1e4)  
    
    title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
end

%%
% Process data:
figure;
for s = 1:length(shotlist)
    
    % Clean data fron NaNs
    MPI_1{s}(isnan(MPI_1{s})) = 0;
    MPI_2{s}(isnan(MPI_2{s})) = 0;
    % Clear data from glitches
    [n1,~] = find(abs(MPI_1{s})>5);MPI_1{s}(n1) = 0;
    [n2,~] = find(abs(MPI_2{s})>5);MPI_2{s}(n2) = 0;
        
    hold on
    plot(-1*MPI_1{s},'k')
    plot(-1*MPI_2{s},'r')
    ylim([-1,5])
    xlim([1,4]*1e4)
        title(['R= ',num2str(x(s))])
    
    % Choose probe that collects the most current (facing flow)
    % Integrate the data: largest integral == largest flow
    S1 = trapz(abs(MPI_1{s})); 
    S2 = trapz(abs(MPI_2{s}));
    
    % rising and faling edges produce delta functions to the peakseek.m
    % function
    dMPI_1{s} = diff(-MPI_1{s});
    dMPI_2{s} = diff(-MPI_2{s});
       
    noffset = 8; 
    if S1 > S2
            h1 = dMPI_1{s};
            h2 = -MPI_1{s};
            MaxI = 1; % probe 1 collects most current
    else
            h1 = dMPI_2{s};
            h2 = -MPI_2{s};
            MaxI = 2; % probe 2 collects most current
    end
            % Rising edge:
            [r,~] = peakseek(h1,200,0.1); % function to find zero crossing, minpeadist, minpkht
            
            % Select peaks based on heigh requirement, this section needs work
            MinPeak = 0.5*mean(h2(r(find(r>25e3)) + noffset));
            rng = find(h2(r+noffset)>=MinPeak); % min height requirement
            
            RisingEdgeMaxLocs{s} = r(rng)+noffset; % location of Rising edge
            RisingEdgeMinLocs{s} = r(rng)-noffset; % Base of Rising edge
            
            % Falling edge:
            [f,~] = peakseek(-h1,200,0.1); % zero crossing
            rng = find(h2(f-noffset)>=MinPeak); % min height requirement
            FallingEdgeMaxLocs{s} = f(rng)-noffset; % location of falling edge
            FallingEdgeMinLocs{s} = f(rng)-noffset; % Base of Falling edge
    
 
    J1_min{s} = -MPI_1{s}(RisingEdgeMinLocs{s}); % define base current
    J1_max{s} = -MPI_1{s}(RisingEdgeMaxLocs{s}); % peak current
    J1{s} = J1_max{s}-J1_min{s};
    
    J2_min{s} = -MPI_2{s}(RisingEdgeMinLocs{s}); % define base current
    J2_max{s} = -MPI_2{s}(RisingEdgeMaxLocs{s}); % peak current
    J2{s} = J2_max{s}-J2_min{s};

    plot(RisingEdgeMaxLocs{s},J1_max{s},'ko')
    plot(RisingEdgeMinLocs{s},J1_min{s},'k.')
    
    plot(RisingEdgeMaxLocs{s},J2_max{s},'ro')
    plot(RisingEdgeMinLocs{s},J2_min{s},'r.')

    plot(FallingEdgeMaxLocs{s},-MPI_1{s}(FallingEdgeMaxLocs{s})  ,'ksq')
    plot(FallingEdgeMaxLocs{s},-MPI_2{s}(FallingEdgeMaxLocs{s})  ,'rsq')

end

%%
% close all
figure; 
for s = 1:3
    subplot(2,2,s); 
    hold on
    plot(RisingEdgeMaxLocs{s},J1{s}*1000/R,'ko')
    plot(RisingEdgeMaxLocs{s},J2{s}*1000/R,'ro')
    title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
    ylim([0,200])
end

figure; 
for s = 1:3
    subplot(2,2,s);hold on
    if shotlist(s) == 13186
            MachNumber{s} = real(log(J2{s}./J1{s}))/1.66;
    else
            MachNumber{s} = log(J1{s}./J2{s})/1.66;
    end
    plot(RisingEdgeMaxLocs{s},MachNumber{s},'ko')
    plot(abs(Isat{s})*0.3)
    ylim([0,1])
    xlim([1,4]*1e4)
    ylabel('Mach number')
    title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
    
    % Gather steady state MachNumber
    rngSS = find(RisingEdgeMaxLocs{s}>=25e3 & RisingEdgeMaxLocs{s}<=31e3);
    MachNumberSS(s) = mean(MachNumber{s}(rngSS));
    dMachNumberSS(s) = std(MachNumber{s}(rngSS));
        if ~isreal(MachNumberSS(s))
                MachNumberSS(s) = NaN;
                dMachNumberSS(s) = NaN;
        end
    
    errorbar(mean([25,31]*1e3),MachNumberSS(s),dMachNumberSS(s),'gsq')
end

% Create table to export to EXCEL
if 0
G = [shotlist',x',MachNumberSS',dMachNumberSS'];
F = {'Shot','R [cm]','MachNum','dMachNum'};
FileName = 'MachNum_Spool_4_5_2017_02_14.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end
