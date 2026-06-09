close all
clear all

shotlist = 24000 + [678];
Rad      =         [0]  ;

% shotlist = 24000 + [697,706];
% Rad      =         [0  ,0  ]  ;

% Acquiring Isat current from MP tip 1 and tip 2
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

DataAddress{1} = [RootAddress,'MP_I1']; % I_senseA Channel 09
DataAddress{2} = [RootAddress,'MP_I2']; % I_senseA Channel 10
DataAddress{3} = [RootAddress,'LP_1']; % I
DataAddress{4} = [RootAddress,'LP_V_RAMP']; % V
% DataAddress{5} = [RootAddress,'PG1']; % V
DataAddress{6} = [RootAddress,'ICH_FWD_PWR']; % V

% Mach Probe tip lengths

MP_type = 'MP1' % MP1 or MP2

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
% [PG,tpg] = my_mdsvalue_v2(shotlist,DataAddress(5));
[RF,trf] = my_mdsvalue_v2(shotlist,DataAddress(6));


% Calibration using 10 Ohm resitor
R = 10.0; % Ohms

    
%% Process data:
C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};
tStart = 4.1; % [s]
tEnd = 4.66;

figure;
for s = 1:length(shotlist)
    
    rng = find(tf1{s}>=tStart & tf2{s}<=tEnd);
%     Select data range:
    MPI_1{s} = f_1{s}(rng); %sgolay_t(f_1{s}(rng),3,41);
    MPI_2{s} = f_2{s}(rng);
    t1{s} = tf1{1}(rng);
    t2{s} = tf2{1}(rng);
    
    % Clean data fron NaNs
    MPI_1{s}(isnan(MPI_1{s})) = 0;
    MPI_2{s}(isnan(MPI_2{s})) = 0;
    % Clear data from glitches
    [n1,~] = find(abs(MPI_1{s})>5);MPI_1{s}(n1) = 0;
    [n2,~] = find(abs(MPI_2{s})>5);MPI_2{s}(n2) = 0;
        
    hold on
    plot(t1{s},-1*MPI_1{s},'k')
    plot(t2{s},-1*MPI_2{s},'r')
    ylim([-1,5])
%     title(['R= ',num2str(x(s))])
    
    % Choose probe that collects the most current (facing flow)
    S1 = trapz(abs(MPI_1{s}));
    S2 = trapz(abs(MPI_2{s}));
    dMPI_1{s} = diff(-MPI_1{s});
    dMPI_2{s} = diff(-MPI_2{s});
       
    noffset = 4;
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
            % Points must be separated by at least 5 ms
            T = 5e-3;
            dt = t2{s}(2)-t1{s}(1);
            dn = floor(0.8*T/dt);
            [r,~] = peakseek(h1,dn,0.1); % zero crossing, minpeadist, minpkht
            
            % Select based on heigh requirement, this section needs work
%             MinPeak = 0.8*mean(h1(r));
            tStart_MinPeak = find(t2{s}>4.32);
            MinPeak = 0.05*mean(h2(r(find(r>6e3)) + noffset));
            rng = find(h2(r+noffset)>=MinPeak); % min height requirement
            
            RisingEdgeMaxLocs{s} = r(rng)+noffset; % location of Rising edge
            RisingEdgeMinLocs{s} = r(rng)-noffset; % Base of Rising edge
            
            % Falling edge:
%             [f,~] = peakseek(-h1,dn,0.1); % zero crossing
%             rng = find(h2(f-noffset)>=MinPeak); % min height requirement
%             FallingEdgeMaxLocs{s} = f(rng)-noffset; % location of falling edge
%             FallingEdgeMinLocs{s} = f(rng)-noffset; % Base of Falling edge
    
            FallingEdgeMaxLocs{s} = RisingEdgeMaxLocs{s} + 14; % location of falling edge
            FallingEdgeMinLocs{s} = RisingEdgeMaxLocs{s} + 22; % Base of Falling edge
            
    fr = 5;
    Y1{s} = sgolay_t(-MPI_1{s},3,fr);     
    J1_min{s} = Y1{s}(RisingEdgeMinLocs{s}); % define base current
    J1_max{s} = Y1{s}(RisingEdgeMaxLocs{s}); % peak current
    J1{s} = J1_max{s}-J1_min{s};
    
    Y2{s} = sgolay_t(-MPI_2{s},3,fr);  
    J2_min{s} = Y2{s}(RisingEdgeMinLocs{s}); % define base current
    J2_max{s} = Y2{s}(RisingEdgeMaxLocs{s}); % peak current
    J2{s} = J2_max{s}-J2_min{s};

    plot(t2{s}(RisingEdgeMaxLocs{s}),J1_max{s},'ko')
    plot(t2{s}(RisingEdgeMinLocs{s}),J1_min{s},'ksq')
    
    plot(t2{s}(RisingEdgeMaxLocs{s}),J2_max{s},'ro')
    plot(t2{s}(RisingEdgeMinLocs{s}),J2_min{s},'rsq')

    plot(t2{s}(FallingEdgeMaxLocs{s}),-MPI_1{s}(FallingEdgeMaxLocs{s})  ,'gsq')
    plot(t2{s}(FallingEdgeMaxLocs{s}),-MPI_2{s}(FallingEdgeMaxLocs{s})  ,'csq')

end

%%
% close all
% figure; 
% for s = 1:length(shotlist)
%     subplot(1,length(shotlist),s);hold on
%     hold on
%     plot(t1{s}(RisingEdgeMaxLocs{s}),J1{s}*1000/R,'ko')
%     plot(t2{s}(RisingEdgeMaxLocs{s}),J2{s}*1000/R,'ro')
%     title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
%     ylim([0,200])
% end

try
    clear MachNumberSS dMachNumberSS
%     close all
catch
end

figure; 
for s = 1:length(shotlist)

    subplot(4,4,s);hold on
    MachNumber{s} = log((L_mp1/L_mp2)*J1{s}./J2{s})/1.66;
    plot(tisat{s}(1:end-1),(Isat{s})*3)
    plot(trf{s}(1:end-1),(RF{s})*0.3)
    plot(t1{s}(RisingEdgeMaxLocs{s}),MachNumber{s},'ko')
    ylim([-0.5,1.5])
    xlim([4,5])
    ylabel('Mach number')
    title(['R= ',num2str(Rad(s)-8.5),', ',num2str(shotlist(s))])
    
    % Gather steady state MachNumber
    tStart_Mach = 4.57;
    tEnd_Mach   = 4.60;
    
%     tStart_Mach = 4.4;
%     tEnd_Mach   = 4.47;
    
    rngSS = find(RisingEdgeMaxLocs{s}>=find(t2{s}>tStart_Mach,1) & RisingEdgeMaxLocs{s}<=find(t2{s}>tEnd_Mach,1));
    MachNumberSS(s) = mean(MachNumber{s}(rngSS));
    dMachNumberSS(s) = std(MachNumber{s}(rngSS));
        if ~isreal(MachNumberSS(s))
                MachNumberSS(s) = NaN;
                dMachNumberSS(s) = NaN;
        end
    
    errorbar(4.57,MachNumberSS(s),dMachNumberSS(s),'gsq')
        if sum(Rad(s) == [9.0,6.0,10,9.5,9.25])>=1
                MachNumberSS(s) = NaN;
                dMachNumberSS(s) = NaN;
        continue
    end
end


figure;
subplot(2,1,1)
[a,b] = sort(Rad);
errorbar(a-1*8.5,MachNumberSS(b),1*dMachNumberSS(b),'r-','LineWidth',2)
xlim([-4,4])
ylim([0,1])
grid on
set(gcf,'color','w')
box on
t = title('Mach number'); t.Interpreter = 'latex'; t.FontSize = 13;
t = xlabel('$PS1$ $[kA]$ '); t.Interpreter = 'latex'; t.FontSize = 11;

%% Convert data into Excel
% Save data to Excel spreadsheet:
% Create table to export to EXCEL
[a,b] = sort(Rad);
if 0
G = [shotlist(b)',Rad(b)',round(MachNumberSS(b),2,'significant')',round(dMachNumberSS(b),2,'significant')'];
F = {'Shot','R [cm]','Mach Number','dM'};
FileName = 'Mach_Spool_2_5_2017_11_21.xlsx';
xlswrite(FileName,[F;num2cell(G)]);
end