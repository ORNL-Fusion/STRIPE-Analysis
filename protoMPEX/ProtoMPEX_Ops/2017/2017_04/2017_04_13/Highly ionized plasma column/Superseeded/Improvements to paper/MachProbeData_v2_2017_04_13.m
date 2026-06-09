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

% shotlist = 13984; % April 13th 2017, MP 10.5 on-axis
% shotlist = [18556,18537]; % Dec 14th 2017, MP 10.5 on-axis
shotlist = 18537; % Dec 14th 2017, MP 10.5 on-axis
shotlist = [13984,18537]; % Dec 14th 2017, MP 10.5 on-axis
x = [0,0];

% shotlist = 18000 + [1033,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051]; 
% x      =         [4.5 ,4.5 ,4.2 ,3.9 ,3.6 ,3.3 ,3.0 ,2.7 ,2.4 ,2.1 ,1.8 ];

% Acquiring Isat current from MP tip 1 and tip 2
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];

DataAddress{1} = [RootAddress,'INT_2MM_1']; % I_senseA Channel 09
DataAddress{2} = [RootAddress,'INT_2MM_2']; % I_senseA Channel 10
DataAddress{3} = [RootAddress,'TARGET_LP']; % I
DataAddress{4} = [RootAddress,'LP_V_RAMP']; % V
DataAddress{5} = [RootAddress,'PG1']; % V
DataAddress{6} = [RootAddress,'RF_FWD_PWR']; % V

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
[PG,tpg] = my_mdsvalue_v2(shotlist,DataAddress(5));
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
    title(['R= ',num2str(x(s))])
    
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
            dn = floor(0.9*T/dt);
            [r,~] = peakseek(h1,dn,0.1); % zero crossing, minpeadist, minpkht
            
            % Select based on heigh requirement, this section needs work
%             MinPeak = 0.8*mean(h1(r));
            tStart_MinPeak = find(t2{s}>4.32);
            MinPeak = 0.2*mean(h2(r(find(r>25e3)) + noffset));
            rng = find(h2(r+noffset)>=MinPeak); % min height requirement
            
            RisingEdgeMaxLocs{s} = r(rng)+noffset; % location of Rising edge
            RisingEdgeMinLocs{s} = r(rng)-noffset; % Base of Rising edge
            
            % Falling edge:
            [f,~] = peakseek(-h1,dn,0.1); % zero crossing
            rng = find(h2(f-noffset)>=MinPeak); % min height requirement
            FallingEdgeMaxLocs{s} = f(rng)-noffset; % location of falling edge
            FallingEdgeMinLocs{s} = f(rng)-noffset; % Base of Falling edge
    
    fr = 11;
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

    plot(t2{s}(FallingEdgeMaxLocs{s}),-MPI_1{s}(FallingEdgeMaxLocs{s})  ,'ksq')
    plot(t2{s}(FallingEdgeMaxLocs{s}),-MPI_2{s}(FallingEdgeMaxLocs{s})  ,'rsq')

end

%%
% close all
figure; 
for s = 1:length(shotlist)
    subplot(1,length(shotlist),s);hold on
    hold on
    plot(t1{s}(RisingEdgeMaxLocs{s}),J1{s}*1000/R,'ko')
    plot(t2{s}(RisingEdgeMaxLocs{s}),J2{s}*1000/R,'ro')
    title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
    ylim([0,200])
end

figure; 
for s = 1:length(shotlist)
    subplot(1,length(shotlist),s);hold on
    if shotlist(s) == 13186
            MachNumber{s} = real(log(J2{s}./J1{s}))/1.66;
    else
            MachNumber{s} = log(J1{s}./J2{s})/1.66;
    end
    plot(t1{s}(RisingEdgeMaxLocs{s}),MachNumber{s},'ko')
    plot(tisat{s}(1:end-1),abs(Isat{s})*0.3)
    ylim([0,1.5])
%     xlim([1,4]*1e4)
    ylabel('Mach number')
    title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
    
    % Gather steady state MachNumber
    rngSS = find(RisingEdgeMaxLocs{s}>=25e3 & RisingEdgeMaxLocs{s}<=34e3);
    MachNumberSS(s) = mean(MachNumber{s}(rngSS));
    dMachNumberSS(s) = std(MachNumber{s}(rngSS));
        if ~isreal(MachNumberSS(s))
                MachNumberSS(s) = NaN;
                dMachNumberSS(s) = NaN;
        end
    
    errorbar(t1{s}(round((mean([25e3,34e3])))),MachNumberSS(s),dMachNumberSS(s),'gsq')
end

%% Figure for paper
figure
s = 1;
tEndData = 4.45;
rng_Data = find(t1{s}(RisingEdgeMaxLocs{s})<=tEndData);
plot(t1{s}(RisingEdgeMaxLocs{s}(rng_Data)),MachNumber{s}(rng_Data),'ko')
ylim([0,1.2])
xlim([4.1,4.5])
set(gcf,'position',[  687.6667  177.6667  325.3333  228.6667],'color','w')
box on
ylabel('Mach number','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')


figure
s = 2;
tEndData = 4.45;
rng_Data = find(t1{s}(RisingEdgeMaxLocs{s})<=tEndData);
plot(t1{s}(RisingEdgeMaxLocs{s}(rng_Data)),MachNumber{s}(rng_Data),'ko')
ylim([0,1])
xlim([4.1,4.5])
box on
ylabel('$M_{\parallel}$','Interpreter','latex','FontSize',14)
xlabel('$time$ [s]','Interpreter','latex')
set(gcf,'position',[  687.6667  177.6667  325.3333  228.6667],'color','w')
