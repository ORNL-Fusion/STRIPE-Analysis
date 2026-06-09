% Mach probe data:
% ECH power scan with downhill configuration

close all
clear all

% Define shotlist:
shotlist = [29749, 29745, 29746,29744,29747,29742,29748,29741];
Rad      = [1    ,2     ,3     ,4    ,5    ,6    ,7    ,8    ];

% Mach Probe tip lengths:
MP_type = 'MP1';
switch MP_type
    case 'MP1'
        L_mp1 = 5.5; % [mm]
        L_mp2 = 5.5; %[mm]
    case 'MP2'
        L_mp1 = 4.0; % [mm]
        L_mp2 = 5.0; %[mm]
end

% Calibration using 10 Ohm resitor:
R = 10.0;

% Get data:
% MPI_1  = probe tip facing target 
% MPI_2  = probe tip facing dump
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RootAddress = [Stem,Branch];
DataAddress{1} = [RootAddress,'MP_I1']; % I_senseA Channel 09
DataAddress{2} = [RootAddress,'MP_I2']; % I_senseA Channel 10
DataAddress{3} = [RootAddress,'LP_1']; % I
DataAddress{4} = [RootAddress,'PWR_28GHZ']; % Isx

[f_1,tf1]    = my_mdsvalue_v2(shotlist,DataAddress(1)); % [V] signal from digitizer from Tip A
[f_2,tf2]    = my_mdsvalue_v2(shotlist,DataAddress(2)); % [V] signal from digitizer from Tip B
[Isat,tisat] = my_mdsvalue_v2(shotlist,DataAddress(3));
[Iech,ech]   = my_mdsvalue_v2(shotlist,DataAddress(4));

C = {'k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:','k','r','bl','g','m','k:','r:','bl:','g:','m:'};

%% figure out ech power level:
figure
hold on
for s = 1:length(shotlist)
    hdum1(s) = plot(ech{s}(1:end-1),Iech{s},C{s},'LineWidth',2);
end
xlim([4,5])
legend(hdum1,num2str(shotlist'))
    
%% Process data:
tStart = 4.1; % [s]
tEnd = 4.66;

figure;
for s = 1:length(shotlist)
    
    % Select data range:
    rng = find(tf1{s}>=tStart & tf2{s}<=tEnd);
    MPI_1{s} = sgolay_t(f_1{s}(rng),3,21);
    MPI_2{s} = sgolay_t(f_2{s}(rng),3,21);
    t1{s} = tf1{1}(rng);
    t2{s} = tf2{1}(rng);
    
    % Clean data fron NaNs:
    MPI_1{s}(isnan(MPI_1{s})) = 0;
    MPI_2{s}(isnan(MPI_2{s})) = 0;
    
    % Clear data from glitches:
    [n1,~] = find(abs(MPI_1{s})>5);MPI_1{s}(n1) = 0;
    [n2,~] = find(abs(MPI_2{s})>5);MPI_2{s}(n2) = 0;
        
    hold on
    plot(t1{s},-1*MPI_1{s},'k')
    plot(t2{s},-1*MPI_2{s},'r')
    ylim([-1,5])
    
    % Choose probe that collects the most current (facing flow):
    S1 = trapz(abs(MPI_1{s}));
    S2 = trapz(abs(MPI_2{s}));
    dMPI_1{s} = diff(-MPI_1{s});
    dMPI_2{s} = diff(-MPI_2{s});
       
    noffset = 15;
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
            [r,~] = peakseek(h1,dn,0.0341); % zero crossing, minpeadist, minpkht
            
            % Select based on heigh requirement, this section needs work
            MinPeak = 4*mean(h1(r));
            tStart_MinPeak = find(t2{s}>4.32);
            %MinPeak = 0.05*mean(h2(r(find(r>6e3)) + noffset));
            rng = find(h2(r+noffset)>=MinPeak); % min height requirement
            
            RisingEdgeMaxLocs{s} = r(rng)+noffset; % location of Rising edge
            RisingEdgeMinLocs{s} = r(rng)-noffset; % Base of Rising edge
            
            % Falling edge:
            [f,~] = peakseek(-h1,dn,0.1); % zero crossing
            rng = find(h2(f-noffset)>=MinPeak); % min height requirement
            FallingEdgeMaxLocs{s} = f(rng)-noffset; % location of falling edge
            FallingEdgeMinLocs{s} = f(rng)-noffset; % Base of Falling edge
    
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

    plot(t2{s}(FallingEdgeMaxLocs{s}),-MPI_1{s}(FallingEdgeMaxLocs{s})  ,'ksq')
    plot(t2{s}(FallingEdgeMaxLocs{s}),-MPI_2{s}(FallingEdgeMaxLocs{s})  ,'rsq')

end

%%
% close all
% figure; 
% for s = 1:length(shotlist)
%     %subplot(1,length(shotlist),s);hold on
%     hold on
%     plot(t1{s}(RisingEdgeMaxLocs{s}),J1{s}*1000/R,'ko')
%     plot(t2{s}(RisingEdgeMaxLocs{s}),J2{s}*1000/R,'ro')
%     title(['R= ',num2str(x(s)),', ',num2str(shotlist(s))])
%     ylim([0,200])
% end

figure('color','w'); 
subplot(2,1,1)
[~,b] = sort(Rad);
for s = 1:length(shotlist)
    hold on
    MachNumber{s} = -log((L_mp1/L_mp2)*J1{s}./J2{s})/1.66;
    %plot(tisat{s}(1:end-1),(Isat{s})*0.3)
    hdum2(s) = plot(t1{s}(RisingEdgeMaxLocs{s}),MachNumber{s},C{b(s)},'LineWidth',2);
    ylim([0,0.6])
    xlim([4.2,4.7])
    box on
    ylabel('Mach number')
    
    % Gather steady state MachNumber
    rngSS = find(RisingEdgeMaxLocs{s}>=find(t2{s}>4.55,1) & RisingEdgeMaxLocs{s}<=find(t2{s}>4.56,1));
    plot(t1{s}(RisingEdgeMaxLocs{s}(rngSS)),MachNumber{s}(rngSS),C{b(s)},'LineWidth',6)
    MachNumberSS(s) = mean(MachNumber{s}(rngSS));
    dMachNumberSS(s) = std(MachNumber{s}(rngSS));
        if ~isreal(MachNumberSS(s))
                MachNumberSS(s) = NaN;
                dMachNumberSS(s) = NaN;
        end
end
s = 1;
plot(t1{s}(RisingEdgeMaxLocs{s}),MachNumber{s},C{b(s)},'LineWidth',3);
ldum = legend(hdum2,num2str(shotlist'));
ldum.Location = 'northwest';

% figure;
subplot(2,1,2)
[~,b] = sort(Rad);
errorbar(Rad(b),MachNumberSS(b),1*dMachNumberSS(b),'r-','LineWidth',2)
xlim([0,10])
ylim([0,1])
grid on
set(gcf,'color','w')
box on
t = title('Mach number'); t.Interpreter = 'latex'; t.FontSize = 13;
t = xlabel('ECH pwr [A.U.]'); t.Interpreter = 'latex'; t.FontSize = 11;

figureName = ['machNumber_EchPwrScan'];
saveas(gcf,figureName,'tiffn')

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