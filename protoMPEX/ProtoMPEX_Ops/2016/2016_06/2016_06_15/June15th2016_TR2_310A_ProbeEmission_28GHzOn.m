% D2 discharge using gas puffing but RF power low enough to prevent the
% formation of the high density mode jump, the ECH (28 GHz) is injected at
% about 4.27 seconds.

% Data for Tim B.

clear all
close all

% =========================================================================
% Connect to Server: 
mdsconnect('mpexserver') 

%Extract data:
Spool = 9.5;
Shots = 8900 + [98:103,112:122,124];
Radius = [0:0.5:2.5,3.0:0.5:5,5.0,-0.5:-0.5:-3.0];
Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RA = [Stem,Branch];

Data{1} = [RA,'LP_V_RAMP'];
Data{2} = [RA,'TARGET_LP'];
Data{3} = [RA,'PWR_28GHz'];
Data{4} = [Stem,'SHOT_NOTE'];

for s = 1:length(Shots)
    [SHOT(s),~]=mdsopen('MPEX',Shots(s)); % see 8636

    [Vm{s},~]  = mdsvalue(Data{1});
    [Im{s},~]  = mdsvalue(Data{2});
    [Gyrotron{s},~]  = mdsvalue(Data{3});
    [NOTE{s},~]= mdsvalue(Data{4});

    [t_1{s},~]  = mdsvalue(['DIM_OF(',Data{1},')']);
    [t_2{s},~]  = mdsvalue(['DIM_OF(',Data{2},')']);
    [t_Gyrotron{s},~] = mdsvalue(['DIM_OF(',Data{3},')']);
end
% Create a function where we only need to give the name of the Stem, Branch
% and data name and it returns both the x and y arrays of the data.
% In addition, we could create a function that evaluates all the data or at
% least all the probes and tells you which datas are populated

mdsclose;
mdsdisconnect;
% =========================================================================
% DLP data analysis: 
switch Spool
    case 6.5
        L_tip = 1.3/1000;  % [m]
    case 9.5
        L_tip = 1.2/1000;  % [m]
    case 4.5
        L_tip = 1.2/1000;  % [m]
end
D_tip = 0.25/1000; % [m]
rp = D_tip/2; % [m]

AreaType = 2;
switch AreaType
    case 1
        Area = L_tip*pi*D_tip + 0.25*pi*D_tip^2;
    case 2
        Area = 2*L_tip*D_tip;
end

e_c = 1.602e-19; 
m_p = 1.6726e-27;
e_0 = 8.854187817*1e-12; % F/m
AMU = 2;

% Break data into cycles: -------------------------------------------------
% tStart = 4.1615; % [s]
% tEnd = 4.3115;
tStart = 4.16; % [s]
tEnd = 4.32;

for s = 1:length(Shots)
    
    rng = find(t_1{s}>=tStart & t_1{s}<=tEnd);
    t{s} = t_1{s}(rng);

    % Apply correct scale factors to data:
    Vm{s} = 22*Vm{s}(rng); 
    if 1
        Im{s} = sgolay_t(Im{s},3,7);
    end
    Im{s} = -2*Im{s}(rng)/150; % 150 Ohms resistor and a 2:1 voltage divider
    tm = 2;

    [locs,~] = peakseek(abs(Vm{s}),211);
    for c = 1:(length(locs)-1)
        I{s}{c} = Im{s}(locs(c):locs(c+1));
        V{s}{c} = Vm{s}(locs(c):locs(c+1));

        switch tm
            case 1
                time{s}(c) = ( t{s}(locs(c)) + t{s}(locs(c+1)) )/2;
            case 2
                time{s}(c) = t{s}(locs(c));
        end

    end

    x0 = [30e-3,2,0,0,0];

    % Fitting routine: --------------------------------------------------------
    for c = 1:length(V{s})
        A1 = V{s}{c};
        A2 = I{s}{c};
        res = @(x)   ( x(1)*tanh((A1-x(4))/x(2)) + (x(3)*(A1-x(4))) + x(5) ) - A2;
        Ifit{s}{c} = @(x)   x(1)*tanh((A1-x(4))/x(2)) + (x(3)*(A1-x(4))) + x(5); 

        [x1{s}(:,c),ssq,cnt] = LMFnlsq(res,x0,'MaxIter',5);
    end

    Te{s} = real(x1{s}(2,:));
    Isat{s} = real(x1{s}(1,:));
    Cs{s} = sqrt(e_c*Te{s}/(AMU*m_p));
    Ni{s} = real(Isat{s}./(0.61*e_c*Area*Cs{s}));
    Ld{s} = real(sqrt(e_0*Te{s}./(Ni{s}*e_c)));% Debye length
    Xi{s} = rp./Ld{s};
end
%%
t_HD_Gather_Start = 4.28;
t_HD_Gather_End   = 4.31;

t_LD_Gather_Start = 4.16;
t_LD_Gather_End   = 4.19;

for s = 1:length(Shots)

    rng2{s} = find(time{s}>=t_HD_Gather_Start & time{s}<=t_HD_Gather_End);
    Ni_HD(s) = mean(Ni{s}(rng2{s}));
    Te_HD(s) = mean(Te{s}(rng2{s}));
    
    rng3{s} = find(time{s}>=t_LD_Gather_Start & time{s}<=t_LD_Gather_End);
    Ni_LD(s) = mean(Ni{s}(rng3{s}));
    Te_LD(s) = mean(Te{s}(rng3{s}));
end

% =========================================================================
%% Plot data:
close all

figure; C = ['k','r','g','k','r','g','k','r','g','k','r','g','k','r','g','k','r','g'];
subplot(3,1,1); hold on
for s = 1:length(Shots)
    plot(time{s},Ni{s},C(s))
end
title('N_i')
ylim([0,8e19])

subplot(3,1,2); hold on
for s = 1:length(Shots)
    Te_f = sgolay_t(Te{s},3,11);
    Te_f = Te{s};
    plot(time{s},Te_f,C(s))
    plot(t_Gyrotron{1},Gyrotron{1})
end
title('Te')
ylim([0,12])
xlim([4.16,4.32])

subplot(3,1,3); hold on
for s = 1:length(Shots)
    plot(time{s},Ni{s}.*Te{s}*e_c,C(s))
end
title('P_e')
set(gcf,'position',[417 106 574 503])
ylim([0,60])

% =========================================================================
% DATA FOR TIM 
% =========================================================================

% BEFORE ECH:==============================================================

figure; hold on
for s = 1:length(Shots)
    rng = find(time{s}>=4.28 & time{s}<=4.31);
    %plot3(time{s}(rng),Radius(s)*ones(size(Ni{s}(rng))),Ni{s}(rng),C(s))
    Ni_mean(s) = mean(Ni{s}(rng));
    dNi_mean(s) = std(Ni{s}(rng));
end

figure; hold on
for s = 1:length(Shots)
    rng = find(time{s}>=4.28 & time{s}<=4.31);
    %plot3(time{s}(rng),Radius(s)*ones(size(Te{s}(rng))),Te{s}(rng),C(s))
    Te_mean(s) = mean(Te{s}(rng));
    dTe_mean(s) = std(Te{s}(rng));
end

figure; 
subplot(2,1,1); hold on
h(1) = errorbar(Radius, Ni_mean,dNi_mean,'ko')
ylim([0,8e19]);
xlim([-5,6])
ylabel('$n_e$ $[m^{-3}]$','interpreter','Latex')
set(gca,'Box','on')
legend(h,'Spool 6.5, TR2 = 310 A, PS1 = 5700 A, No helicon mode "jump"','Location','NorthWest')
set(gca,'XTick',[-5:1:6])

subplot(2,1,2); hold on
errorbar(Radius, Te_mean,dTe_mean,'ko')
ylim([0,20])
xlim([-5,6])
ylabel('$T_e$ $[V]$','interpreter','Latex')
xlabel('Vertical location [cm]')
set(gca,'XTick',[-5:1:6])

set(gcf,'color','w')
set(gca,'Box','on')

% Export data to Excel:
[~,n] = sort(Radius);
D = [Radius(n)',Ni_mean(n)',Te_mean(n)'];
F = {'R [cm]','n_i [m^-3]','T_e [eV]'};
FileName = 'During_ECH_NiTe_Spool_6_5.xlsx';
xlswrite(FileName,[F;num2cell(D)]);


% DURING ECH:==============================================================

%figure; hold on
for s = 1:length(Shots)
    rng = find(time{s}>=4.25 & time{s}<=4.265);
    %plot3(time{s}(rng),Radius(s)*ones(size(Ni{s}(rng))),Ni{s}(rng),C(s))
    Ni_mean(s) = mean(Ni{s}(rng));
    dNi_mean(s) = std(Ni{s}(rng));
end

%figure; hold on
for s = 1:length(Shots)
    rng = find(time{s}>=4.25 & time{s}<=4.265);
    %plot3(time{s}(rng),Radius(s)*ones(size(Te{s}(rng))),Te{s}(rng),C(s))
    Te_mean(s) = mean(Te{s}(rng));
    dTe_mean(s) = std(Te{s}(rng));
end

figure; 
subplot(2,1,1); hold on
h(1) = errorbar(Radius, Ni_mean,dNi_mean,'ko')
ylim([0,8e19]);
xlim([-5,6])
ylabel('$n_e$ $[m^{-3}]$','interpreter','Latex')
set(gca,'Box','on')
legend(h,'Spool 6.5, TR2 = 310 A, PS1 = 5700 A, No helicon mode "jump"','Location','NorthWest')
set(gca,'XTick',[-5:1:6])

subplot(2,1,2); hold on
errorbar(Radius, Te_mean,dTe_mean,'ko')
ylim([0,20])
xlim([-5,6])
ylabel('$T_e$ $[V]$','interpreter','Latex')
xlabel('Vertical location [cm]')
set(gca,'XTick',[-5:1:6])

set(gcf,'color','w')
set(gca,'Box','on')

% Export data to Excel:
[~,n] = sort(Radius);
D = [Radius(n)',Ni_mean(n)',Te_mean(n)'];
F = {'R [cm]','n_i [m^-3]','T_e [eV]'};
FileName = 'Before_ECH_NiTe_Spool_6_5.xlsx';
xlswrite(FileName,[F;num2cell(D)]);

%% Fit check:
if 1
    s = 1;
    figure; hold on
    plot(t{s},Vm{s})
    plot(t{s},Im{s}*1000)

    figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(V{s}{c},I{s}{c},'k')
        plot(V{s}{c},Ifit{s}{c}(x1{s}(:,c)),'r')
    end
end
