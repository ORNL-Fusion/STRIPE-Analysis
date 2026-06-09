% Radial DLP scan at spool 6.5 
% May 19th 2016, 28 GHz injection with mode converter installed into a high
% density mode jump D2 discharge.

% Logbook page 256 and 257
close all
clear all

mdsconnect('mpexserver')

L_tip = 1.2/1000;  % [m]
D_tip = 0.25/1000; % [m]
rp = D_tip/2; % [m]

AreaType = 2;
switch AreaType
    case 1
        Area = L_tip*pi*D_tip + 0.25*pi*D_tip^2;
    case 2
        Area = 2*L_tip*D_tip;
end
AMU = 2; % D2 gas

Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RA = [Stem,Branch];
Data{1} = [RA,'LP_V_RAMP']; % V
Data{2} = [RA,'TARGET_LP']; % I
Data{3} = [Stem,'FSCOPE:','TUBE07:','PMT_VOLT']; % 28 GHz pulse


Shots = 8800 + [66,67,70,72,73,74,76,77,78,79,80];
R = [0,1,2,3,3,4,5,-4,-2,-3,-3];

[Vm,t_v] = my_mdsvalue_v2(Shots,Data(1));
[Im,t_i] = my_mdsvalue_v2(Shots,Data(2));
[Tube,t_t] = my_mdsvalue_v2(Shots,Data(3));

tStart = 4.15; % [s]
tEnd = 4.32;
for s = 1:length(Shots)
    
    rng = find(t_v{s}>=tStart & t_v{s}<=tEnd);
    t{s} = t_v{s}(rng);

    % Apply correct scale factors to data:
    Vm{s} = 22*Vm{s}(rng); 
    if 0
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
% Plot data:
figure;
subplot(3,1,1); hold on
for s = 1:length(Shots)
    h(s) = plot(time{s},Ni{s})
end
ylim([0,5e19]); xlim([tStart,tEnd])
legend(h,'0 cm','+4 cm','+5 cm','-4 cm','location','NorthWest')
ylabel('n_i [m^{-3}]')

subplot(3,1,2); hold on
for s = 1:length(Shots)
    h(s) = plot(time{s},sgolay_t(Te{s},3,7));
end
ylim([0,18]); xlim([tStart,tEnd])
ylabel('T_e [eV]')

subplot(3,1,3); hold on
for s = 1:length(Shots)
    h(s) = plot(t_t{s},Tube{s})
end
xlim([tStart,tEnd])
ylabel('28 GHz pulse')

set(gcf,'color','w')

% Fit Test:
figure;
k = 3;
n = 45;
for s = 1:16
subplot(4,4,s); hold on
plot(V{k}{s+n},I{k}{s+n},'k')
plot(V{k}{s+n},Ifit{k}{s+n}(x1{k}(:,s+n)),'r')
title(num2str(time{k}(s+n)))
end

% Plotting radial profiles:
t_HD_Gather_Start = 4.295;
t_HD_Gather_End   = 4.305;

t_LD_Gather_Start = 4.18;
t_LD_Gather_End   = 4.2;

for s = 1:length(Shots)
    rng2{s} = find(time{s}>=t_HD_Gather_Start & time{s}<=t_HD_Gather_End);
    Ni_HD(s) = mean(Ni{s}(rng2{s}));
    Te_HD(s) = mean(Te{s}(rng2{s}));
    
    rng3{s} = find(time{s}>=t_LD_Gather_Start & time{s}<=t_LD_Gather_End);
    Ni_LD(s) = mean(Ni{s}(rng3{s}));
    Te_LD(s) = mean(Te{s}(rng3{s}));
end

%%
figure; 
subplot(2,1,1); hold on
plot(R,Ni_LD,'ko')
ylim([0,4e19])
subplot(2,1,2); hold on
plot(R,Te_LD,'ko')

figure; 
subplot(2,1,1); hold on
plot(R,Ni_HD,'ko')
ylim([0,4e19])
subplot(2,1,2); hold on
plot(R,Te_HD,'ko')
set(gcf,'color','w')

%% Export data to Excel
