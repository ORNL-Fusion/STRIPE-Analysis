clear all
close all

% =========================================================================
% Connect to Server: 
mdsconnect('mpexserver') 

%Extract data:
Shots = [8516,8517,8518];
Spool = [9.5,10.5,6.5];
z =     [3.55,3.9,2.7]; % approximate location, spool 4.5 is at 2.1 m

Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RA = [Stem,Branch];

Data{1} = [RA,'LP_V_RAMP'];
Data{2} = [RA,'TARGET_LP'];
Data{3} = [RA,'PG1'];
Data{4} = [Stem,'SHOT_NOTE'];

for s = 1:length(Shots)
    [SHOT(s),~]=mdsopen('MPEX',Shots(s)); % see 8636

    [Vm{s},~]  = mdsvalue(Data{1});
    [Im{s},~]  = mdsvalue(Data{2});
    [Pg{s},~]  = mdsvalue(Data{3});
    [NOTE{s},~]= mdsvalue(Data{4});

    [t_1{s},~]  = mdsvalue(['DIM_OF(',Data{1},')']);
    [t_2{s},~]  = mdsvalue(['DIM_OF(',Data{2},')']);
    [t_Pg{s},~] = mdsvalue(['DIM_OF(',Data{3},')']);
end
% Create a function where we only need to give the name of the Stem, Branch
% and data name and it returns both the x and y arrays of the data.
% In addition, we could create a function that evaluates all the data or at
% least all the probes and tells you which datas are populated

mdsclose;
mdsdisconnect;
% =========================================================================
% DLP data analysis: 
L_tip = 1.2/1000;  % [m]
D_tip = 0.254/1000; % [m]
rp = D_tip/2; % [m]

AreaType = 1;
switch AreaType
    case 1
        Area = L_tip*pi*D_tip + 0.25*pi*D_tip^2; % surface area
    case 2
        Area = 2*L_tip*D_tip; % projection of area on to B0
end

e_c = 1.602e-19; 
m_p = 1.6726e-27;
e_0 = 8.854187817*1e-12; % F/m
AMU = 2;

% Break data into cycles: -------------------------------------------------
% tStart = 4.1615; % [s]
% tEnd = 4.3115;
tStart = 4.15; % [s]
tEnd = 4.32;

for s = 1:length(Shots)
    
    rng = find(t_1{s}>=tStart & t_1{s}<=tEnd);
    t{s} = t_1{s}(rng);

    % Apply correct scale factors to data:
    Vm{s} = 2*12.05*Vm{s}(rng); 
    if 1
        Im{s} = sgolay_t(Im{s},3,7);
    end
    Im{s} = -2*Im{s}(rng)/142.5; % 150 Ohms resistor and a 2:1 voltage divider
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

    Te{s} = real(x1{s}(2,:))/2;
    Isat{s} = real(x1{s}(1,:));
    Cs{s} = sqrt(e_c*Te{s}/(AMU*m_p));
    Ni{s} = real(Isat{s}./(0.61*e_c*Area*Cs{s}));
    Ld{s} = real(sqrt(e_0*Te{s}./(Ni{s}*e_c)));% Debye length
    Xi{s} = rp./Ld{s};
end
% =========================================================================
%% Plot data:
close all

figure; C = ['k','r','g'];
subplot(3,1,1); hold on
for s = 1:3
    plot(time{s},Ni{s},C(s))
end
title('N_i')
ylim([0,8e19])

subplot(3,1,2); hold on
for s = 1:3
    plot(time{s},Te{s},C(s))
end
title('Te')
ylim([0,12])

subplot(3,1,3); hold on
for s = 1:3
    plot(time{s},Ni{s}.*Te{s}*e_c,C(s))
end
title('P_e')
set(gcf,'position',[417 106 574 503])
ylim([0,60])

tLD = 4.2;
tHD = 4.3;
% Spool 6.5;
n = find(time{3}>=tHD); n = n(1);
N_e_HD(1) = Ni{3}(n); T_e_HD(1) = Te{3}(n); 
n = find(time{3}>=tLD); n = n(1);
N_e_LD(1) = Ni{3}(n); T_e_LD(1) = Te{3}(n); 

% Spool 9.5;
n = find(time{1}>=tHD); n = n(1);
N_e_HD(2) = Ni{1}(n); T_e_HD(2) = Te{1}(n); 
n = find(time{1}>=tLD); n = n(1);
N_e_LD(2) = Ni{1}(n); T_e_LD(2) = Te{1}(n);

% Spool 10.5;
n = find(time{2}>=tHD); n = n(1);
N_e_HD(3) = Ni{2}(n); T_e_HD(3) = Te{2}(n); 
n = find(time{2}>=tLD); n = n(1);
N_e_LD(3) = Ni{2}(n); T_e_LD(3) = Te{2}(n); 

P_e_HD = e_c.*N_e_HD.*T_e_HD;
P_e_LD = e_c.*N_e_LD.*T_e_LD;
z = [2.7, 3.55, 3.9];

figure; 
subplot(3,1,1); hold on
plot(z,N_e_HD,'ko-')
plot(z,N_e_LD,'ro-')
line([1.6,1.9],[0,0],'color','k','LineWidth',4)
xlim([0,5])
title('ne')

subplot(3,1,2); hold on
plot(z,T_e_HD,'ko-')
plot(z,T_e_LD,'ro-')
line([1.6,1.9],[0,0],'color','k','LineWidth',4)
xlim([0,5])
title('Te')

subplot(3,1,3); hold on
plot(z,P_e_HD,'ko-')
plot(z,T_e_LD,'ro-')
line([1.6,1.9],[0,0],'color','k','LineWidth',4)
xlim([0,5])
title('Pe')
xlabel('z [cm]')

% Spool 4.5
line([2.1,2.1],[0,100],'LineStyle',':','color','k','LineWidth',1)
% Spool 2.5
line([1.4,1.4],[0,100],'LineStyle',':','color','k','LineWidth',1)
% Spool 1.5
line([1.1,1.1],[0,100],'LineStyle',':','color','k','LineWidth',1)

% Target location:
line([4.2,4.2],[0,10],'LineStyle','-','color','k','LineWidth',2)


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
