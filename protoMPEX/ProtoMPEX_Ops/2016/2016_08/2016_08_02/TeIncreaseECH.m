clear all
close all

% =========================================================================
% Connect to Server: 
mdsconnect('mpexserver') 

%Extract data:
Spool = 9.5;

%ECH On and Off: 30:34
%ECH change timing: 33:36
Shots = 9600 + [30:36];

Stem = '\MPEX::TOP.';
Branch = 'MACHOPS1:';
RA = [Stem,Branch];

Data{1} = [RA,'LP_V_RAMP'];
Data{2} = [RA,'TARGET_LP'];
Data{3} = [RA,'PWR_28GHz'];
Data{4} = [Stem,'SHOT_NOTE'];

for s = 1:length(Shots)
    [SHOT(s),~]=mdsopen('MPEX',Shots(s)); %

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
    case 10.5
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

    Te{s} = real(x1{s}(2,:))/2;
    Isat{s} = real(x1{s}(1,:));
    Cs{s} = sqrt(e_c*Te{s}/(AMU*m_p));
    Ni{s} = real(Isat{s}./(0.61*e_c*Area*Cs{s}));
    Ld{s} = real(sqrt(e_0*Te{s}./(Ni{s}*e_c)));% Debye length
    Xi{s} = rp./Ld{s};
end

%% Plot data:

C = ['k','r','r','k','r','g','k','r','g','k','r','g'];
TimeStart = 4.16; 
TimeEnd = 4.32;

mode = 2; 
% ECH On and Off
figure; 
subplot(2,2,1); hold on
switch mode
    case 1
for s = 1:length(Shots(1:4))
    h(s) = plot(time{s},Ni{s},C(s))
end
    case 2
    h(1) = plot(time{1},(Ni{1}+Ni{4})/2,C(1));
    h(2) = plot(time{2},(Ni{2}+Ni{3})/2,C(2));
    plot(t_Gyrotron{1},Gyrotron{1}*0.7e19,C(1))
end     
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
xlim([TimeStart,TimeEnd])
ylim([0,6e19])
ax(1) = gca;
legend([h(1:2)],'ECH On, y = +1.5 cm','ECH Off, y = +1.5 cm','ECH Off','location','NorthWest')

subplot(2,2,2); hold on
switch mode
    case 1
for s = 1:length(Shots(1:4))
    Te_f = Te{s};
    plot(time{s},Te_f,C(s))
    plot(t_Gyrotron{s},Gyrotron{s},C(s))
end
    case 2
    h(1) = plot(time{1},(Te{1}+Te{4})/2,C(1));
    h(2) = plot(time{2},(Te{2}+Te{3})/2,C(2));
    plot(t_Gyrotron{1},Gyrotron{1},C(1))
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,8])
xlim([TimeStart,TimeEnd])
ax(2) = gca;

subplot(2,2,3); hold on
switch mode
    case 1
for s = 1:length(Shots(1:4))
    plot(time{s},Ni{s}.*Te{s}*e_c,C(s))
end
    case 2
    Pe{1} = Ni{1}.*Te{1}*e_c;
    Pe{2} = Ni{2}.*Te{2}*e_c;
    Pe{3} = Ni{3}.*Te{3}*e_c;
    Pe{4} = Ni{4}.*Te{4}*e_c;
    
    h(1) = plot(time{1},(Pe{1}+Pe{4})/2,C(1));
    h(2) = plot(time{2},(Pe{2}+Pe{3})/2,C(2));
    plot(t_Gyrotron{1},3*Gyrotron{1},C(1))
end
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
set(gcf,'position',[417 106 574 503])
ylim([0,30])
xlim([TimeStart,TimeEnd])
ax(3) = gca;

set(gcf,'color','w')
set(ax,'Box','on')

%% Moving ECG timing

TimeStart = 4.16; 
TimeEnd = 4.32;

mode = 2; 
% ECH On and Off
figure; 
subplot(2,2,1); hold on
switch mode
    case 1
for s = 4:7
    h(s) = plot(time{s},Ni{s},C(s));
end
    case 2
    h(1) = plot(time{4},(Ni{4}+Ni{5})/2,C(1));
    plot(t_Gyrotron{4},Gyrotron{4}*0.7e19,C(1))

    h(2) = plot(time{6},(Ni{6}+Ni{7}(1:62))/2,C(2));
    plot(t_Gyrotron{7},Gyrotron{7}*0.7e19,C(2))
end     
title('$ n_e $ $ [m^{-3}] $','interpreter','Latex','FontSize',13,'Rotation',0)
xlim([TimeStart,TimeEnd])
ylim([0,6e19])
ax(1) = gca;
legend([h(1:2)],'ECH On, y = +1.5 cm','ECH Off, y = +1.5 cm','ECH Off','location','NorthWest')

subplot(2,2,2); hold on
switch mode
    case 1
for s = 2:length(Shots(2:6))
    Te_f = Te{s};
    plot(time{s},Te_f,C(s))
    plot(t_Gyrotron{s},Gyrotron{s},C(s))
end
    case 2    
    h(1) = plot(time{4},(Te{4}+Te{5})/2,C(1));
    plot(t_Gyrotron{4},Gyrotron{4},C(1))

    h(2) = plot(time{6},(Te{6}+Te{7}(1:62))/2,C(2));
    plot(t_Gyrotron{7},Gyrotron{7},C(2))
end
title('$ T_e $ $ [eV] $','interpreter','Latex','FontSize',13,'Rotation',0)
ylim([0,8])
xlim([TimeStart,TimeEnd])
ax(2) = gca;

subplot(2,2,3); hold on
switch mode
    case 1
for s = 2:length(Shots(2:6))
    plot(time{s},Ni{s}.*Te{s}*e_c,C(s))
end
    case 2
    Pe{4} = Ni{4}.*Te{4}*e_c;
    Pe{5} = Ni{5}.*Te{5}*e_c;
    Pe{6} = Ni{6}.*Te{6}*e_c;
    Pe{7} = Ni{7}(1:62).*Te{7}(1:62)*e_c;
    
        
    h(1) = plot(time{4},(Pe{4}+Pe{5})/2,C(1));
    plot(t_Gyrotron{4},3*Gyrotron{4},C(1))

    h(2) = plot(time{6},(Pe{6}+Pe{7})/2,C(2));
    plot(t_Gyrotron{7},3*Gyrotron{7},C(2))
end
title('$ P_e $ $ [Pa] $','interpreter','Latex','FontSize',13,'Rotation',0)
set(gcf,'position',[417 106 574 503])
ylim([0,40])
xlim([TimeStart,TimeEnd])
ax(3) = gca;

set(gcf,'color','w')
set(ax,'Box','on')

%% Fit check:
if 1
    figure;     hold on
    for s = 1:length(Shots);
    plot(t{s},Vm{s})
    plot(t{s},Im{s}*1000)
    end
    grid on
    ylim([-100,100])
end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(V{s}{c},I{s}{c},'k')
        plot(V{s}{c},Ifit{s}{c}(x1{s}(:,c)),'r')
    end
end
