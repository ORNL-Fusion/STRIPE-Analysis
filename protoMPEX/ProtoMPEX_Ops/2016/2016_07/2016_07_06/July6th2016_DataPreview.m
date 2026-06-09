clear all
close all

% =========================================================================
% Connect to Server: 
mdsconnect('mpexserver') 

%Extract data:
Spool = 6.5;
Shots = 9200 + [90, 97, 98, 99, 100,103,106,106,107,108,109,111,113,114,115,116];
R     =        [0 ,0.5, 1 ,1.5,  2 ,2.5,2.5,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5]; 
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
    case 10.5
        L_tip = 1.2/1000;  % [m]
end
D_tip = 0.25/1000; % [m]
rp = D_tip/2; % [m]

AreaType = 1;
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
%%
% Before ECH pulse
t_A_Gather_Start = 4.24;
t_A_Gather_End   = 4.264;

% During ECH pulse
t_B_Gather_Start = 4.27;
t_B_Gather_End   = 4.292;

for s = 1:length(Shots)

    rngA{s} = find(time{s}>=t_A_Gather_Start & time{s}<=t_A_Gather_End);
    Ni_A(s) = mean(Ni{s}(rngA{s}));
    d_Ni_A(s) = std(Ni{s}(rngA{s}),1);

    Te_A(s) = mean(Te{s}(rngA{s}));
    d_Te_A(s) = std(Te{s}(rngA{s}),1);

    
    rngB{s} = find(time{s}>=t_B_Gather_Start & time{s}<=t_B_Gather_End);
    Ni_B(s) = mean(Ni{s}(rngB{s}));
    d_Ni_B(s) = std(Ni{s}(rngB{s}),1);

    Te_B(s) = mean(Te{s}(rngB{s}));
    d_Te_B(s) = std(Te{s}(rngB{s}),1);
end

close all

figure;
subplot(2,1,1)
s = find(R==2.5); s = s(1);
plot(time{s},Ni{s})
subplot(2,1,2); hold on
plot(time{s},Te{s})
plot(t_Gyrotron{s},Gyrotron{s})
ylim([0,12])
xlim([4.16,4.32])

figure; 
subplot(2,1,1); hold on
h(1) = errorbar(R,Ni_B,d_Ni_B,'r')
h(2) = errorbar(R,Ni_A,d_Ni_A,'bl')
set(gca,'box','on')
xlim([-0.5,7])
ylabel('$ n_e $','Interpreter','latex','Rotation',0,'HorizontalAlignment','right','FontSize',12)
legend(h,'During ECH','Before ECH')

subplot(2,1,2); hold on
h(1) = errorbar(R,Te_B,d_Te_B,'r')
h(2) = errorbar(R,Te_A,d_Te_A,'bl')
set(gca,'box','on')
ylim([0,20]); xlim([-0.5,7])
ylabel('$ T_e $','Interpreter','latex','Rotation',0,'HorizontalAlignment','right','FontSize',12)

set(gcf,'color','w','Position',[360 107 444 510])

% Export data to Excel:
% Before the ECH pulse
[~,n] = sort(R);
D_A = [R(n)',Ni_A(n)',d_Ni_A(n)',Te_A(n)',d_Te_A(n)'];
F = {'R [cm]','n_i [m^-3]','d_n_i [m^-3]','T_e [eV]','d_T_e [eV]'};
FileName = 'Before_ECH_NiTe_Spool_6_5.xlsx';
xlswrite(FileName,[F;num2cell(D_A)]);

% During the ECH pulse:
D_B = [R(n)',Ni_B(n)',d_Ni_B(n)',Te_B(n)',d_Te_B(n)'];
F = {'R [cm]','n_i [m^-3]','d_n_i [m^-3]','T_e [eV]','d_T_e [eV]'};
FileName = 'During_ECH_NiTe_Spool_6_5.xlsx';
xlswrite(FileName,[F;num2cell(D_B)]);

% =========================================================================
%% Plot data:

figure; C = ['k','r','g','k','r','g','k','r','g','k','r','g','k','r','g','k','r','g'];
subplot(3,1,1); hold on
for s = 1:length(Shots)
    plot(time{s},Ni{s},C(s))
end
title('N_i')
ylim([0,8e19])

subplot(3,1,2); hold on
for s = 1:length(Shots)
    %Te_f = sgolay_t(Te{s},3,11);
    Te_f = Te{s};
    plot(time{s},Te_f,C(s))
    plot(t_Gyrotron{end},Gyrotron{end})
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


%% Fit check:
if 1
    s = 1;
    figure; hold on
    plot(t{s},Vm{s})
    plot(t{s},Im{s}*1000)

end

if 1
    
     figure;
    for c = 1:25;
        subplot(5,5,c); hold on
        plot(V{s}{c},I{s}{c},'k')
        plot(V{s}{c},Ifit{s}{c}(x1{s}(:,c)),'r')
    end
end
