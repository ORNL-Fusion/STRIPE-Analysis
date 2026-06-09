% Calibration of ORNL puffer at spool 1.5
close all
clear all

% Need to classify data into two classes

CMPT = 0;
V_chamber = 6.616; % [m^3] Volume of ProtoMPEX including balast tank
tEnd = 30; 

if CMPT == 1
    shotlist = 18000 + [759,760,761,762,763,764,766,767,768,769,770,771,772 ];
    X =                [3.5,4.0,4.5,5.0,5.5,6.0,100, 50, 25, 10,200,500,1000];
    % Gas pressure at spool 2.5
    address{1} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
    address{2} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5
    address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
    address{4} = '\MPEX::TOP.MACHOPS1:PG1'; % PG10.5

    [s2,t_s2] = my_mdsvalue_v2(shotlist,address(1));
    [s4,t_s4] = my_mdsvalue_v2(shotlist,address(2));
    [s6,t_s6] = my_mdsvalue_v2(shotlist,address(3));
    [s10,t_s10] = my_mdsvalue_v2(shotlist,address(4));
    
    save('Pressure_2017_12_20.mat')
else
    load('Pressure_2017_12_20.mat')
end

% Arrange data
n0 = 1e3;
for s = 1:length(shotlist)
    PG2{s} = (s2{s}-mean(s2{s}(1:n0)))*2/7.5;
    t_PG2{s} = t_s2{s}(1:end-1);
    pg2(s) = mean(PG2{s}(end-500:end));
    
    PG4{s} = (s4{s}-mean(s4{s}(1:n0)))*10/7.5;
    t_PG4{s} = t_s4{s}(1:end-1);
    pg4(s) = mean(PG4{s}(end-500:end));

    
    PG6{s} = (s6{s}-mean(s6{s}(1:n0)))*2/7.5;
    t_PG6{s} = t_s6{s}(1:end-1);
    pg6(s) = mean(PG6{s}(end-500:end));

        
    PG10{s} = (s10{s}-mean(s10{s}(1:n0)))*2/7.5;
    t_PG10{s} = t_s10{s}(1:end-1);
    pg10(s) = mean(PG10{s}(end-500:end));
end
%% Plot data
close all
figure; hold on
s = 13;
plot(t_PG2{s},PG2{s},'k')
plot(t_PG4{s},PG4{s},'r')
plot(t_PG6{s},PG6{s},'bl')
plot(t_PG10{s},PG10{s},'g')
xlim([3.5,tEnd])
ylim([-0.05,0.7]); 
ylabel('$[Pa]$','Interpreter','latex','FontSize',13)
title(num2str(shotlist(s)))
set(gcf,'color','w')
box on
%% Arrange into two groups
% Voltage scan:
rng = find(X<= 5 | X == 100);
Volts = [3.5,4.0,4.5,5.0,6.0];
p2{1} = pg2(rng);
p4{1} = pg4(rng);
p6{1} = pg6(rng);
p10{1} = pg10(rng);

% Time scan:
rng = find(shotlist>=18766);
dT = [100, 50, 25, 10,200,500,1000];
p2{2} = pg2(rng);
p4{2} = pg4(rng);
p6{2} = pg6(rng);
p10{2} = pg10(rng);

%% Plotting data
close all

% =========================================================================
% Voltage scan:
% Fit line:
RngData = 1:5; 
p_v =  (p2{1}(RngData) + p4{1}(RngData)  + p6{1}(RngData) )/3;
c_v = polyfit(Volts(RngData),p_v,1);
x_v = linspace(0,10);

figure; hold on
hV(6) = plot(x_v,polyval(c_v,x_v),'g-','lineWidth',2)
hV(1) = plot(Volts,p2{1},'ko-')
hV(2) = plot(Volts,p4{1},'ro-')
hV(3) = plot(Volts,p6{1},'blo-')
hV(4) = plot(Volts,p10{1},'go-')
hV(5) = plot(Volts(RngData),p_v,'msq-')

xlim([0,10])
ylim([0,0.1])
ylabel('$[Pa]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
box on
legend(hV,{'PG 2.5','PG 4.5','PG 6.5','PG 10.5','Mean',['slope: ',num2str(c_v(1),3),' [V/Pa]']},'location','NorthWest')
grid on

% =========================================================================
% Time scan:
% Fit line:
RngData = 1:7; 
p_t =  (p2{2}(RngData) + p4{2}(RngData)  + p6{2}(RngData) + p10{2}(RngData) )/4;
c_t = polyfit(dT(RngData),p_t,1);
x_t = linspace(0,2000);
[a,b] = sort(dT);

figure; hold on
hT(6) = plot(x_t,polyval(c_t,x_t),'g-','lineWidth',2)
hT(1) = plot(dT(b),p2{2}(b),'k.-')
hT(2) = plot(dT(b),p4{2}(b),'r.-')
hT(3) = plot(dT(b),p6{2}(b),'bl.-')
hT(4) = plot(dT(b),p10{2}(b),'g.-')
hT(5) = plot(dT,p_t,'mo')


ylabel('$[Pa]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
set(gca,'Yscale','lin','Xscale','log')
box on
legend(hT,{'PG 2.5','PG 4.5','PG 6.5','PG 10.5','Mean',['slope: ',num2str(c_t(1),3),' [s/Pa]']},'location','NorthWest')
grid on
xlim([0,1000])
ylim([0,0.6])

%% Plotting the D2 flow rate [D2/s] at the 6 V level
close all
% V = 6 V, various dT
Tg = 300; % D2 gas at room temp 300 K
n_dot_D2_6V = (p_t./(dT*1e-3))*V_chamber/(k_B*Tg);
n_dot_atoms_6V = 2*n_dot_D2_6V;

% dT = 100 ms, various Voltages
n_dot_D2 = (p_v./(0.1))*V_chamber/(k_B*Tg);
n_dot_atoms = 2*n_dot_D2;

figure; 
subplot(1,2,1);hold on
hF(1) = plot(dT,n_dot_D2_6V*1e-20,'ko')
hF(2) = plot(dT,n_dot_atoms_6V*1e-20,'ro')

ylim([0,20])
title('Puffer at 6 V, $T_{D_2}$ at 300 K','Interpreter','latex','FontSize',13)
ylabel('$[s^{-1}]$ $\times 10^{20}$','Interpreter','latex','FontSize',13)
xlabel('$\Delta{t} [ms]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
set(gca,'Yscale','lin','Xscale','log')
box on
legend(hF,'D_2','Atoms','location','NorthWest')
grid on

subplot(1,2,2);hold on
hG(1) = plot(Volts,n_dot_D2*1e-20,'ko')
hG(2) = plot(Volts,n_dot_atoms*1e-20,'ro')

ylim([0,20])
title('100 ms puffs, $T_{D_20}$ at 300 K','Interpreter','latex','FontSize',13)
ylabel('$[s^{-1}]$ $\times 10^{20}$','Interpreter','latex','FontSize',13)
xlabel('$\Delta{V} [V]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
set(gca,'Yscale','lin','Xscale','lin')
box on
legend(hG,'D_2','Atoms','location','NorthWest')
grid on


% note:
% During normal operations, we puff at a voltage of 3.5 V at the end of the
% pulse during the steady state part. we have 4.3-3.6 V at the start of the
% pulse and 3.5 V at the end of the pulse.
% this means that during steady state we are injecting about 6e20 atoms/sec
% For the low pressure case: upstream puffing, ORNL spool 1.5 puffer,
% lowest pressure case.

% Suppose we have a flow of mach 0.5, Te is about 2 eV and ne is triangular
% with a max ne of 6e19 m-3 we have 
ne_mean = 3e19;
Te_mean = 2; 
Uz = 0.55*C_s(Te_mean,2);
PlasmaRadius = 3/100; 
Area = pi*PlasmaRadius^2;
FluxDensity = ne_mean*Uz;
PlasmaFlux = FluxDensity*Area

%% Figure for TOFE 2018 paper

figure;
hold on
hG(1) = plot(Volts,n_dot_D2*1e-20,'ko-')
hG(2) = plot(Volts,n_dot_atoms*1e-20,'ro-')

ylim([0,20])
title('100 ms puffs, $T_{D_2}$ at 300 K','Interpreter','latex','FontSize',13)
ylabel('$[s^{-1}]$ $\times 10^{20}$','Interpreter','latex','FontSize',13)
xlabel('$\Delta{V} [V]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
set(gca,'Yscale','lin','Xscale','lin')
box on
legend(hG,'D_2','Atoms','location','NorthWest')
grid on

