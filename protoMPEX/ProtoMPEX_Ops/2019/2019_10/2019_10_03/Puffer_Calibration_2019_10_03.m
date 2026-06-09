% Calibration of ORNL puffer at spool 1.5
close all
clear all

V = [0.489875]; % [m^3] Volume of ProtoMPEX including balast tank

ShotType = 5;

switch ShotType
    case 1
        % =========================================================================
        % First pulse, all pumps ON
        shotlist = 28000 + [18  ]; 
        X =                [6.5];
        tEnd = 20;
    
    case 2
        % =========================================================================
        % Piezo 2.5 calibration shots, ALL pumps OFF
        shotlist = 28000 + [4  ]; 
        X =                [0.0];
        tEnd = 30; 
        
    case 3
        % =========================================================================
        % Piezo 2.5 calibration shots, ALL pumps OFF
        shotlist = 28000 + [4  ,3  ,8  ,9  ,10 ,11 ,12 ,13 ,14 ,15 ]; 
        X =                [0.0,6.0,4.0,4.0,3.0,5.0,7.0,4.5,3.5,2.5];
        tEnd = 30; 
        
    case 4
        % =========================================================================
        % Piezo 1.5 calibration shots, ALL pumps OFF
        shotlist = 28000 + [4  ,19 ,20 ,21 ,22 ,23 ,24  ,25 ,26 ]; 
        X =                [0.0,6.0,7.0,8.0,9.0,9.5,10.0,5.0,4.0];
        tEnd = 30; 
        
    case 5
        % Both Piezo 1.5 and 2.5 calibration shots, ALL pumps OFF
        % =========================================================================
        shotlist = 28000 + [4,27 ]; 
        X =                [0,6.5];
        Y =                [0,9.5];
        tEnd = 30; 
end

% Gas pressure at spool 2.5
address{1} = '\MPEX::TOP.MACHOPS1:PG2'; % PG8.5
address{2} = '\MPEX::TOP.MACHOPS1:PG4'; % PG2.5
address{3} = '\MPEX::TOP.MACHOPS1:PG3'; % PG6.5
address{4} = '\MPEX::TOP.MACHOPS1:PG1'; % PG10.5

[s2,t_s2] = my_mdsvalue_v2(shotlist,address(1));
[s4,t_s4] = my_mdsvalue_v2(shotlist,address(2));
[s6,t_s6] = my_mdsvalue_v2(shotlist,address(3));
[s10,t_s10] = my_mdsvalue_v2(shotlist,address(4));

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
s = 2;
for s = 1:numel(shotlist);
    figure;
    hold on
    h{s}(1) = plot(t_PG2{s},PG2{s},'k')
    h{s}(2) = plot(t_PG4{s},PG4{s},'r')
    h{s}(3) = plot(t_PG6{s},PG6{s},'bl')
    h{s}(4) = plot(t_PG10{s},PG10{s},'g')
    legend(h{s},'8.5','4.5','6.5','10.5')
    legend(h{s},'PG2','PG1','PG3','PG4')
    xlim([3.5,tEnd])
    ylim([-0.05,2]); 
    ylabel('$[Pa]$','Interpreter','latex','FontSize',13)
    title(['Shot: ',num2str(shotlist(s)),' ,V = ',num2str(X(s))])
    set(gcf,'color','w')
    box on
    grid on
    xlabel('time [sec]')
end
%% Take steady state values

[a,b] = sort(X);
figure; hold on
h{s}(1) = plot(X(b),pg2(b),'ko-')
h{s}(2) = plot(X(b),pg4(b),'ro-')
h{s}(3) = plot(X(b),pg6(b),'blo-')
h{s}(4) = plot(X(b),pg10(b),'go-')

ylabel('$[Pa]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
box on
legend(h{s},'PG 2.5','PG 4.5','PG 6.5','PG 10.5')
legend(h{s},'PG2'   ,'PG1'   ,'PG3'   ,'PG4'    )
grid on
xlabel('Request Voltage [V]')

[a,b] = sort(X);
figure; hold on
dt = 0.1;
h{s}(1) = plot(X(b),2*(pg2(b) -pg2(1) )*V/(dt*k_B*300),'ko-')
h{s}(2) = plot(X(b),2*(pg4(b) -pg4(1) )*V/(dt*k_B*300),'ro-')
h{s}(3) = plot(X(b),2*(pg6(b) -pg6(1) )*V/(dt*k_B*300),'blo-')
h{s}(4) = plot(X(b),2*(pg10(b)-pg10(1))*V/(dt*k_B*300),'go-')
grid on
ylabel('$D_0$ $[s^{-1}]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
box on
legend(h{s},'PG 2.5','PG 4.5','PG 6.5','PG 10.5')
legend(h{s},'PG2'   ,'PG1'   ,'PG3'   ,'PG4'    )
ylim([0,3e21])
xlabel('Request Voltage [V]')
title('Gas flow rate: $D_0$ $[s^{-1}]$','Interpreter','latex','FontSize',13)
