% Calibration of ORNL puffer at spool 1.5
close all
clear all

V = 6.616; % [m^3] Volume of ProtoMPEX including balast tank

% =========================================================================
shotlist = 18000 + [754];
tEnd = 10;
% =========================================================================
shotlist = 18000 + [759,760,761,762,763,764]; % 3.0 to 6.0 V, 100 ms
X =                [3.5,4.0,4.5,5.0,5.5,6.0];
tEnd = 30; 
% Time recording time was increased to 30 seconds to enable the neutral gas
% pressure to equilibriate over the entire vacuum chamber

% The pressure increases linearly with puff voltage however it levels off
% after 5V due to depleting of the gas inside the reservour. every time we
% isolate the turbo pumps, the D2 bottle is isolated from the chamber so we
% cannot maintain constant pressure in the ORNL puffer.
% =========================================================================
shotlist = 18000 + [766,767,768,769,770,771,772 ]; 
X =                [100, 50, 25, 10,200,500,1000];
tEnd = 30; 
% =========================================================================
% shotlist = 18000 + 775;
% X = 100;

% Gas pressure at spool 2.5
address{1} = '\MPEX::TOP.MACHOPS1:PG2'; % PG2.5
address{2} = '\MPEX::TOP.MACHOPS1:PG4'; % PG4.5
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
close all
figure; hold on
s = 1;
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
%% Take steady state values
[a,b] = sort(X);
figure; hold on
h(1) = plot(X(b),pg2(b),'ko-')
h(2) = plot(X(b),pg4(b),'ro-')
h(3) = plot(X(b),pg6(b),'blo-')
h(4) = plot(X(b),pg10(b),'go-')

ylabel('$[Pa]$','Interpreter','latex','FontSize',13)
set(gcf,'color','w')
box on
legend(h,'PG 2.5','PG 4.5','PG 6.5','PG 10.5')

% xlim([0,100])
% ylim([0,0.06])
