% Exploring the gas-dynamic effect in Proto-MPEX
clear all 
close all

% The data for this code can be found on 2018_11_26

% Source reference density
n1_0 = 7.47e19;
% Target reference density
n2_0 = 2.8e19;

% Magnetic field at Target DLP
B2 = 0.3;
dB2 = 0.05;
% Magnetic field at Dump DLP
B1 = 0.23;
dB1 = 0.1;
% Magnetic field at coils 7-8
Bm = [0.5, 0.72, 0.94, 1.17, 1.39, 1.61];
% Magnetic field at the source:
Bs = 0.05; 
% Mirror ratio:
R = Bm/Bs;

% Plasma density measured at Target DLP
n2_i_xp = [2.8,2.5,2.2,1.8,1.5,1.3]*1e19;

figure
plot(R,n2_i_xp,'ro')
ylim([0,0.5e20])
xlim([0,50])

% Sound speed 
% Target:
c2      = 1.0*[  1,  1,  1,  1,  1,  1];
dc1 = 0.1; 
c2_0 = c2(1);
% Source:
c1      = 0.5*[  1,  1,  1,  1,  1,  1];
dc2 = 0.1;

% Gas dynamic prediction of source plasma density
n1_i = n1_0 + (B1/B2).*((n2_0*c2_0 - n2_i_xp.*c2)./c1);

dn1_i = (n1_i - n1_0).*sqrt( (dB1/B1)^2 + (dB2/B2)^2 + (dc1/c1(1))^2 + (dc2/c2(1))^2 );

% Plasma density measured at the source
n1_i_xp = [7.47, 7.73, 8.07, 8.52, 8.72, 9.27]*1e19;


% Figure of raw data:
% Here plot the density measurements 

FontSize = 16;
figure; 
hold on
hn(2) = errorbar(R,n1_i   ,dn1_i,'r');
set(hn(2),'LineWidth',2)
hn(1) = plot(R,n1_i_xp,'ksq-');
set(hn(1),'LineWidth',1)
set(hn(1),'MarkerSize',9)
set(gca,'FontName','Times','FontSize',FontSize)
ylim([0,1.3e20])
xlim([5,35])
xlabel('Mirror ratio R','Interpreter','latex','FontSize',FontSize)
ylabel('$n_e$ $[m^{-3}]$','Interpreter','latex','FontSize',FontSize+3)
box on
grid on
L = legend(hn,'Experiment, Source','Theory');
L.Box = 'off';
set(gcf,'position',[360.3333  279.0000  482.0000  338.6667],'color','w')
L.Location = 'SouthEast';