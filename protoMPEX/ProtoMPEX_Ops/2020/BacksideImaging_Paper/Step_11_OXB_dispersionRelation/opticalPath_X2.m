% X-mode 2nd harmonic absorption length
% Based on reference:
% Electron cyclotron resonance heating of weakly relativistic plasmas
%K. R. Chu, and B. Hui
%Citation: The Physics of Fluids 26, 69 (1983); doi: 10.1063/1.863986

clear all
close all

% Magnetic field:
B0 = 0.4; 

% Temperature:
Te = logspace(1,4);

% Frequency:
freq = 28E9;
w = 2*pi*freq;

% Dimensionless parameters:
X = 0.5;
Y = w_ce(B0)/w;

% Characteristic length:
Rp = 2/100;
legnthCase = 1;
switch legnthCase
    case 1
        dB = B0/10;
        L  = Rp*B0/dB;
    case 2
        L = 2*Rp;
end

% Derived quantities:
W = X*(Y^(-2));
E1 = (6-W)^(5/2);
E2 = (2-W)^(1/2);
E3 = (3-W)^(-5/2);
G  = (e_c*Te/(m_e*c_light^2));
k0 = w/c_light;

% Absorption factor:
eta_X2 = k0*L*(pi/8)*W*Y*E1*E2*E3.*G;

% Absorbed fraction:
f_transmitted = exp(-eta_X2);

% Transmission factor:
f_absorbed = 1 - f_transmitted;

% Plot results:
figure('color','w')
plot(Te,f_transmitted)
set(gca,'Xscale','log')
set(gca,'Yscale','lin')
title('Transmitted power')

figure('color','w')
plot(Te,eta_X2)
set(gca,'Xscale','log')
set(gca,'Yscale','log')
title('Absorption coefficient vs T_e')
