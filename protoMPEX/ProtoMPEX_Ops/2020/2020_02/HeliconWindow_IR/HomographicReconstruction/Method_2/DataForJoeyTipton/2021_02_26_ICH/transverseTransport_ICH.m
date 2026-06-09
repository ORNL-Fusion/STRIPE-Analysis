% Calculate the cross field transport under the ICH window;

clear all
close all

% Plasma conditions:
Ti = [25,30]; 
A = 2;
mi = A*m_p;
ni = 6e19;
B0 = 0.1; 

% Geometry:
L = 0.25;  % [m]
R = 3/100; % [m]

% Derived quantities:
uTi = sqrt(e_c*Ti/mi);

% Frequencies:
nuie = nu_ie(ni,mi,uTi);
wci = w_ci(B0,A);

% Diffusion coeficients:
ratio = wci./nuie;
D_par = e_c*Ti./(mi*nuie);
D_per = D_par.*ratio.^(-2);

% Particle fluxes:
S_per = (ni/R)*4*pi*R*L*D_per; % [s^-1]

% Power flux:
Ei = e_c.*Ti;
Q_per = Ei.*S_per

%% functions:
function y = nu_ie(ni,mi,ui)
% Assume singly ionized, Z = 1;

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
e_0 = 8.8542e-12;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c = 3e8;

logA = 12; 

y = ni*(e_c^4)*logA./(2*pi*((mi*e_0)^2)*ui.^3);

end