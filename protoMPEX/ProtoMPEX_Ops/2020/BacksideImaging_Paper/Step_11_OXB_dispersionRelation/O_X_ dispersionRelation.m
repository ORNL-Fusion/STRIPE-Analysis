% Understanding the O-X mode conversion process:
clear all
close all

% Define refractive index
N = 1;

% Define incidence angle:
theta = 30*pi/180;
N_par = N*cos(theta);

% Define dimensionless parameters:
X = linspace(0,2);
Y = 0.5;

% Assemble Stix terms:
S = 1 - X*(1/(1 - Y^2));
P = 1 - X;
R = 1 - X*(1/(1 - Y));
L = 1 - X*(1/(1 + Y));

% Assemble terms of biquadratic dispersion relation:
Aperp = S;
Bperp = R*L + R*S - (N_par^2)*(P + S); 
Cperp = P*(N_par^2 - R)*(N_par^2 - L);

% Solve squared refractive index:
n_perp_1 = (-Bperp + sqrt(Bperp^2 - 4*Aperp*Cperp))/(2*Aperp);