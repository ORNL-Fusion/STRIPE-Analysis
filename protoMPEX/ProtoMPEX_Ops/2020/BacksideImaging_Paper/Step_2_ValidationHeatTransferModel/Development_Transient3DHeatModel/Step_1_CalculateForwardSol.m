% Test forward solution method:
% =========================================================================

clear all
close all
clc

% Computation flags:
% =========================================================================
saveFlag = 1;
fileName = 'forwardSol_1.mat';

% Create input heat flux field:
% =========================================================================
% Geometry:
Lx       = 90E-3;
Ly       = 90E-3;

% Characteristic value for input field f:
q_star   = 20E6;
q_0      = 0;

% Information on input field data:
I        = 300;
J        = 300;
R        = 100;

% Simulation time:
t_end   = 1.5;

% Gaussian function:
g = @(x,mean_x,sigma_x) exp(-0.5*((x - mean_x)/sigma_x).^2);

% Coordinates:
t = linspace(0,1,R)'*t_end;
x = linspace(-0.5,0.5,J)'*Lx;
y = linspace(-0.5,0.5,I)'*Ly;

% Spreads:
sigma_x = 27.5E-3;
sigma_y = 27.5E-3;
sigma_t = 0.07;

% Mean locations:
mean_x = 0;
mean_y = 0;
mean_t = 0.5;

% Create profiles:
gx = g(x,mean_x,sigma_x);
gy = g(y,mean_y,sigma_y);
gt = g(t,mean_t,sigma_t);
gt2 = g(t,mean_t + 0.5,sigma_t);

% Assemble heat flux field:
for r = 1:R
    q.f(:,:,r) = (gx*gy').*(gt(r) + gt2(r))*q_star + q_0;
end
q.x = x;
q.y = y;
q.t = t;

figure('color','w')
nI = round(I/2);
nJ = round(J/2);
plot(q.t,permute(q.f(nI,nJ,:),[3,1,2]),'k')

% Define input field:
% =========================================================================
inputField = q;

% Define simulation parameters:
% =========================================================================
% Material properties:
params.material = 'W';
params.kt       = 134;
params.rho      = 19300; 
params.cp       = 134;

% Geometry:
params.Lx       = max(inputField.x)*2;
params.Ly       = max(inputField.y)*2;
params.Lz       = 3E-3;

% Characteristic value for input field f:
params.f_star   = 20E6;
params.f_0      = 0;

% Initial condition for output field:
params.T_0      = 300; 

% Information on input field data:
params.I        = size(inputField.f,1);
params.J        = size(inputField.f,2);
params.R        = size(inputField.f,3);
params.solType  = 'forward';

% Position where temperature and heat flux fields are located:
params.z_q     = 0;
params.z_T     = params.Lz;

% Simulation time:
params.t_start = inputField.t(1);
params.t_end   = inputField.t(end);

% Fourier-cosine modes:
params.M        = 50;
params.N        = 50;
params.L        = 100;

% Solve problem:
% =========================================================================
t0 = tic;
disp('Calculating solution ...')

[params,inputField_norm,outputField_norm] = Slab3dTransientThermalModel(params,inputField);

t1 = toc(t0);
disp(['Solution computed in ', num2str(t1),' seconds'])

% Save data:
% =========================================================================
if saveFlag
    % Assemble dataset:
    fwdSol.inputField_norm  = inputField_norm;
    fwdSol.outputField_norm = outputField_norm;
    fwdSol.params           = params;
   
    % Save dataset:
    disp('Saving data ...')
    save(fileName,'fwdSol')
    disp('Saving data completed !')
end

%% Plot data:
% =========================================================================

% Extract data from solution:
T_backSurface = permute(outputField_norm.f_IJR(nI,nJ,:),[3,1,2])*params.T_star;
t_backSurface = outputField_norm.t;

% Plot temperature:
figure('color','w')
plot(t_backSurface,T_backSurface,'k')
ylim([0,800])
xlim([0,params.t_end])
box on
grid on