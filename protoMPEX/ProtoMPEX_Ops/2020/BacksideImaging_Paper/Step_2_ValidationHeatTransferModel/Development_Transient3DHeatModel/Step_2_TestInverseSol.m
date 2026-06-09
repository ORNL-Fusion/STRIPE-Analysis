% Test forward solution method:
% =========================================================================

clear all
close all
clc

% Input normalized data set:
% =========================================================================
fileName = 'forwardSol_1.mat';
load(fileName);

% Convert data into physical units:
% =========================================================================
T_star = fwdSol.params.T_star;
T_0    = fwdSol.params.T_0;

inputField.f = fwdSol.outputField_norm.f_IJR*T_star + T_0;
inputField.x = fwdSol.outputField_norm.x;
inputField.y = fwdSol.outputField_norm.y;
inputField.t = fwdSol.outputField_norm.t;

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
params.Lz       = fwdSol.params.Lz;

% Characteristic value for input field f:
params.f_star   = T_star;
params.f_0      = T_0;

% Information on input field data:
params.I        = size(inputField.f,1);
params.J        = size(inputField.f,2);
params.R        = size(inputField.f,3);
params.solType  = 'inverse';

% Position where temperature and heat flux fields are located:
params.z_q     = 0;
params.z_T     = params.Lz;

% Simulation time:
params.t_start = inputField.t(1);
params.t_end   = inputField.t(end);

% Fourier-cosine modes:
params.M        = 10;
params.N        = 10;
params.L        = 1e2;

% Solve problem:
% =========================================================================
t0 = tic;
disp('Calculating solution ...')

[params,inputField_norm,outputField_norm] = Slab3dTransientThermalModel(params,inputField);

t1 = toc(t0);
disp(['Solution computed in ', num2str(t1),' seconds'])

%% Plot data:
% =========================================================================

% Extract data from solution:
nI = round(params.I/2);
nJ = round(params105J/2);
q2   = permute(outputField_norm.f_IJR(nI,nJ,:),[3,1,2])*params.q_star;
t_q2 = outputField_norm.t;

% Real solution:
q_star = fwdSol.params.q_star;
q1   = permute(fwdSol.inputField_norm.f_IJR(nI,nJ,:),[3,1,2])*q_star;
t_q1 = fwdSol.inputField_norm.t;

% Plot temperature:
figure('color','w')
hold on
hq(1) = plot(t_q1,q1,'k');
hq(2) = plot(t_q2,q2,'r');
ylim([0,20E6])
xlim([0,params.t_end])
box on
grid on
set(hq,'lineWidth',2)