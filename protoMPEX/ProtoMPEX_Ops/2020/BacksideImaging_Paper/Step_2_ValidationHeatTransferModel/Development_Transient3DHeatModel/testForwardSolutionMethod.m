% Test forward solution method:
% =========================================================================

clear all
close all
clc

% Define simulation parameters:
% =========================================================================
% Material properties:
params.material = 'W';
params.kt       = 134;
params.rho      = 19300; 
params.cp       = 134;

% Geometry:
params.Lx       = 90E-3;
params.Ly       = 90E-3;
params.Lz       = 10E-3;

% Characteristic value for input field f:
params.f_star   = 20E6;
params.f_0      = 0;

% Information on input field data:
params.I        = 200;
params.J        = 200;
params.R        = 100;
params.solType  = 'forward';

% Position where temperature and heat flux fields are located:
params.z_q     = 0;
params.z_T     = 0*params.Lz;

% Simulation time:
params.t_end   = 1.5;

% Fourier-cosine modes:
params.M        = 20;
params.N        = 20;
params.L        = 100;

% Input field:
% =========================================================================
g = @(x,mean_x,sigma_x) exp(-0.5*((x - mean_x)/sigma_x).^2);

t = linspace(0,1,params.R)'*params.t_end;
x = linspace(-0.5,0.5,params.J)'*params.Lx;
y = linspace(-0.5,0.5,params.I)'*params.Ly;

sigma_x = 27.5E-3;
sigma_y = 27.5E-3;
sigma_t = 0.07;

mean_x = 0;
mean_y = 0;
mean_t = 0.25;

gx = g(x,mean_x,sigma_x);
gy = g(y,mean_y,sigma_y);
gt = g(t,mean_t,sigma_t);

for r = 1:params.R
    inputField.f(:,:,r) = (gx*gy').*gt(r)*params.f_star + params.f_0;
end
inputField.x = x;
inputField.y = y;
inputField.t = t;

% Solve problem:
% =========================================================================
% Select locations for temperature fields:
z = linspace(0,params.Lz,10)';

% Calculate solution:
t0 = tic;
disp('Calculating solution ...')

for k = 1:numel(z)
    
    params.z_T = z(k);
    
    [params,inputField,outputField] = Slab3dTransientThermalModel(params,inputField);
    
    T(:,:,:,k) = outputField.f_IJR*params.T_star;
end

T = permute(T,[1,2,4,3]);
t1 = toc(t0);

disp(['Solution computed in ', num2str(t1),' seconds'])

%% Plot data:

% Plot front and back temperatures:
% =========================================================================
figure('color','w')
hold on
plot(outputField.t,permute(T(100,100,1,:),[4,3,1,2]),'r')
plot(outputField.t,permute(T(100,100,end,:),[4,3,1,2]),'k')
ylim([0,800])
xlim([0,params.t_end])
box on
grid on

% Plot x-z slice:
% =========================================================================
ii = round(params.I/2);
tt = find(outputField.t >= 0.5, 1);
T_XZ = permute(T(ii,:,:,tt),[3,2,1,4]);

figure('color','w')
contourf(outputField.x,z,flipud(T_XZ),20,'LineStyle','none')
axis image
caxis([30,200])
colormap('hot')
colorbar