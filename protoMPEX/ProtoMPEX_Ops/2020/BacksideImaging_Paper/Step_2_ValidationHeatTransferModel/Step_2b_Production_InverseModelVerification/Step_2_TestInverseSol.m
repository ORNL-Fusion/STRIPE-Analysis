% Test forward solution method:
% =========================================================================

% Need to add labels and text to figures and automatically generate .tiff
% and .pdf files

clear all
close all
clc

% Flags:
% =========================================================================
saveFig = 1;

% Select case:
% =========================================================================
dataset = 5:8;

% Plate thickness per case:
% =========================================================================
Lz = [0.75,1.75,3,4,0.75,1.75,3,4]*1E-3;

for jj = 1:numel(dataset)
% Load data:
% =========================================================================
matFileName = ['fwdSol_comsol_case_'   ,num2str(dataset(jj)),'.mat']; 
load(matFileName);

% Create inputField variable:
% =========================================================================
inputField.f = T;
inputField.x = xx';
inputField.y = yy';
inputField.t = tt';

% Calculate the time averaged thermal conductivity:
% =========================================================================
% Thermal conductivity of W:
k_W = @(X) 149.441 - ((45.466*1E-3)*X) + ((13.193*1E-6)*(X.^2)) - ((1.484*1E-9)*(X.^3)) + ((3.866*1E6)*(X.^(-2)));

% Mean space-integrated temperature:
T_mean = zeros(size(1:numel(tt)))';

% Mean space-integrated thermal conductivity:
k_mean = zeros(size(T_mean));
Lx = max(inputField.x)*2;
Ly = max(inputField.y)*2;
dx = inputField.x(2) - inputField.x(1);
dy = inputField.y(2) - inputField.y(1);

for ii = 1:numel(tt)
    T_mean(ii) = sum(sum(T(:,:,ii)))*dx*dy/(Lx*Ly);
    k_mean(ii) = k_W(T_mean(ii));
end

% Define simulation parameters:
% =========================================================================
% Material properties:
params.material = 'W';
params.kt       = 173;
params.kt       = mean(k_mean);
params.rho      = 19300; 
params.cp       = 134;

% Geometry:
params.Lx       = Lx;
params.Ly       = Ly;
params.Lz       = Lz(dataset(jj));

% Characteristic value for input field f:
params.f_star   = 500;
params.f_0      = 300;

% Information on input field data:
params.I        = size(inputField.f,1);
params.J        = size(inputField.f,2);
params.R        = size(inputField.f,3);
params.solType  = 'inverse';
params.num_iterations = 150;

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

% Extract data from solution:
% =========================================================================
nI = round(params.I/2);
nJ = round(params.J/2);
q2   = permute(outputField_norm.f_IJR(nI,nJ,:),[3,1,2])*params.q_star;
t_q2 = outputField_norm.t;

% Real solution:
% =========================================================================
q1   = permute(q(nI,nJ,:),[3,1,2]);
t_q1 = tt;

%% Plot data:
% =========================================================================

% Plot temperature:
figure('color','w')
hold on
hq(1) = plot(t_q1,q1*1E-6,'k');
hq(2) = plot(t_q2,q2*1E-6,'r');

% Formatting:
set(gca,'FontName','Times','FontSize',14)
xlabel('time [s]','interpreter','latex','fontSize',16)
ylabel('q [MWm$^{-2}$]','interpreter','latex','fontSize',16)
ylim([0,round(max(q1*1E-6))*1.5])
xlim([0,1])
box on
grid on
set(hq,'lineWidth',3)
hL = legend(hq,'Synthetic $q(t)$','Reconstructed $q(t)$');
set(hL,'interpreter','latex','fontSize',13)
title(['dataset ',num2str(dataset(jj))],'interpreter','latex','fontSize',15)
title(['$L_z$ = ', num2str(params.Lz*1000), ' [mm]' ])
set(gcf,'position',[360   260   452   358])

switch dataset(jj)
    case 5
        hT = text(0.03,round(max(q1*1E-6))*1.37,'(a)');
    case 6
        hT = text(0.03,round(max(q1*1E-6))*1.37,'(b)');
    case 7
        hT = text(0.03,round(max(q1*1E-6))*1.37,'(c)');
    case 8
        hT = text(0.03,round(max(q1*1E-6))*1.37,'(d)');
end

try
    set(hT,'fontSize',18,'interpreter','latex')
catch
end

% Save figure:
% =========================================================================
if saveFig
    figureName = ['Step_2_InverseSolution_case_',num2str(dataset(jj))];

    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

end
