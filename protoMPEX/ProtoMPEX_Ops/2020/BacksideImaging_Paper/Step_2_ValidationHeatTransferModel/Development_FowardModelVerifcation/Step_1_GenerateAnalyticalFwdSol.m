% Step 1: Generate analytical foward solution:
% =========================================================================

clear all
close all
clc

% Select analytical case:
% =========================================================================
caseSelection = 11;
computeData   = 1;
saveData      = 0;
saveFig       = 0;

% Select comsol case:
% =========================================================================
caseSelection_comsol = 6;

% Read parameter spreadsheet:
% =========================================================================
disp('Reading parameter spreadsheet...');
fileName = 'parameters_analyticalSol.xlsx';
[~,~,paramsData] = xlsread(fileName,'parameters');
ii = caseSelection + 2;
disp('Reading complete!');

if computeData
ta = tic;
% Assemble params object:
% =========================================================================
params.case               = paramsData{ii,1};
params.fileName           = [paramsData{ii,2},'case_',num2str(params.case),'.mat'];
params.material           = paramsData{ii,3};
params.rho                = paramsData{ii,4};
params.kt                 = paramsData{ii,5};
params.kt_type            = paramsData{ii,6};
params.cp                 = paramsData{ii,7};
params.q_star             = paramsData{ii,8}*1e6;
params.T_init             = paramsData{ii,9};
params.Lx                 = paramsData{ii,10}*1e-3;
params.Ly                 = paramsData{ii,11}*1e-3;
params.Lz                 = paramsData{ii,12}*1e-3;
params.t_end              = paramsData{ii,13};
params.spatialVarHeatFlux = paramsData{ii,14};
params.timeVarHeatFlux    = paramsData{ii,15};
params.mean_x             = paramsData{ii,16}*1e-3;
params.mean_y             = paramsData{ii,17}*1e-3;
params.sigma_x            = paramsData{ii,18}*1e-3;
params.sigma_y            = paramsData{ii,19}*1e-3;
params.sigma_t            = paramsData{ii,20};
params.mean_t             = paramsData{ii,21};
params.I                  = paramsData{ii,22};
params.J                  = paramsData{ii,23};
params.K                  = paramsData{ii,24};
params.R                  = paramsData{ii,25};
params.M                  = paramsData{ii,26};
params.N                  = paramsData{ii,27};
params.L                  = paramsData{ii,28};

% Initialize solid object:
% =========================================================================
solid = InitSolid(params);

% Initialize heat flux object:
% =========================================================================
heatFlux = InitHeatFlux(params,solid);

% Calculate Fourier-cosine series of heat flux profile:
% =========================================================================
disp('Computing heat flux fourier-cosine series ...');
heatFlux = CalculatefourierCosineHeatFlux(params,solid,heatFlux);
disp('Calculation complete!');
disp('')

% Initialize Toeplitz matrix:
% =========================================================================
P_RR = cell(params.M,params.N,params.K);

% Calculate Toeplitz matrix:
% =========================================================================
disp('Computing Toeplitz matrix ...');
t0 = tic;
P_RR = CalculateToeplitzMatrix(params,solid,P_RR);
t1 = toc(t0);
disp(['Calculation completed in ',num2str(t1),' seconds']);
disp('')

% Calculate foward solution:
% =========================================================================
disp('Computing foward solution ...');
t0 = tic;
solid = CalculateForwardSolution(params,solid,heatFlux,P_RR);
t1 = toc(t0);
disp(['Calculation completed in ',num2str(t1),' seconds']);
disp('')
tb = toc(ta);

% Display total calculation time:
% =========================================================================
disp(['Total calculation completed in ',num2str(tb),' seconds']);

% Save data:
% =========================================================================
if saveData
    varList = {'params','solid','heatFlux'};
    fileName = params.fileName;
    root = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\BacksideImaging_Paper\Step_2_ValidationHeatTransferModel\fwdSol_analytical_data\';
    
    disp(['file name: ',fileName]);
    disp('Save data ...');
    save([root,fileName],varList{:});
    disp(['Saving completed!']);
end

else
    % Load data:
    % =====================================================================
    fileName = [paramsData{ii,2},'case_',num2str(paramsData{ii,1})];
    root = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\BacksideImaging_Paper\Step_2_ValidationHeatTransferModel\fwdSol_analytical_data\';

    disp(['file name: ',fileName]);
    disp('loading data ...');
    load([root,fileName]);
    disp(['loading completed!']);
end

% Load COMSOL data:
% =========================================================================
% Read spreadsheet:
fileName = 'parameters_comsolSol.xlsx';
[~,~,paramsData_comsol] = xlsread(fileName,'parameters');

% Assemble fileName:
root = 'C:\Users\nfc\Documents\ProtoMPEX_Ops\2020\BacksideImaging_Paper\Step_2_ValidationHeatTransferModel\fwdSol_comsol_data\';
ii = caseSelection_comsol + 2;
params_comsol.case = paramsData_comsol{ii,1};
fileNameV{1} = ['T_back_case_',num2str(params_comsol.case),'.txt'];
fileNameV{2} = ['T_front_case_',num2str(params_comsol.case),'.txt'];
fileNameV{3} = ['T_slice_case_',num2str(params_comsol.case),'.txt'];

% load data:
T_back = load([root,fileNameV{1}]);
T_front = load([root,fileNameV{2}]);

d = readcell([root,fileNameV{3}]);
T2D_comsol = str2num(cell2mat(d(5:size(d,1))));
x_comsol   = str2num(cell2mat(d(2)));
z_comsol   = str2num(cell2mat(d(3)));

if 1
    % Demonstrate how to read time dependent surface data:
    % this will be needed for the full inverse test between COMSOL produced
    % Temperature and the fourier-inverse method to recover the applied
    % heat flux.
    
    d = readcell('T_slice_all.txt','Delimiter','');
    xx = str2num(cell2mat(d(2)));
    zz = str2num(cell2mat(d(3)));
    kk = 5;
    for jj = 1:150
        for ss = 1:30
            TT(ss,:,jj) = str2num(cell2mat(d(kk)))';
            kk = kk + 1;
        end
        kk = kk + 1;
    end
    
    figure; 
    for jj = 1:1:150
        mesh(TT(:,:,jj));
        caxis([30,1000])
        zlim([30,1000])
        view([30,60])
        drawnow
    end
end

%% Plot data:

% Back and front temperature temporal evolution:
% =========================================================================
figure('color','w')
hold on
hT(1) = plot(solid.t,permute(solid.u(round(params.I/2),round(params.J/2),1  ,:)*solid.T_star,[4,3,2,1]));
hT(2) = plot(solid.t,permute(solid.u(round(params.I/2),round(params.J/2),end,:)*solid.T_star,[4,3,2,1]));
hT(3) = plot(T_front(:,1),T_front(:,2),'r.-','lineWidth',1);
hT(4) = plot(T_back(:,1),T_back(:,2),'k.-','lineWidth',1);

legText{1} = ['Front (Analytical case:',num2str(params.case),')'];
legText{2} = ['Back (Analytical case:',num2str(params.case),')'];
legText{3} = ['Front (COMSOL case:',num2str(params_comsol.case),')'];
legText{4} = ['Back (COMSOL case:',num2str(params_comsol.case),')'];

% Formatting:
% -----------
set(gca,'fontName','Times','fontSize',11)
set(hT(1),'Color','r','lineWidth',2)
set(hT(2),'Color','k','lineWidth',2)
grid on
box on
title(['k = ', num2str(params.kt),' [Wm$^{-1}$K$^{-1}$]'],'interpreter','latex')
ylim([0,800])
xlim([0,1.5])
set(gca,'YTick',[0:50:800])
xlabel('time [s]','interpreter','latex','fontSize',12)
ylabel('$\Delta T$ [K]','interpreter','latex','fontSize',12)

% Legend:
% -------
hLeg = legend(hT(1:4),legText);
set(hLeg,'interpreter','latex','fontSize',12);

% Save figure:
% =========================================================================
figureName = ['T_Back_Front_COMSOL_',num2str(params_comsol.case),'_Analytic_',num2str(params.case)];
if saveFig
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

% 2D slice contour:
% =========================================================================
figure('color','w')

% Assemble dataset:
tt = find(solid.t >= 0.5,1);
T2D = permute(solid.u(:,round(params.J/2),:,tt),[1,3,2])'*solid.T_star + params.T_init;

% Create coordinates:
T2D_x = solid.x - params.Lx/2;
T2D_z = solid.z;

% Create plot:
subplot(2,1,1)
contourf(T2D_x*100, T2D_z*100, T2D - params.T_init,30,'LineStyle','none')

% Formatting:
set(gca,'fontName','Times','fontSize',11)
set(gca,'YTick',[0,params.Lz*100])
view([-180,90])
colormap(flipud('hot'))
axis('image')
set(gca,'PlotboxAspectRatio',[1 1/10 1])
colorbar
caxis([30,200])
title(['Analytic case: ', num2str(params.case)],'interpreter','latex')

% Create plot:
subplot(2,1,2);
contourf(x_comsol*100,(z_comsol + params.Lz/2)*100,flipud(T2D_comsol),30,'LineStyle','none')

% Formatting:
set(gca,'fontName','Times','fontSize',11)
set(gca,'YTick',[0,params.Lz*100],'YTickLabel',[0,1])
view([-180,90])
colormap(flipud('hot'))
axis('image')
% set(gca,'PlotboxAspectRatio',[1 1/10 1])
colorbar
caxis([30,200])
title(['FEM case: ',num2str(paramsData_comsol{3,1})],'interpreter','latex')

% Save figure:
% =========================================================================
figureName = ['T_slice_COMSOL_',num2str(params_comsol.case),'_Analytic_',num2str(params.case)];
if saveFig
    % PDF figure:
    exportgraphics(gcf,[figureName,'.pdf'],'Resolution',300) 

    % TIFF figure:
    exportgraphics(gcf,[figureName,'.tiff'],'Resolution',600) 
end

%% Plotting the thermal conductivity of W:

% Assemble the modified Hust-Lankford fit:
A0 = 149.441;
A1 = -45.466*1E-3;
A2 = 13.193*1E-6;
A3 = -1.484*1E-9;
A4 = 3.866*1E6;

k_W = @(x) A0 + A1*x + A2*x.^2 + A3*x.^3 + A4./(x.^2);

% Plot fit:
T1 = 300 + 280;
T2 = 300 + 0;
T = linspace(300,3000);
figure('color','w')
hold on
hk(1) = plot(T,k_W(T),'k','lineWidth',2);
plot(T1,k_W(T1),'ko')
plot(T2,k_W(T2),'ko')

% Formatting:
set(gca,'fontName','Times','fontSize',11)
box on 
grid on
title('W thermal conductivity','interpreter','latex')
xlabel('T [K]','interpreter','latex','fontSize',12)
ylabel('[W/(m K)]','interpreter','latex','fontSize',12)
ylim([0,200])


%% Functions:

% #########################################################################
% Function 1:
% #########################################################################
function solid = InitSolid(params)

% Assign temporary variables:
% ===========================
kt     = params.kt;
rho    = params.rho;
cp     = params.cp;
Ly     = params.Ly;
Lx     = params.Lx;
Lz     = params.Lz;
q_star = params.q_star;
I      = params.I;
J      = params.J;
K      = params.K;
R      = params.R;
M      = params.M;
N      = params.N;
L      = params.L;

% Thermal diffusivity:
% ====================
solid.a = kt/(rho*cp);

% Characteristic dimensions:
% ==========================
% time:
solid.t_star = (Lz^2)/solid.a;
% temperature:
solid.T_star = (q_star*Lz/kt);

% "x" coordinate:
% ===============
% "i" is the index
dx      = Lx/(I-1);
solid.x = ((1:I)-1)'*dx;
% Dimensionless space
dx_I      = dx/Lx;
solid.x_I = solid.x/Lx;

% "y" coordinate:
% ===============
% "j" is the index
dy      = Ly/(J-1);
solid.y = ((1:J)-1)'*dy;
% Dimensionless space
dy_J      = dy/Ly;
solid.y_J = solid.y/Ly;

% "z" coordinate:
% ===============
% "k" is the index
dz      = Lz/(K-1);
solid.z = ((1:K)-1)'*dz;
% Dimensionless space
dz_K      = dz/Lz;
solid.z_K = solid.z/Lz;

% "t" coordinate: 
% ===============
% "r" is the index
dt      = params.t_end/(R-1);
solid.t = ((1:R)-1)'*dt;
% Dimensionless time
dt_R      = dt/solid.t_star;
solid.t_R = solid.t/solid.t_star; 

% Fourier indices:
% =================
% "m" is the "x" Fourier series index
solid.m_M = [0:1:M-1]';
% "n" is the "y" Fourier series index
solid.n_N = [0:1:N-1]';
% "l" is the index for cos(pi*l*z)
solid.l_L = [0:1:L-1]';

% Allocate memory:
% ================
solid.u    = zeros(I,J,K,R);
end

% #########################################################################
% Function 2:
% #########################################################################
function heatFlux = InitHeatFlux(params,solid)

% Gaussian based surface heat flux:
% =================================
g = @(Y,dY,A) A*exp(-0.5*(Y/dY).^2);


% Heat flux spatial variation:
% ============================
heatFlux.q_star = params.q_star;

switch params.spatialVarHeatFlux
    case 1
        % x component:
        f  = (solid.x - params.Lx/2)/params.Lx; 
        df = params.sigma_x/params.Lx;
        heatFlux.hx = g(f,df,1);  
        
        % y component:
        f  = (solid.y - params.Ly/2)/params.Ly; 
        df = params.sigma_y/params.Ly;
        heatFlux.hy = g(f,df,1);
                
    case 2
        
    case 3
end

% Heat flux temporal variation:
% =============================
switch params.timeVarHeatFlux
    case 1                 
        f  = (solid.t - params.mean_t)/solid.t_star;
        df = params.sigma_t/solid.t_star;
        heatFlux.ht = g(f,df,1);
    case 2
end

% Assemble heat flux:
% ===================
for r = 1:params.R
      surface_xy = heatFlux.hx*heatFlux.hy';
      heatFlux.h_IJ{r} = heatFlux.ht(r)*surface_xy;
end

end

% #########################################################################
% Function 3:
% #########################################################################
function heatFlux = CalculatefourierCosineHeatFlux(params,solid,heatFlux)
% h_IJ{R}
% H_R{M,N}
% Phi_IM
% Phi_JM

[H_R,Phi_IM,Phi_JN] = Create_H_R(solid.x_I,solid.y_J,heatFlux.h_IJ,params.M,params.N);

heatFlux.H_R    = H_R;
heatFlux.Phi_IM = Phi_IM;
heatFlux.Phi_JN = Phi_JN;
end

% #########################################################################
% Function 4:
% #########################################################################
function P_RR = CalculateToeplitzMatrix(params,solid,P_RR)
% P_RR{M,N,K}
% This cell has MxNxK matrices.
% This consumes the largest memory in the calculation
% By evaluating only one "z_k" location we can reduce the size to MxN
% P_RR_new = P_RR{:,:,k};

[P_RR] = Toeplitz_MN(solid.m_M,solid.n_N,solid.l_L,params.Lx,params.Ly,params.Lz,solid.z_K,solid.t_R);
end

% #########################################################################
% Function 5:
% #########################################################################
function solid = CalculateForwardSolution(params,solid,heatFlux,P_RR)
% P_RR{M,N,K}
% H_R{M,N}
% U_MN{R,K}

solid.U_MN = Compute_U_MN(P_RR,heatFlux.H_R);

% U_MN{R,K}
% Phi_IM
% Phi_JM
% u_IJKR

solid.u = Assemble_U_MN(solid.U_MN,heatFlux.Phi_IM,heatFlux.Phi_JN);
end

% #########################################################################
% Function 6:
% #########################################################################
function [H_R,Phi_IM,Phi_JN] = Create_H_R(x_I,y_J,h_IJ,M,N)
% Takes in the surface heat flux 2D transient data and calculates the
% Fourier cosine series coeficients and eigenfunction matrices

R = size(h_IJ,2);

% Fourier cosine series representation of input signal:
% =====================================================
for r = 1:R
    [~,H_MN{r},~,~] = FourierCosine_2D_v2(x_I,y_J,h_IJ{r},M,N);
end

% Fourier-cosine eigenfunction matrices:
% ======================================
[~,~,Phi_IM,Phi_JN] = FourierCosine_2D_v2(x_I,y_J,h_IJ{R},M,N);

% Convert data into format needed:
% ================================
for m = 1:M
    for n = 1:N
        for r = 1:R
            H_R{m,n}(r,1) = H_MN{r}(m,n);
        end
    end
end

end

% #########################################################################
% Function 7:
% #########################################################################
function [P_RR] = Toeplitz_MN(m_M,n_N,l_L,Lx,Ly,Lz,z_K,t_R)
% Toeplitz_MN computes the Toeplitz matrix associated with the 3D
% transient heat conduction problem. 

M = length(m_M);
N = length(n_N);
L = length(l_L);
K = length(z_K);
R = length(t_R);

%==========================================================================
% Create the 1D Green's function for the 1D finite slab heat conduction
% problem with a surface heat flux.
% =========================================================================

% The normalized 1D finite slab heat conduction problem with an inpulse
% response heat flux (dirac delta function) is the following: 
% PDE: u_t - u_zz = 0 
% Boundary conditions: 
% u_z = -q_0*dirac(t-0) at z = 0 
% u_z = 0               at z = 1
% Initial conditions: 
% u = 0 for all "z" at t = 0;


% Calculate fundamental solution:
% ===============================
% We begin by computing the unit step response and then taking the time
% derivative to calculate the impulse response which is the 1D green's
% function

ll = l_L(2:end);
lambda_l = pi*ll;
dt_R = t_R(2) - t_R(1);

for k = 1:K
    % Eigenfunction:
    p_l = cos(lambda_l*z_K(k)); 
    cosTerm = p_l./(lambda_l.^2);
    
    % Step response:
    expTerm = exp( -(lambda_l.^2)*t_R');
    StepResponse = t_R + 2*(1-expTerm')*cosTerm;

    % Impulse response:
    ImpulseResponse = diff(StepResponse)/dt_R;
    
    % Assemble Green's function:
    G_R{k} = [ImpulseResponse;ImpulseResponse(end)];
     
end

% =========================================================================
% Create Toeplitz matrix
% =========================================================================
% It turns out that in a 3D slab heat conduction problem, the solution can
% be expressed in terms of the 1D Green's function solution "G_R"
% multiplied by some factors that depend on the "m" and "n" modes used on
% the Fourier-cosine representation of the temperature fields.

for k = 1:K
    for m = 1:M
        for n = 1:N
                % Eigenvalues:
                lambda_m = pi*m_M(m)*Lz/Lx;
                lambda_n = pi*n_N(n)*Lz/Ly;
                
                % Kernel:
                E_R{m,n} = exp(-(lambda_m.^2)*t_R).*exp(-(lambda_n.^2)*t_R);
                K_R{m,n,k} = E_R{m,n}.*G_R{k};
                
                % Toeplitz matrix:
                P_RR{m,n,k} = toeplitz(K_R{m,n,k},[K_R{m,n,k}(1);zeros(size(K_R{m,n,k})-[1,0])])*dt_R;
        end
    end
end
end

% #########################################################################
% Function 8:
% #########################################################################
function [U_MN] = Compute_U_MN(P_RR,H_R)
% Given the Fourier components of the heat flux, compute the fourier
% components of the temperature:

M = size(P_RR,1);
N = size(P_RR,2);
K = size(P_RR,3);
R = size(H_R{1,1});

% =========================================================================
% Fourier coefficients of solution
% =========================================================================
for m = 1:M
    for n = 1:N
        for k = 1:K
            U_R{m,n,k} = P_RR{m,n,k}*H_R{m,n};
        end
    end
end

% =========================================================================
% Arrange U_R in a manner that can be used for the final solution
% =========================================================================
for m = 1:M
    for n = 1:N
        for k = 1:K
            for r = 1:R
                U_MN{r,k}(m,n) = U_R{m,n,k}(r);
            end
        end
    end
end

end

% #########################################################################
% Function 9:
% #########################################################################
function [u] = Assemble_U_MN(U_MN,Phi_IM,Phi_JN)
% Given the Fourier components of the temperature and the eigenfunction,
% assemble the temperature fields in physical space:

R = size(U_MN,1);
K = size(U_MN,2);

for r = 1:R
    for k = 1:K
        f = Phi_IM*U_MN{r,k}*Phi_JN';
        u(:,:,k,r) = f;
    end
end

end