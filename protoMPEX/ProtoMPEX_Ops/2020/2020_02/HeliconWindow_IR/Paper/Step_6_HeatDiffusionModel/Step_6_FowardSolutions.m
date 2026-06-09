% =========================================================================
% VARIATION SOLUTION TO THE INVERSE HEAT CONDUCTION PROBLEM, 1D FINITE
% LENGH SLAB
% 
% Here we implement the variational calculus method to minimize the least
% square expression
% =========================================================================
clear all 
% close all

figureName = 'Step_6_ForwardSolution';
VideoName = 'Step_6_InverseSolution';

% Input:
% =========================================================================
% Slab length:
L1  = 6e-3; % [mm]
L2  = 600-3; % [mm]
% Total simulation time:
T  = 1;    % [sec]
% IC:
T0 = 0;
% Characteristic heat flux:
q0 = 1e6;
% Material:
slabMaterial = 'AlN';
% Solution location:
surfaceLocation = 0;

% Slab properties
% =========================================================================
switch slabMaterial
    case 'SS'
        rho = 8000; 
        k = 15;
        cp = 500;
    case 'W'
        rho = 19300; 
        k = 173;
        cp = 134;
    case 'AlN'
        rho = 3000; 
        k = 180;
        cp = 740;
end

% Derived quantities:
% =========================================================================
a = k/(rho*cp);
T1_star = q0*L1/k;
T2_star = q0*L2/k;

% Discretization
% =========================================================================
% Spatial
% "i" is the index
Nx = 100;
dx = L1/(Nx-1);
x = ((1:Nx)-1)'*dx;
% Dimensionless space
xx   = x/L1;
dxx  = dx/L1;

% Temporal
% "j" is the index
Nt = 1e4; 
dt = T/(Nt-1);
t = ((1:Nt)-1)'*dt;
% Dimensionless time
tt  = a*t/(L1^2);
dtt = a*dt/(L1^2);

% Partial sum
% "n" is the index
Ns = 1000;

% Define applied heat flux:
% =========================================================================
% Define dimensionless time-dependent heat flux
HeatFluxType = 3;
switch HeatFluxType
    case 1
    case 2
    case 3
      t_qStart = 0.2;
      t_qEnd = 0.8;
      qs = zeros(size(t));
      rng = find(t>t_qStart & t<t_qEnd);
      qs(rng) = +1; 
    case 4
end

% Create problem condition structures:
% =========================================================================
% Problem 1:
C1.L = L1;
C1.slabMaterial = slabMaterial;
C1.rho = rho;
C1.k = k;
C1.cp = cp;
C1.qs = qs;
C1.q0 = q0;
C1.T0 = T0;
C1.T = T;
C1.Nx = Nx;
C1.Nt = Nt;
C1.Ns = Ns;
C1.surfaceLocation = surfaceLocation;
C1.T_star = T1_star;

% Solve 1D heat diffusion problem:
% =========================================================================
C1 = SolveHeatDiffusion_1D(C1);

% Infinite 1D slab heat diffusion solution:
% =========================================================================
a = k/(rho*cp);
t0 = 0.105;
t_Temp = t;
Temp = T0 + (q0/k)*(2*sqrt(a*(t_Temp-t_qStart)/pi));
Temp(rng(end):end) = [];
t_Temp(rng(end):end) = [];

figure('color','w');
hold on
plot(C1.t_Temp,C1.Temp,'k')
plot(t_Temp,Temp,'r')

%% Plot results:
% =========================================================================
close all

% Foward solution:
figure('color','w'); 
fontSize1 = 12;
fontSize2 = 12;
fontSize3 = 12;

% Applied heat flux profile:
subplot(1,2,1)
hold on
if 0
    for i = 1:1:length(q)
        plot(t,q{i},'m:')
    end
end
plot(t,q0*qs*1e-6,'k.-','linewidth',2)
ylim([0,1.2])
xlim([0,1])
grid on
box on
set(gcf,'color','w')
xlabel('time [sec]','interpreter','latex','fontsize',fontSize1)
ylabel('Heat flux [MWm$^{-2}$]','interpreter','latex','fontsize',fontSize1)
set(gca,'fontName','Times','Fontsize',fontSize2)
axis square
% Annotation:
text(0.45,1.08,'$q(t)$','interpreter','Latex','FontSize',fontSize3)
ax(1) = gca;

% Surface temperature profile:
subplot(1,2,2)
hold on
if 0
    for i = 1:10:length(u)
        plot(C1.t_Temp,C1.Temp,'m:')
    end
end
hM (1) = plot(C1.t_Temp,C1.Temp,'k','linewidth',2)
hM (2) = plot(t_Temp,Temp,'r','linewidth',2)
xlabel('time [sec]','interpreter','latex','fontsize',fontSize1)
ylabel('$\Delta$T [K]','interpreter','latex','fontsize',fontSize1)
set(gca,'fontName','Times','Fontsize',fontSize2)
xlim([0,1])
box on
grid on
title(['Material: ',slabMaterial],'interpreter','latex','fontsize',fontSize3)
axis square
% Annotation:
text(0.08,52,'$T(x = 0,t) - T_0$','interpreter','Latex','FontSize',fontSize3)
hLeg = legend(hM,'Finite','Infinite');
hLeg.Interpreter = 'Latex';
hLeg.FontSize = 12;
hLeg.Location = 'southeast';
ax(2) = gca;

% Final formatting:
% =========================================================================
set(ax,'FontName','times','FontSize',11)
text(ax(1),0.02,1.3,'(a)','Interpreter','latex','FontSize',13)
text(ax(2),0.02,65,'(b)','Interpreter','latex','FontSize',13)

% Save figure:
% =========================================================================
saveas(gcf,figureName,'tiffn')

 
 %% Functions:
 function C = SolveHeatDiffusion_1D(C)
 % Solves for the 1D heat diffusion problem associated with the helicon
 % window
 % C: structure that defines the problem to solve
 % P: Toeplitz matrix 
 % u: dimensionless temperature profile
 
 % Input data representing the problem conditions:
 % ========================================================================
 % Geometry:
 L = C.L;
 
 % Material properties:
 slabMaterial = C.slabMaterial;
 rho = C.rho;
 k = C.k;
 cp = C.cp;
 
 % Boundary and initial conditions:
 qs = C.qs;
 q0 = C.q0;
 T0 = C.T0;
 
 % Solution domain:
 T = C.T;
 Nx = C.Nx;
 Nt = C.Nt; 
 Ns = C.Ns;
 surfaceLocation = C.surfaceLocation;
 
 % Derived quantities:
% =========================================================================
a = k/(rho*cp);
t_star = L*L/a;
T_star = q0*L/k;

% Discretization
% =========================================================================
% Spatial
% "i" is the index
dx = L/(Nx-1);
x = ((1:Nx)-1)'*dx;
% Dimensionless space
xx   = x/L;
dxx  = dx/L;

% Temporal
% "j" is the index
dt = T/(Nt-1);
t = ((1:Nt)-1)'*dt;
% Dimensionless time
tt  = t/t_star;
dtt = dt/t_star;

% Create the convolution operator:
% =========================================================================
SurfaceLocation = 0;
switch SurfaceLocation 
    case 0 % Front Surface 
        xx0 = 0.005;
        % To capture the details of the front surface we need to use 1e4
        % partial sum terms and use xx0 = 0.01 to 0.05
    case 1 % Back Surface
        xx0 = 1;
end
 for s = 1:length(tt)
       G(s) = GreensFunction(xx0,tt(s),Ns);
 end
P = toeplitz(G,zeros(size(G)));
 
% Compute temperature profile (dimensionless):
u = P*qs*dtt; 

% To convert to physical units:
Temp = u*T_star + T0;
x_Temp = xx*L;
t_Temp = tt*t_star;

% Output data:
% =========================================================================
% Normalized data:
C.u = u;
C.tt = tt;
C.xx = xx;

% Data in physical units:
C.Temp = Temp;
C.x_Temp = x_Temp;
C.t_Temp = t_Temp;

% Convolution operator:
C.P = P;

% Derived quantities:
C.a = a;
C.T_star = T_star;
C.t_star = t_star;
 end
 
 function [G] = GreensFunction(x,t,Ns)
%==========================================================================
% 1D-TRANSIENT HEAT CONDUCTION SOLUTION WITH IMPULSE SURFACE HEAT FLUX
% Independent variables: x and t
%==========================================================================
% GreensFunction is the solution of the 1D transient heat diffusion equation subject
% to an IMPULSE surface heat flux at x = 0. It can be viewed as the Green's
% function solution to the problem:

% u_t = u_xx
%
% u(x,t=0) = 0
% ux @ (x = 0) = -delta(t)
% ux @ (x = 1) = 0

% Where u is the dimensionless temperature u = (T-T0)/(q0*L/k), t is the
% dimensionless time (alpha*t/L^2) and x is the dimensionless time (x/L)
% subscripts t and xx represent first temporal and second spatial derivatives.

% x: axial location normalized to L, the axial length of slab, R^1
% t: dimensionless time value, R^1
% Ns: total number of terms in the infinite sum

if (t ==0 && x == 1)
    G  = 0;
else
    C = exp(-(([1:Ns]'*pi).^2)*t);
    E = cos([1:Ns]'*pi*x);
    G  = 1 + 2*C'*E;
end

end