% =========================================================================
% DIRECT COMPARISON BETWEEN FINITE AND INFINITE LENGTH SLAB WITH PRESCRIBED
% SURFACE HEAT FLUX
% 
% For time less than ts/3, where ts = L*L/alpha, both finite and infinite
% models are the same. under these conditions the thickness of the plate is
% not important. For a constant surface heat flux, after t>ts/3, the main
% characteristic of the finite slab model is that T scales linearly with t
% while in the infinite model T scales linearly with t^0.5

% How do we assess the importance of 3D effects?
% I think is related to the spatial fourier content of the surface
% temperature distribution. if it can be well described by a single cosine
% function, then a the 1D model is adequate
% how to we express the 1D model in the context of the full 3D solution?
% From my notes I can see that the 1D model is actually part of the 3D
% solution and that higher order modes are added to reconstruct the 2D
% nature of the surface temperature
% =========================================================================

clear all
close all

% Slab properties

Material = 'ALN';
L  = (0.25*25.4)*1e-3;
t = linspace(0,0.7,1e2);
q0 = 0.2e6; % Characteristic heat flux

switch Material
    case 'SS'
        rho = 8000; 
        k   = 15;
        cp  = 500;
    case 'W'
        rho = 19300; 
        k   = 173;
        cp  = 134;
    case 'ALN'
        rho = 3300; 
        k   = 180 ;
        cp  = 740 ;
end

a = k/(rho*cp);
T0 = 0;

x = linspace(0,L);

Fo = a*t/(L^2);
X  = x/L;

N = 500;
n = 1:N;
f = @(i,x,fo) (1./(i.^2)).*cos(i*pi*x).*(1-exp(-((i*pi).^2).*fo));

% Theta is the normalized surface temperature
Theta = zeros(size(X));
for nt = 1:length(Fo)
    for nx = 1:length(X)
        S = sum(f(n,X(nx),Fo(nt)));
        Theta(nt,nx) = Fo(nt) + (2/(pi^2))*S;
        
    end
end
T = T0 + (q0*L/k)*Theta;

% Infinite slab solution
eta = @(x,t) x./sqrt(4*a*t);
T1 = @(x,t) T0 + (q0/k)*(sqrt(4*a*t/pi).*exp(-eta(x,t).^2) - x.*erfc(eta(x,t)));

for nt = 1:length(Fo)
    for nx = 1:length(X)
        T_InfiniteSlab(nt,nx) = T1(x(nx),t(nt));
    end
end
T_InfiniteSlab(1,1) = 0;
%         T_InfiniteSlab = T0 + (q0/k)*sqrt(4*a/pi)*sqrt(t);

%% Animated figure
figure
set(gcf,'color','w')
box on
rng = find(t<5);
for nt = 1:1:length(t(rng))
    plot(x*1e3,T(nt,:),'k','LineWidth',2)
    hold on
    plot(x*1e3,T_InfiniteSlab(nt,:),'g','LineWidth',2)

%     ylim([0,max(max(T))*1.2])
    ylim([0,2*max(max(T))])
    xlim([0,L*1e3])
    ylabel('{\Delta}T')
    xlabel('[mm]')
    title(['t: ',num2str(t(nt)*1e3),' [ms]',', Absolute T(z) vs time'])
    drawnow
    hold off
end

%%
% Front surface temperature vs time
figure; 
hold on
hS(1) = plot(Fo,T(:,1),'k.-')
hS(2) = plot(Fo,T_InfiniteSlab(:,1),'r.-')
title('Time evolution of T surface')
legend(hS,'FiniteSlab','InfiniteSlab')
box on
set(gcf,'color','w')
grid on
ylabel('T [C]')
xlabel('t/{t^*}')

% Checking linearity of (T-Tmin)^2
Tmin{1} = min(T(:,1)             );
Tmin{2} = min(T_InfiniteSlab(:,1));
t0 = t(1);

T2{1} = (             T(:,1)-Tmin{1}).^2;
T2{2} = (T_InfiniteSlab(:,1)-Tmin{2}).^2;


figure
hold on
hS(1) = plot(t-t0,T2{1},'k.-');
hS(2) = plot(t-t0,T2{2},'r.-');
title('Time evolution of (T-Tmin)^2 surface')
legend(hS,'FiniteSlab','InfiniteSlab')
xlabel('time [s]')
box on
grid on
set(gcf,'color','w')

% Reconstructing heat flux from synthetic data

dT2dt{1} = diff(T2{1}')./diff(t-t0);
dT2dt{2} = diff(T2{2}')./diff(t-t0);
ts = L*L/a;

figure; 
hold on
hS(1) = plot(t(1:end-1)-t0,sqrt((pi/4)*k*rho*cp*dT2dt{1})*1e-3,'k.-');
hS(2) = plot(t(1:end-1)-t0,sqrt((pi/4)*k*rho*cp*dT2dt{2})*1e-3,'r.-');
ylim([0,300])
ylabel('kWm^{-2}')
title('Heat flux reconstruction from (T-Tmin)^2')
legend(hS,'FiniteSlab','InfiniteSlab')
xlabel('time [s]')
box on
grid on
set(gcf,'color','w')